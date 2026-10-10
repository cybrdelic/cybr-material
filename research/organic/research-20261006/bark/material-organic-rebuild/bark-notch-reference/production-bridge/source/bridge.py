"""Constrain real production prisms; solve/evaluate their native energy."""
from pathlib import Path
from types import SimpleNamespace
import ctypes
import sys
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT/'frozen-production'
sys.path.insert(0, str(FROZEN))
from quadratic_prism import QuadraticPrism, flat_reference, shapes
from native_volume import NativePrisms
from native_cohesion_r2 import NativeCohesion

DP = ctypes.POINTER(ctypes.c_double)
def ptr(a):
    return a.ctypes.data_as(DP)


class Bridge:
    def __init__(self, n, cfg):
        self.cfg, self.n = cfg, n
        self.width, self.height, self.a, self.t = .04, .04, .02, .001
        self.E, self.nu = cfg['material']['E_Pa'], cfg['material']['nu']
        self.h = self.width/n
        xx, yy = np.meshgrid(np.linspace(0, self.width, n+1), np.linspace(0, self.height, n+1))
        self.vertices = np.column_stack([xx.ravel(), yy.ravel()])
        ids = np.arange((n+1)**2).reshape(n+1, n+1)
        lo, lr, hi, ul = ids[:-1,:-1].ravel(), ids[:-1,1:].ravel(), ids[1:,1:].ravel(), ids[1:,:-1].ravel()
        self.tri = np.stack([np.column_stack([lo, lr, hi]), np.column_stack([lo, hi, ul])], axis=1).reshape(-1, 3)
        self.cells = [QuadraticPrism(flat_reference(self.vertices[tri], self.t), self.E, self.nu,
                       cfg['material']['density_kg_m3'], nq=cfg['bulk_quadrature_nq']) for tri in self.tri]
        self.N = len(self.cells)
        self.X = np.stack([c.X for c in self.cells])
        self.xy = self.X[:, :6, :2]
        self.bulk = NativePrisms(self.cells)
        self.H0 = np.zeros((self.N, len(self.cells[0].weights), 3, 3))
        self.ops_unit, self.facet_point_rows = self.make_interfaces()
        pair_X = np.stack([np.vstack([self.X[i], self.X[j]]) for i,j in self.ops_unit['pairs']])
        pair_X -= pair_X[:, :1]
        self.tref = np.ascontiguousarray(np.column_stack([
            np.einsum('pi,pia->pa', self.ops_unit['T1'], pair_X),
            np.einsum('pi,pia->pa', self.ops_unit['T2'], pair_X)]))
        self.vlib = ctypes.CDLL(str(FROZEN/'libvolume_incremental.so'))
        self.vfun = self.vlib.volume_incremental
        self.vfun.argtypes = [ctypes.c_int, ctypes.c_int]+[DP]*10
        self.vfun.restype = ctypes.c_int
        self.clib = ctypes.CDLL(str(FROZEN/'libcohesion_incremental.so'))
        self.cfun = self.clib.cohesive_incremental
        self.cfun.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_int)]+[DP]*13
        self.cfun.restype = ctypes.c_int
        self.Kbulk_local = self.bulk_tangent()
        self.Kinterface_local = self.interface_tangent()

    def make_interfaces(self):
        rows = {k: [] for k in ('pairs','J','T1','T2','area','kt','kn','normal_sign','strength','Gc')}
        edges, facets = {}, []
        for c, tri in enumerate(self.tri):
            for ia,ib in ((0,1),(1,2),(2,0)):
                key = tuple(sorted((int(tri[ia]),int(tri[ib]))))
                if key in edges:
                    d,ja,jb = edges.pop(key)
                    ca = list(tri).index(int(self.tri[d,ja])); cb = list(tri).index(int(self.tri[d,jb]))
                    facets.append((d,ja,jb,c,ca,cb))
                else: edges[key] = (c,ia,ib)
        g,w = leggauss(self.cfg['facet_quadrature_order'])
        point_groups = []
        for i,ia,ib,j,ja,jb in facets:
            group = []
            edge = self.vertices[self.tri[i,ib]]-self.vertices[self.tri[i,ia]]
            normal = np.array([-edge[1],edge[0]])/np.linalg.norm(edge)
            hnormal = abs(normal @ (self.vertices[self.tri[j]].mean(0)-self.vertices[self.tri[i]].mean(0)))
            K = self.E/hnormal
            e1,e2 = np.zeros(3),np.zeros(3)
            e1[ia],e1[ib],e2[ja],e2[jb] = -1,1,-1,1
            d1,d2 = np.r_[e1[1:],0],np.r_[e2[1:],0]
            X = np.vstack([self.X[i],self.X[j]])
            for s,ws in zip((g+1)/2,w/2):
                p1,p2 = np.zeros(3),np.zeros(3)
                p1[ia],p1[ib],p2[ja],p2[jb] = 1-s,s,1-s,s
                for z,wz in zip(g,w):
                    N1,D1 = shapes(p1[1],p1[2],z); N2,D2 = shapes(p2[1],p2[2],z)
                    J = np.r_[-N1,N2]
                    A = np.r_[.5*(D1@d1),.5*(D2@d2)]
                    B = np.r_[.5*D1[:,2],.5*D2[:,2]]
                    J -= J.mean(); A -= A.mean(); B -= B.mean()
                    cross = np.cross(A@X,B@X)
                    sign = 1. if cross@(X[18:].mean(0)-X[:18].mean(0))>0 else -1.
                    group.append(len(rows['area']))
                    vals = ((i,j),J,A,B,ws*wz*np.linalg.norm(cross),K,K,sign,180000.,10.)
                    for key,value in zip(rows,vals): rows[key].append(value)
            point_groups.append(group)
        return {k:np.asarray(v,dtype=np.int32 if k=='pairs' else float) for k,v in rows.items()}, point_groups

    @staticmethod
    def assemble(blocks, ids, size):
        d = ids.shape[1]
        return coo_matrix((np.asarray(blocks).ravel(),
                          (np.repeat(ids,d,axis=1).ravel(),np.tile(ids,(1,d)).ravel())),
                         shape=(size,size)).tocsr()

    def bulk_tangent(self):
        D = np.diag([self.E,self.E,self.E/2])
        blocks = []
        for cell in self.cells:
            G = cell.grad.reshape(-1,3,6,3).sum(axis=1)
            B = np.zeros((len(G),3,12))
            B[:,0,0::2] = G[:,:,0]; B[:,1,1::2] = G[:,:,1]
            B[:,2,0::2] = G[:,:,1]; B[:,2,1::2] = G[:,:,0]
            blocks.append(np.einsum('qai,ab,qbj,q->ij',B,D,B,cell.weights))
        return self.assemble(blocks,np.arange(self.N*12).reshape(self.N,12),self.N*12)

    def interface_tangent(self):
        blocks, dofs = [], []
        o = self.ops_unit
        for group in self.facet_point_rows:
            i,j = o['pairs'][group[0]]
            jc = o['J'][group].reshape(-1,2,3,6).sum(axis=2).reshape(-1,12)
            scalar = np.einsum('qi,qj,q->ij',jc,jc,o['area'][group]*o['kt'][group])
            blocks.append(np.kron(scalar,np.eye(2)))
            dofs.append(np.r_[np.arange(i*12,(i+1)*12),np.arange(j*12,(j+1)*12)])
        return self.assemble(blocks,np.array(dofs),self.N*12)

    def prepare(self, eta):
        # eta=None identifies coincident nodes exactly. Same production energies.
        self.eta = eta
        if eta is None:
            grid = np.rint(self.xy/(self.h/2)).astype(int)
            self.node_ids = grid[:,:,1]*(2*self.n+1)+grid[:,:,0]
            nnodes = (2*self.n+1)**2
        else:
            self.node_ids = np.arange(self.N*6).reshape(self.N,6)
            nnodes = self.N*6
        self.dof_ids = np.stack([2*self.node_ids,2*self.node_ids+1],axis=-1).reshape(self.N,12)
        self.map = self.dof_ids.ravel()
        self.size = 2*nnodes
        def condense(K):
            c=K.tocoo()
            return coo_matrix((c.data,(self.map[c.row],self.map[c.col])),shape=(self.size,self.size)).tocsr()
        self.Kb = condense(self.Kbulk_local)
        self.Kc = condense(self.Kinterface_local)/(eta if eta else .0005)
        # Identified traces have mathematically zero jump. Avoid cancellation
        # from assembling a null penalty after exact node identification.
        self.K = self.Kb if eta is None else self.Kb+self.Kc
        tol = self.h*1e-8
        top_nodes = np.unique(self.node_ids[np.abs(self.xy[:,:,1]-self.height)<tol])
        bottom_nodes = np.unique(self.node_ids[(abs(self.xy[:,:,1])<tol)&(self.xy[:,:,0]>=self.a-tol)])
        anchor = self.node_ids[(abs(self.xy[:,:,0]-self.width)<tol)&(abs(self.xy[:,:,1]-self.height)<tol)][0]
        self.top_dofs = 2*top_nodes+1; self.bottom_dofs = 2*bottom_nodes+1; self.anchor = 2*anchor
        self.fixed = np.unique(np.r_[self.top_dofs,self.bottom_dofs,self.anchor])
        mask=np.ones(self.size,dtype=bool);mask[self.fixed]=False
        self.free=np.flatnonzero(mask)
        self.lu=splu(self.K[self.free][:,self.free].tocsc())
        # The validated production constructor makes its own owned copies.
        # Avoid a redundant third full copy of immutable interpolation arrays.
        op=dict(self.ops_unit)
        op['kt']=op['kt']/(eta if eta else .0005)
        op['kn']=op['kn']/(eta if eta else .0005)
        self.bonds=NativeCohesion(self.N,op)

    def lift(self,u):
        local=u[self.dof_ids].reshape(self.N,6,2)
        out=np.zeros((self.N,3,6,3));out[:,:,:,:2]=local[:,None]
        return np.ascontiguousarray(out.reshape(self.N,18,3))

    def restrict(self,g):
        local=g.reshape(self.N,3,6,3).sum(axis=1)[:,:,:2]
        return np.bincount(self.map,weights=local.ravel(),minlength=self.size)

    def evaluate(self,u):
        displacement=self.lift(u)
        gb=np.zeros_like(displacement); Eb=np.zeros(self.N); V=np.zeros(self.N); minJ=np.zeros(1)
        rc=self.vfun(self.N,self.bulk.weights.shape[1],ptr(displacement),ptr(self.H0),ptr(self.bulk.grad),ptr(self.bulk.weights),ptr(self.bulk.mu),ptr(self.bulk.lam),ptr(gb),ptr(Eb),ptr(V),ptr(minJ))
        if rc:raise RuntimeError('Production bulk status '+str(rc))
        gc=np.zeros_like(displacement); Ec=np.zeros(1); opening=np.zeros(self.bonds.N); o=self.bonds.ops
        self.bonds._check_state()
        rc=self.cfun(self.N,self.bonds.N,o['pairs'].ctypes.data_as(ctypes.POINTER(ctypes.c_int)),ptr(displacement),ptr(self.tref),ptr(o['J']),ptr(o['T1']),ptr(o['T2']),ptr(o['area']),ptr(o['kt']),ptr(o['kn']),ptr(self.bonds.damage),ptr(o['normal_sign']),ptr(Ec),ptr(gc),ptr(opening))
        if rc:raise RuntimeError('Production cohesive status '+str(rc))
        if not all(np.isfinite(a).all() for a in (gb,Eb,V,minJ,gc,Ec,opening)):raise FloatingPointError('Nonfinite production result')
        return {'energy_half_J':float(Eb.sum()+Ec[0]),'bulk_energy_half_J':float(Eb.sum()),
                'interface_energy_half_J':float(Ec[0]),'g':self.restrict(gb+gc),
                'gb':self.restrict(gb),'gc':self.restrict(gc),
                'max_opening_over_initiation':float(np.max(opening/self.bonds.delta0)),
                'max_opening_m':float(opening.max()),'min_det_F':float(minJ[0]),
                'max_omitted_z_force_N':float(np.abs((gb+gc)[:,:,2]).max()),
                'full_nodal_bulk_gradient':gb,'full_nodal_interface_gradient':gc}

    def solve(self,delta):
        u=np.zeros(self.size);u[self.top_dofs]=delta
        u[self.free]=self.lu.solve(-(self.K@u)[self.free])
        linear=u.copy();trace=[]
        for it in range(20):
            ev=self.evaluate(u);P=float(ev['g'][self.top_dofs].sum())
            residual=float(np.linalg.norm(ev['g'][self.free])/abs(P))
            trace.append(residual)
            if residual < 1e-10:break
            u[self.free]+=self.lu.solve(-ev['g'][self.free])
        else:raise RuntimeError('Native equilibrium failed: '+str(trace))
        ev.update({'u':u,'linear_u':linear,'P_N':P,'Delta_m':2*delta,
                   'compliance_m_N':2*delta/P,'energy_full_J':2*ev['energy_half_J'],
                   'free_residual_relative':residual,'iteration_residuals':trace,
                   'bottom_half_reaction_N':float(ev['g'][self.bottom_dofs].sum()),
                   'horizontal_anchor_reaction_N':float(ev['g'][self.anchor]),
                   'reaction_balance_relative':abs(float(ev['g'][self.bottom_dofs].sum())+P)/abs(P)})
        return ev


def serial(ev):
    return {k:v for k,v in ev.items() if not isinstance(v,np.ndarray)}
