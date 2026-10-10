"""A single prescribed longitudinal crack surface; all other faces stay welded.

No physical fracture coefficient is selected here. Binding requires an explicit
conditional contract. Interface evaluations never mutate accepted history.
"""
from dataclasses import dataclass
from copy import copy
import ctypes,json
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.sparse import coo_matrix,kron,eye
from periodic_sector import PeriodicSector,rotation_y,ROOT
from volume_coupon import shapes,ORGANIC
from native_cohesion_r2 import NativeCohesion


@dataclass(frozen=True)
class ConditionalFracture:
    peak_Pa: float
    Gc_J_m2: float
    penalty_Pa_per_m: float
    rationale: str
    evidence_status: str = 'conditional scenario, not measured cohesive cork parameters'

    def validate(self):
        values=np.array([self.peak_Pa,self.Gc_J_m2,self.penalty_Pa_per_m])
        if not np.isfinite(values).all() or np.any(values<=0):raise ValueError('Positive finite explicit scenario values required')
        if not self.rationale.strip():raise ValueError('A written scenario rationale is required')
        if self.delta_f<=self.delta_0:raise ValueError('The bilinear law requires 2 Gc/T0 > T0/Kpen')
        return self
    @property
    def delta_0(self):return self.peak_Pa/self.penalty_Pa_per_m
    @property
    def delta_f(self):return 2*self.Gc_J_m2/self.peak_Pa
    def report(self):
        self.validate()
        return dict(peak_Pa=self.peak_Pa,Gc_J_m2=self.Gc_J_m2,penalty_Pa_per_m=self.penalty_Pa_per_m,
            delta_0_m=self.delta_0,delta_f_m=self.delta_f,bilinear_integrated_work_J_m2=.5*self.peak_Pa*self.delta_f,
            rationale=self.rationale,evidence_status=self.evidence_status,
            crack_system=dict(normal='tangential',advance='radial',front='axial'),
            limitations='A prescribed repeated crack path. No identified mixed-mode law, friction, virgin-cork calibration or crack-history qualification.')


class SingleSeam:
    def __init__(self,control=None,face_order=4,allow_underintegrated_control=False):
        if face_order<3 and not allow_underintegrated_control:
            raise ValueError('P2 interface needs at least 3x3 quadrature; order 2 is reserved for the explicit failed control')
        self.control=PeriodicSector() if control is None else control
        self.face_order=face_order
        b=self.control.base;self.base=copy(b);self.config=b.config
        center=np.isclose(b.X[:,0],0.,atol=1e-12,rtol=0);ids=np.flatnonzero(center)
        duplicate=np.full(len(b.X),-1,int);duplicate[ids]=len(b.X)+np.arange(len(ids))
        self.base.X=np.vstack((b.X,b.X[ids]));self.base.map=b.map.copy()
        positive_cell=b.Xcell.mean(axis=1)[:,0]>0
        for cell in np.flatnonzero(positive_cell):
            old=self.base.map[cell];on=center[old];self.base.map[cell,on]=duplicate[old[on]]
        self.cell_map=self.base.map;self.duplicated_nodes=ids
        tags=np.zeros(len(self.base.X),int);tags[len(b.X):]=1
        X=self.base.X;r=np.linalg.norm(X[:,[0,2]],axis=1);sx=np.arctan2(X[:,0],X[:,2])*self.config['outer_radius_m'];y=X[:,1]
        right=np.isclose(sx,self.config['width_m']/2,atol=1e-11,rtol=0);upper=np.isclose(y,self.config['height_m']/2,atol=1e-11,rtol=0)
        keys=np.column_stack((r,np.where(right,-self.config['width_m']/2,sx),np.where(upper,-self.config['height_m']/2,y),tags))
        _,first,self.groups=np.unique(np.round(keys,12),axis=0,return_index=True,return_inverse=True)
        self.X=X[first].copy();self.rotations=np.broadcast_to(np.eye(3),(len(X),3,3)).copy();self.rotations[right]=rotation_y(self.control.alpha)
        rows=np.broadcast_to(np.arange(len(X))[:,None,None]*3+np.arange(3)[None,:,None],self.rotations.shape)
        cols=np.broadcast_to(self.groups[:,None,None]*3+np.arange(3)[None,None,:],self.rotations.shape)
        self.T=coo_matrix((self.rotations.ravel(),(rows.ravel(),cols.ravel())),shape=(len(X)*3,len(self.X)*3)).tocsr()
        self.base.fixed_nodes=np.unique(self.cell_map[:b.T,:6]);self.base.free_nodes=np.setdiff1d(np.arange(len(X)),self.base.fixed_nodes)
        self.fixed_nodes=np.unique(self.groups[self.base.fixed_nodes]);self.free_nodes=np.setdiff1d(np.arange(len(self.X)),self.fixed_nodes)
        M=np.array([cell.mass_scalar for cell in b.cells]);mp=self.cell_map
        self.base.mass_matrix=coo_matrix((M.ravel(),(np.broadcast_to(mp[:,:,None],M.shape).ravel(),np.broadcast_to(mp[:,None,:],M.shape).ravel())),shape=(len(X),len(X))).tocsr()
        self.mass_matrix=(self.T.T@kron(self.base.mass_matrix,eye(3),format='csr')@self.T).tocsr()
        self.mass=float(self.base.mass_matrix.sum());self.base.mass=self.mass
        self.geometry,self.faces=self.build_interface_geometry()
        self.lib=ctypes.CDLL(str(ORGANIC/'causal-bark/volume-core/source/libcohesion_incremental.so'))
        P=ctypes.POINTER(ctypes.c_double);self.function=self.lib.cohesive_incremental
        self.function.argtypes=[ctypes.c_int,ctypes.c_int,ctypes.POINTER(ctypes.c_int)]+[P]*13;self.function.restype=ctypes.c_int

    def expand(self,u):return np.asarray(self.T@u.ravel()).reshape(self.base.X.shape)
    def bulk(self,u,**kwargs):
        E,g,r=self.base.evaluate(self.expand(u),**kwargs)
        return E,np.asarray(self.T.T@g.ravel()).reshape(self.X.shape),r
    def kinetic(self,v):return .5*float(v.ravel()@(self.mass_matrix@v.ravel()))
    def build_interface_geometry(self):
        b=self.base;edge_map={};faces=[]
        for k,tri in enumerate(b.tri):
            for ia,ib in ((0,1),(1,2),(2,0)):
                if not np.all(abs(b.points[tri[[ia,ib]],0])<1e-12):continue
                key=tuple(sorted((int(tri[ia]),int(tri[ib]))))
                if key not in edge_map:edge_map[key]=(k,ia,ib)
                else:
                    j,ja,jb=edge_map.pop(key)
                    if b.xy[j,:,0].mean()<0:neg=(j,ja,jb);pos=k
                    else:neg=(k,ia,ib);pos=j
                    ni,na,nb=neg;pa=list(b.tri[pos]).index(int(b.tri[ni,na]));pb=list(b.tri[pos]).index(int(b.tri[ni,nb]))
                    faces.append((ni,na,nb,pos,pa,pb))
        assert not edge_map
        data={k:[] for k in ('pairs','J','T1','T2','area','normal_sign','tref')}
        gauss,weights=leggauss(self.face_order)
        for layer in range(3):
            for i,ia,ib,j,ja,jb in faces:
                ci=layer*b.T+i;cj=layer*b.T+j;X=np.vstack((b.Xcell[ci],b.Xcell[cj]))
                e1=np.zeros(3);e2=np.zeros(3);e1[ia]=-1;e1[ib]=1;e2[ja]=-1;e2[jb]=1
                d1=np.r_[e1[1:],0.];d2=np.r_[e2[1:],0.]
                for t,wt in zip((gauss+1)/2,weights/2):
                    p1=np.zeros(3);p2=np.zeros(3);p1[ia]=1-t;p1[ib]=t;p2[ja]=1-t;p2[jb]=t
                    for z,wz in zip(gauss,weights):
                        N1,D1=shapes(p1[1],p1[2],z);N2,D2=shapes(p2[1],p2[2],z)
                        J=np.r_[-N1,N2];A=np.r_[.5*(D1@d1),.5*(D2@d2)];B=np.r_[.5*D1[:,2],.5*D2[:,2]]
                        J-=J.mean();A-=A.mean();B-=B.mean();xc=X-X[:1];t1=A@xc;t2=B@xc;cross=np.cross(t1,t2)
                        sign=1. if cross@(X[18:].mean(axis=0)-X[:18].mean(axis=0))>0 else -1.
                        row=dict(pairs=(ci,cj),J=J,T1=A,T2=B,area=wt*wz*np.linalg.norm(cross),normal_sign=sign,tref=np.r_[t1,t2])
                        for key,value in row.items():data[key].append(value)
        return {key:np.array(value,dtype=np.int32 if key=='pairs' else float) for key,value in data.items()},faces

    def bind(self,contract):
        contract.validate();n=len(self.geometry['area']);ops={k:v for k,v in self.geometry.items() if k!='tref'}
        ops.update(kt=np.full(n,contract.penalty_Pa_per_m),kn=np.full(n,contract.penalty_Pa_per_m),strength=np.full(n,contract.peak_Pa),Gc=np.full(n,contract.Gc_J_m2))
        return NativeCohesion(len(self.base.cells),ops)

    def interface(self,u,bonds,damage=None):
        # Caller-owned trial history: accepted bonds history is never changed.
        bonds._check_state();damage=bonds.damage.copy() if damage is None else np.array(damage,dtype=float,copy=True,order='C')
        if damage.shape!=bonds.damage.shape or np.any(damage<0) or np.any(damage>1):raise ValueError('Invalid trial damage')
        local=np.ascontiguousarray(self.expand(u)[self.cell_map]);o=bonds.ops;n=bonds.N
        E=np.zeros(1);g=np.zeros_like(local);opening=np.zeros(n);P=ctypes.POINTER(ctypes.c_double)
        ptr=lambda a:a.ctypes.data_as(P)
        rc=self.function(len(self.base.cells),n,o['pairs'].ctypes.data_as(ctypes.POINTER(ctypes.c_int)),ptr(local),ptr(self.geometry['tref']),ptr(o['J']),ptr(o['T1']),ptr(o['T2']),ptr(o['area']),ptr(o['kt']),ptr(o['kn']),ptr(damage),ptr(o['normal_sign']),ptr(E),ptr(g),ptr(opening))
        if rc:raise ValueError(f'Incremental cohesive status {rc}')
        assembled=np.zeros_like(self.base.X);np.add.at(assembled,self.cell_map,g)
        return float(E[0]),np.asarray(self.T.T@assembled.ravel()).reshape(self.X.shape),opening

    @staticmethod
    def propose_history(bonds,accepted_maximum,opening):
        maximum=np.maximum(np.array(accepted_maximum,copy=True),opening)
        return maximum,bonds._law(maximum)
