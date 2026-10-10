"""Conditional single-path energy/force with a consistent irreversible merit."""
from dataclasses import dataclass
import json,time
import numpy as np
from scipy.sparse import coo_matrix,kron,eye
from scipy.sparse.linalg import splu
from periodic_sector import PeriodicSector,ROOT
from single_seam import SingleSeam,ConditionalFracture


@dataclass(frozen=True)
class History:
    maximum: np.ndarray
    damage: np.ndarray
    @classmethod
    def make(cls,k,d):
        k=np.array(k,dtype=float,copy=True);d=np.array(d,dtype=float,copy=True)
        k.flags.writeable=False;d.flags.writeable=False
        return cls(k,d)


class OpeningModel:
    def __init__(self,radial_pitch=.001,quadrature=4,penalty=None,Gc=None):
        self.parameters=json.loads((ROOT/'config/conditional_scenario.json').read_text())
        p=self.parameters
        if penalty is not None:p['penalty_Pa_per_m']=penalty
        if Gc is not None:p['Gc_J_m2']=Gc
        self.control=PeriodicSector(radial_pitch=radial_pitch)
        self.seam=SingleSeam(self.control,face_order=quadrature);b=self.seam.base
        self.base=b;self.X=self.seam.X;self.fixed_nodes=self.seam.fixed_nodes;self.free_nodes=self.seam.free_nodes
        self.free=(self.free_nodes[:,None]*3+np.arange(3)).ravel()
        self.contract=ConditionalFracture(p['peak_Pa'],p['Gc_J_m2'],p['penalty_Pa_per_m'],p['central_provenance'])
        self.bonds=self.seam.bind(self.contract)
        o=self.seam.geometry
        X=b.Xcell[o['pairs'][:,0]]
        self.interface_reference=np.einsum('qi,qia->qa',-o['J'][:,:18],X)
        self.starter=self.interface_reference[:,2]>b.config['outer_radius_m']-p['prepared_outer_slit_m']-1e-12
        maximum=np.where(self.starter,self.contract.delta_f,0.)
        self.initial_history=History.make(maximum,self.bonds._law(maximum))
        self.area=o['area'];self.initial_phi=float(self.area@self.phi(maximum))
        # Reassemble the unchanged per-cell bulk preconditioner on cut topology.
        K=self.control.base.cell_K;cell_dof=(self.seam.cell_map[...,None]*3+np.arange(3)).reshape(-1,54)
        rows=np.broadcast_to(cell_dof[:,:,None],K.shape).ravel();cols=np.broadcast_to(cell_dof[:,None,:],K.shape).ravel()
        full=coo_matrix((K.ravel(),(rows,cols)),shape=(len(b.X)*3,len(b.X)*3)).tocsc()
        self.bulk_K=(self.seam.T.T@full@self.seam.T).tocsc()
        pairs=o['pairs'];self.interface_nodes=np.concatenate((self.seam.cell_map[pairs[:,0]],self.seam.cell_map[pairs[:,1]]),axis=1)
        self.last_precondition_damage=None;self.factor=None

    def phi(self,k):
        return self.contract.Gc_J_m2*np.clip((k-self.contract.delta_0)/(self.contract.delta_f-self.contract.delta_0),0.,1.)
    def set_load(self,s):self.base.set_load(s)
    @property
    def load(self):return self.base.load
    def evaluate(self,u,accepted,work=False):
        if accepted.maximum.shape!=(self.bonds.N,) or not np.isfinite(accepted.maximum).all() or np.any(accepted.maximum<0):raise ValueError('Invalid immutable accepted history')
        if not np.allclose(accepted.damage,self.bonds._law(accepted.maximum),atol=2e-12,rtol=0):raise ValueError('Inconsistent immutable accepted history')
        Eb,gb,r=self.seam.bulk(u,work=work)
        _,_,opening=self.seam.interface(u,self.bonds,accepted.damage)
        k=np.maximum(accepted.maximum,opening);damage=self.bonds._law(k)
        Ei,gi,_=self.seam.interface(u,self.bonds,damage)
        additional=float(self.area@(self.phi(k)-self.phi(accepted.maximum)))
        total_new=float(self.area@self.phi(k))-self.initial_phi
        stored=Eb+Ei;g=gb+gi
        r.update(stored_energy_J=stored,bulk_energy_J=Eb,interface_energy_J=Ei,
                 additional_dissipation_J=additional,new_dissipation_from_initial_J=total_new,
                 incremental_potential_J=stored+additional,maximum_opening_m=float(opening.max()),
                 maximum_new_damage=float(np.max(damage[~self.starter])),
                 maximum_free_force_N=float(np.max(np.abs(g[self.free_nodes]))),prepared_slit_phi_J=self.initial_phi)
        r['fracture_dissipation_J']=total_new
        r.pop('damage',None)
        return stored+additional,g,History.make(k,damage),r,opening

    def precondition(self,damage):
        if self.factor is not None and np.max(abs(damage-self.last_precondition_damage))<1e-4:return self.factor
        J=self.seam.geometry['J'];weight=self.area*self.contract.penalty_Pa_per_m*(1-damage)
        K=weight[:,None,None]*J[:,:,None]*J[:,None,:]
        nodes=self.interface_nodes;rows=np.broadcast_to(nodes[:,:,None],K.shape).ravel();cols=np.broadcast_to(nodes[:,None,:],K.shape).ravel()
        scalar=coo_matrix((K.ravel(),(rows,cols)),shape=(len(self.base.X),len(self.base.X))).tocsr()
        interface=self.seam.T.T@kron(scalar,eye(3),format='csr')@self.seam.T
        full=(self.bulk_K+interface).tocsc();self.precondition_matrix=full[self.free,:][:,self.free]
        self.factor=splu(self.precondition_matrix);self.last_precondition_damage=damage.copy()
        return self.factor

    def equilibrate(self,start,accepted,deadline=None):
        self.last_trial=None
        p=self.parameters;u=start.copy();merit,g,history,r,opening=self.evaluate(u,accepted)
        counters=dict(evaluations=1,factorizations=0,line_search_rejections=0);pairs=[]
        for iteration in range(p['maximum_equilibrium_iterations']):
            self.last_trial=dict(displacement=u.copy(),maximum=history.maximum.copy(),damage=history.damage.copy(),report=dict(r),opening=opening.copy())
            if deadline is not None and time.monotonic()>deadline:raise TimeoutError('Station solve reached process budget; accepted state remains unchanged')
            oldfactor=self.factor;factor=self.precondition(history.damage)
            counters['factorizations']+=int(factor is not oldfactor)
            gf=g.ravel()[self.free];steepest=factor.solve(gf);metric=float(gf@steepest)
            relative=np.sqrt(max(metric,0.)/max(2*(r['stored_energy_J']+r['new_dissipation_from_initial_J']),1e-20))
            if relative<p['equilibrium_relative_energy_norm'] and r['maximum_free_force_N']<p['equilibrium_force_N']:
                merit,g,history,r,opening=self.evaluate(u,accepted,work=True)
                r.update(iterations=iteration,relative_energy_residual=relative,**counters)
                return u,history,r,opening
            q=gf.copy();alphas=[]
            for s,y,rho in reversed(pairs):
                a=rho*float(s@q);alphas.append(a);q-=a*y
            z=factor.solve(q)
            for (s,y,rho),a in zip(pairs,reversed(alphas)):
                z+=s*(a-rho*float(y@z))
            slope=float(gf@z)
            if slope<=0 or not np.isfinite(slope):pairs=[];z=steepest;slope=metric
            direction=np.zeros(u.size);direction[self.free]=-z;direction=direction.reshape(u.shape)
            step=1.;last_error=None
            while step>=1e-8:
                try:
                    trial=u+step*direction
                    value,gg,hh,rr,oo=self.evaluate(trial,accepted);counters['evaluations']+=1
                    if value<=merit-1e-4*step*slope:
                        s=step*direction.ravel()[self.free];y=gg.ravel()[self.free]-gf;sy=float(s@y)
                        if sy>1e-10*np.linalg.norm(s)*np.linalg.norm(y):
                            pairs.append((s.copy(),y.copy(),1/sy));pairs=pairs[-p['lbfgs_memory']:]
                        u=trial;merit=value;g=gg;history=hh;r=rr;opening=oo;break
                except ValueError as exc:last_error=str(exc)
                counters['line_search_rejections']+=1;step*=.5
            else:raise RuntimeError(f'Incremental-potential line search failed; residual={relative}, force_N={r["maximum_free_force_N"]}, domain={last_error}')
        raise RuntimeError(f'Equilibrium iteration limit; relative residual={relative}, force_N={r["maximum_free_force_N"]}, new damage={r["maximum_new_damage"]}')
