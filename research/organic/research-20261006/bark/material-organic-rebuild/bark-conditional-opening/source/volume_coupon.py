"""Conforming, fixed-reference retained-mismatch surrogate; no hidden damage."""
from pathlib import Path
import json, sys, time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[1]
ORGANIC = ROOT.parent
sys.path.insert(0, str(ORGANIC / 'causal-bark/volume-core/source'))
from quadratic_prism import QuadraticPrism, shapes
from birth_cell import cylindrical_nodes
sys.path.insert(0, str(ORGANIC / 'bark-directional-native-port/source'))
from native_directional import NativePrisms
from directional_candidate import Branch, DirectionalEnvelope


def configuration():
    return json.loads((ROOT/'config/experiment.json').read_text())


class DomainError(ValueError):
    pass


class Coupon:
    def __init__(self, pitch=None, precondition=True, config_override=None, factorize=True):
        started = time.monotonic()
        c = self.config = configuration() if config_override is None else config_override
        pitch = c['baseline_pitch_m'] if pitch is None else pitch
        self.pitch = pitch
        nx = round(c['width_m']/pitch); ny = round(c['height_m']/pitch)
        self.x = np.linspace(-c['width_m']/2, c['width_m']/2, nx+1)
        self.y = np.linspace(-c['height_m']/2, c['height_m']/2, ny+1)
        xx, yy = np.meshgrid(self.x, self.y)
        self.points = np.column_stack((xx.ravel(), yy.ravel()))
        tri = []
        for j in range(ny):
            for i in range(nx):
                a=j*(nx+1)+i; b=a+1; d=a+nx+1; e=d+1
                tri.extend(([a,b,e],[a,e,d]) if (i+j)%2==0 else ([a,b,d],[b,e,d]))
        self.tri = np.array(tri); self.T = len(tri); self.xy = self.points[self.tri]
        stocks = np.array(c['cohort_thickness_m'])
        subdivisions=np.rint(stocks/c.get('radial_pitch_m',min(stocks))).astype(int) if 'radial_pitch_m' in c else np.ones(3,int)
        thickness = np.repeat(stocks/subdivisions,subdivisions)
        self.subdivisions=subdivisions;self.numeric_layers=len(thickness)
        offsets=np.r_[0,np.cumsum(subdivisions)]*self.T
        self.material_slices=[slice(int(offsets[i]),int(offsets[i+1])) for i in range(3)]
        inner = c['outer_radius_m']-thickness.sum()
        mid = inner + thickness.cumsum()-thickness/2
        self.cells = []
        for layer in range(self.numeric_layers):
            for a in self.xy:
                material = a.copy(); material[:,0] /= c['outer_radius_m']
                X = cylindrical_nodes(material, mid[layer], thickness[layer])
                self.cells.append(QuadraticPrism(X, 1., 0., c['reference_density_kg_m3'], c['quadrature_order']))
        self.Xcell = np.array([cell.X for cell in self.cells])
        # Exact topological identification before forces and mass are assembled.
        _, first, inverse = np.unique(np.round(self.Xcell.reshape(-1,3),12), axis=0, return_index=True, return_inverse=True)
        self.X = self.Xcell.reshape(-1,3)[first].copy()
        self.map = inverse.reshape(-1,18)
        self.mapping_error = float(np.abs(self.X[self.map]-self.Xcell).max())
        assert self.mapping_error < 1e-14
        self.grad = np.array([cell.grad for cell in self.cells])
        self.weights = np.array([cell.weights for cell in self.cells])
        qX = np.einsum('cqi,cia->cqa', np.array([cell.N for cell in self.cells]), self.Xcell)
        theta = np.arctan2(qX[...,0], qX[...,2])
        radial=np.stack((np.sin(theta), np.zeros_like(theta), np.cos(theta)),axis=-1)
        tangent=np.stack((np.cos(theta), np.zeros_like(theta), -np.sin(theta)),axis=-1)
        axial=np.broadcast_to(np.array([0.,1.,0.]),radial.shape)
        self.axes=np.stack((radial,tangent,axial),axis=-1)
        assert np.max(np.abs(np.linalg.det(self.axes)-1)) < 1e-12
        self.Aprime=np.repeat(np.array(c['elastic_tangential_mismatch_at_unit_load']),subdivisions*self.T)[:,None,None,None]*tangent[..., :,None]*tangent[...,None,:]
        self.fixed_nodes=np.unique(self.map[:self.T,:6])
        self.free_nodes=np.setdiff1d(np.arange(len(self.X)),self.fixed_nodes)
        self.free=(self.free_nodes[:,None]*3+np.arange(3)).ravel()
        self.cell_dof=(self.map[...,None]*3+np.arange(3)).reshape(-1,54)
        # Condensed consistent mass M_global = S.T M_cell S, not averaging.
        M=np.array([cell.mass_scalar for cell in self.cells])
        self.mass_matrix=coo_matrix((M.ravel(),(np.broadcast_to(self.map[:,:,None],M.shape).ravel(), np.broadcast_to(self.map[:,None,:],M.shape).ravel())),shape=(len(self.X),len(self.X))).tocsr()
        self.mass=float(self.mass_matrix.sum())
        self.volume=float(self.weights.sum())
        self.laws=[]
        for layer in range(3):
            E=c['tangential_initial_modulus_Pa'][layer]; knee=c['tangential_knee_strain_assumed']
            H=(c['tangential_peak_nominal_Pa'][layer]-E*knee)/(c['tangential_measured_endpoint_strain'][layer]-knee)
            linear=lambda modulus: Branch(modulus,modulus,.2)
            tension=(linear(c['radial_tension_modulus_Pa_assumed']),Branch(E,H,knee),linear(c['axial_tension_modulus_Pa_assumed']))
            compression=tuple(Branch(v,c['compression_plateau_tangent_Pa_assumed'],c['compression_knee_strain_assumed'],c['compression_densification_modulus_Pa_assumed'],c['compression_densification_onset_assumed']) for v in c['compression_initial_modulus_Pa'])
            self.laws.append(DirectionalEnvelope(tension,compression,c['shear_modulus_Pa_assumed'],minimum_stretch=c['minimum_stretch_numerical'],maximum_stretch=c['maximum_stretch_numerical']))
        self.load=None
        self.set_load(0.)
        if precondition:
            normal_E=np.repeat(np.array([[c['radial_tension_modulus_Pa_assumed'],v,c['axial_tension_modulus_Pa_assumed']] for v in c['tangential_initial_modulus_Pa']]),subdivisions*self.T,axis=0)
            K=np.zeros((len(self.cells),54,54))
            local_grad=np.einsum('cqia,cqab->cqib', self.grad,self.axes)
            for i in range(3):
                B=(local_grad[...,i,None]*self.axes[...,None,:,i]).reshape(len(self.cells),-1,54)
                K+=np.einsum('cqi,cqj,cq,c->cij',B,B,self.weights,normal_E[:,i],optimize=True)
            for i,j in ((0,1),(0,2),(1,2)):
                B=(local_grad[...,j,None]*self.axes[...,None,:,i]+local_grad[...,i,None]*self.axes[...,None,:,j]).reshape(len(self.cells),-1,54)
                K+=c['shear_modulus_Pa_assumed']*np.einsum('cqi,cqj,cq->cij',B,B,self.weights,optimize=True)
            self.cell_K=K
            rows=np.broadcast_to(self.cell_dof[:,:,None],K.shape).ravel(); cols=np.broadcast_to(self.cell_dof[:,None,:],K.shape).ravel()
            self.K=coo_matrix((K.ravel(),(rows,cols)),shape=(3*len(self.X),3*len(self.X))).tocsc()
            if factorize:self.factor=splu(self.K[self.free,:][:,self.free])
        self.build_s=time.monotonic()-started

    def set_load(self,load):
        if self.load == load: return
        self.operators=[]
        for layer in range(3):
            sl=self.material_slices[layer]
            H0=load*self.Aprime[sl]
            A=np.eye(3)+H0
            transformed=np.einsum('cqia,cqab->cqib',self.grad[sl],A)
            self.operators.append(NativePrisms(transformed,self.weights[sl],H0,self.axes[sl],self.config['reference_density_kg_m3'],self.laws[layer]))
        self.load=load

    def evaluate(self,u,check_domain=True,work=False):
        local=u[self.map]
        energy=0.; grad=np.zeros_like(u); minimum=1.; maxima=[]; bulk_J=[]
        for layer,op in enumerate(self.operators):
            sl=self.material_slices[layer]
            E,g,d=op.evaluate(local[sl])
            stretches=d['stretches']
            cap=self.config['tangential_measured_endpoint_strain'][layer]
            tangential=float(np.max(stretches[...,1]-1))
            if check_domain and tangential > cap+1e-12:
                raise DomainError(f'cohort {layer}: tangential strain {tangential:.9g} exceeds reported endpoint {cap:.9g}')
            if check_domain and np.max(stretches[...,[0,2]]-1)>self.config['other_tension_domain_strain_assumed']:
                raise DomainError(f'cohort {layer}: unmeasured radial/axial tension {float(np.max(stretches[...,[0,2]]-1)):.9g} exceeds declared assumed scope {self.config["other_tension_domain_strain_assumed"]}')
            energy+=E; np.add.at(grad,self.map[sl],g)
            minima=float(d['minimum_det_F']); minimum=min(minimum,minima); maxima.append(tangential)
            bulk_J.append(d['J'])
        physical_H=np.einsum('cia,cqib->cqab',local-local.mean(axis=1,keepdims=True),self.grad)
        physical_F=physical_H+np.eye(3)
        physical_J=np.linalg.det(physical_F)
        if physical_J.min()<=0: raise DomainError('inverted physical deformation')
        report=dict(load=self.load,stored_energy_J=energy,maximum_tangential_strain_by_cohort=maxima,
                    minimum_elastic_J=minimum,minimum_physical_J=float(physical_J.min()),
                    physical_volume_m3=float(np.sum(self.weights*physical_J)),reference_mass_kg=self.mass,
                    maximum_free_force_N=float(np.max(np.abs(grad[self.free_nodes]))),damage=0.,fracture_dissipation_J=0.,damping_J=0.,kinetic_J=0.)
        if work:
            conjugate=0.
            for layer,op in enumerate(self.operators):
                sl=self.material_slices[layer]
                A=np.eye(3)+self.load*self.Aprime[sl]
                elastic_H=np.einsum('cqab,cqbd->cqad',physical_H[sl],A)+self.load*self.Aprime[sl]
                _,P,_=op.material.evaluate_incremental(elastic_H,self.axes[sl])
                dF=np.einsum('cqab,cqbd->cqad',physical_F[sl],self.Aprime[sl])
                conjugate+=float(np.sum(self.weights[sl]*np.einsum('cqab,cqab->cq',P,dF)))
            report['eigenstrain_conjugate_J']=conjugate
            report['mean_current_density_kg_m3']=self.mass/report['physical_volume_m3']
        return energy,grad,report

    def equilibrate(self,u):
        c=self.config; u=u.copy(); trace=[]
        E,g,r=self.evaluate(u)
        for iteration in range(c['maximum_equilibrium_iterations']):
            gf=g.ravel()[self.free]; z=self.factor.solve(gf)
            metric=float(gf@z)
            relative=np.sqrt(max(metric,0.)/max(2*E,1e-20))
            if relative<c['equilibrium_relative_energy_norm_tolerance'] and r['maximum_free_force_N']<c['equilibrium_maximum_force_N']:
                return u,dict(iterations=iteration,relative_energy_residual=relative,line_search_evaluations=sum(trace),**self.evaluate(u,work=True)[2])
            direction=np.zeros(u.size); direction[self.free]=-z; direction=direction.reshape(u.shape)
            step=1.; accepted=False
            ntrial=0; last_invalid=None
            while step>=c['line_search_minimum_step']:
                ntrial+=1
                try:
                    trial=u+step*direction
                    e,gg,rr=self.evaluate(trial)
                    if e<=E-c['line_search_armijo']*step*metric:
                        u=trial; E=e; g=gg; r=rr; accepted=True; break
                except ValueError as error:
                    last_invalid=str(error)
                step*=.5
            trace.append(ntrial)
            if not accepted:
                raise RuntimeError(f'line search failed at {iteration}, residual {relative}, force {r["maximum_free_force_N"]}; last invalid trial: {last_invalid}')
        raise RuntimeError(f'equilibrium iteration limit, residual {relative}, force {r["maximum_free_force_N"]}')

    def kinetic(self,v):
        return .5*float(np.sum(v*(self.mass_matrix@v)))


if __name__=='__main__':
    m=Coupon(); print(json.dumps(dict(nodes=len(m.X),cells=len(m.cells),mass=m.mass,build_s=m.build_s)))
