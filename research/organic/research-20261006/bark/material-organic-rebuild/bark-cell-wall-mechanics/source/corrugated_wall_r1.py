"""Planar extensible discrete elastic wall strip, in SI units.

A resolved stress-free corrugation is shortened or lengthened between end plates.
Axial and angle-bending energies share this same material reference. End tangent
energies are the half-cell quadrature boundary terms for clamped wall normals.
This is a bounded mechanism test, not a calibrated three-dimensional cork law.
"""
from dataclasses import dataclass
import numpy as np
from scipy.linalg import eigh,solve

@dataclass(frozen=True)
class Parameters:
    height_m:float=43e-6
    wall_thickness_m:float=1e-6
    strip_width_m:float=1e-6
    corrugation_amplitude_m:float=2.8e-6
    corrugation_cycles:int=3
    wall_modulus_Pa:float=1e9
    wall_density_kg_m3:float=1150.

def wrap(a):return np.arctan2(np.sin(a),np.cos(a))

class CorrugatedWall:
    def __init__(self,segments=48,parameters=Parameters()):
        if segments<8 or parameters.height_m<=0 or parameters.wall_thickness_m<=0:raise ValueError('Invalid wall geometry')
        self.p=parameters;self.n=segments;self.L=parameters.height_m
        z=np.linspace(0,1,segments+1)
        self.X=np.c_[parameters.corrugation_amplitude_m/self.L*np.sin(2*np.pi*parameters.corrugation_cycles*z),z]
        d=np.diff(self.X,axis=0);self.l0=np.linalg.norm(d,axis=1);self.theta0=np.arctan2(d[:,1],d[:,0]);self.turn0=wrap(np.diff(self.theta0))
        self.A=parameters.wall_thickness_m*parameters.strip_width_m;self.I=parameters.strip_width_m*parameters.wall_thickness_m**3/12
        self.EA=parameters.wall_modulus_Pa*self.A;self.EI=parameters.wall_modulus_Pa*self.I;self.unit=self.EI/self.L
        self.axial=self.EA*self.L**2/self.EI/self.l0
        self.bending=1/((self.l0[:-1]+self.l0[1:])/2)
        self.end=2/np.array([self.l0[0],self.l0[-1]])
        self.volume_m3=self.A*self.L*self.l0.sum();self.mass_kg=self.volume_m3*parameters.wall_density_kg_m3
        # Incidence maps node variations to segment-vector variations.
        self.D=np.zeros((2*segments,2*(segments+1)))
        for i in range(segments):self.D[2*i:2*i+2,2*i:2*i+2]=-np.eye(2);self.D[2*i:2*i+2,2*i+2:2*i+4]=np.eye(2)
        self.C=np.zeros((segments,segments))
        for i,k in enumerate(self.bending):self.C[i:i+2,i:i+2]+=k*np.array([[1.,-1.],[-1.,1.]])
        self.C[0,0]+=self.end[0];self.C[-1,-1]+=self.end[1]

    def evaluate(self,x,rotation_rad=0.,hessian=False):
        x=np.asarray(x,dtype=float).reshape(self.n+1,2);d=np.diff(x,axis=0);length=np.linalg.norm(d,axis=1)
        if np.any(length<1e-7):raise ValueError('Degenerate wall segment')
        t=d/length[:,None];theta=np.arctan2(d[:,1],d[:,0]);turn=wrap(np.diff(theta)-self.turn0)
        ends=wrap(theta[[0,-1]]-self.theta0[[0,-1]]-rotation_rad)
        stretch=length-self.l0
        energy=.5*np.sum(self.axial*stretch**2)+.5*np.sum(self.bending*turn**2)+.5*np.sum(self.end*ends**2)
        q=np.zeros(self.n);q[:-1]-=self.bending*turn;q[1:]+=self.bending*turn;q[0]+=self.end[0]*ends[0];q[-1]+=self.end[1]*ends[1]
        jt=np.c_[-d[:,1],d[:,0]]/length[:,None]**2
        gd=self.axial[:,None]*stretch[:,None]*t+q[:,None]*jt
        g=(self.D.T@gd.ravel()).reshape(x.shape)
        parts=dict(axial_J=float(.5*np.sum(self.axial*stretch**2)*self.unit),bending_J=float((.5*np.sum(self.bending*turn**2)+.5*np.sum(self.end*ends**2))*self.unit),max_axial_strain=float(abs(stretch/self.l0).max()),max_bending_outer_fibre_increment=float(abs(turn/(self.L*(self.l0[:-1]+self.l0[1:])/2)*self.p.wall_thickness_m/2).max()),reference_volume_m3=self.volume_m3)
        if not hessian:return float(energy),g,parts
        Hd=np.zeros((2*self.n,2*self.n));Jt=np.zeros((self.n,2*self.n))
        for i in range(self.n):
            a,b=d[i];l=length[i]
            Htheta=np.array([[2*a*b,b*b-a*a],[b*b-a*a,-2*a*b]])/l**4
            Ha=self.axial[i]*((1-self.l0[i]/l)*np.eye(2)+self.l0[i]/l*np.outer(t[i],t[i]))
            Hd[2*i:2*i+2,2*i:2*i+2]=Ha+q[i]*Htheta;Jt[i,2*i:2*i+2]=jt[i]
        H=self.D.T@(Hd+Jt.T@self.C@Jt)@self.D
        return float(energy),g,H,parts

    def solve(self,engineering_strain,initial=None,maxiter=100):
        x=self.X.copy() if initial is None else np.asarray(initial).copy();old=x[-1,1]-x[0,1];target=1+engineering_strain
        x[:,1]=x[0,1]+(x[:,1]-x[0,1])*(target/old);x[0]=self.X[0];x[-1]=self.X[-1]+[0,engineering_strain]
        free=np.arange(2,2*self.n);trace=[];passed=False
        for it in range(maxiter):
            e,g,H,parts=self.evaluate(x,hessian=True);gf=g.ravel()[free];Hf=H[np.ix_(free,free)];scale=max(float(abs(np.diag(Hf)).max()),1.)
            eig=eigh(Hf,eigvals_only=True,subset_by_index=[0,0])[0];res=float(abs(gf).max());trace.append(dict(iteration=it,dimensionless_force_residual=res,minimum_Hessian_eigenvalue=float(eig)))
            if res<2e-8 and eig>0:passed=True;break
            # Regularize only the numerical search metric; physical law unchanged.
            shift=max(0.,-eig+1e-8*scale);dx=solve(Hf+shift*np.eye(len(free)),-gf,assume_a='pos');descent=float(gf@dx);alpha=min(1.,.02/max(float(abs(dx).max()),1e-30))
            accepted=False
            for ls in range(35):
                trial=x.copy();trial.ravel()[free]+=alpha*dx
                try:et,gt,_=self.evaluate(trial)
                except ValueError:et=np.inf
                if et<=e+1e-4*alpha*descent+1e-15:accepted=True;break
                alpha*=.5
            trace[-1].update(metric_shift=float(shift),line_search_steps=ls+1,step_fraction=float(alpha))
            if not accepted:break
            x=trial
        e,g,H,parts=self.evaluate(x,hessian=True);freeH=H[np.ix_(free,free)];mineig=eigh(freeH,eigvals_only=True,subset_by_index=[0,0])[0]
        reaction=float(g[-1,1]*self.unit/self.L)
        receipt=dict(engineering_strain=float(engineering_strain),passed=bool(passed),iterations=len(trace),energy_J=e*self.unit,end_reaction_N=reaction,wall_nominal_stress_Pa=reaction/self.A,force_residual_N=float(abs(g.ravel()[free]).max()*self.unit/self.L),dimensionless_force_residual=float(abs(g.ravel()[free]).max()),minimum_Hessian_eigenvalue=float(mineig),mass_kg=self.mass_kg,**parts,trace=trace)
        return x,receipt
