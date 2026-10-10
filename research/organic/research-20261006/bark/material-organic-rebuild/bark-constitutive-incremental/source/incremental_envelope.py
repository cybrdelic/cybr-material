"""Cancellation-resistant evaluation of the same synthetic directional energy.

H=F-I is the explicit input so tiny perturbations are not lost by reconstructing
world positions. The normalized material Gram determinant supplies angular
energy without subtracting two almost-equal logarithms. No material law changes.
"""
from pathlib import Path
import sys,numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'bark-constitutive-upgrade/source'))
from directional_candidate import DirectionalEnvelope

class IncrementalEnvelope(DirectionalEnvelope):
    def evaluate_incremental(self,H):
        H=np.asarray(H,dtype=float)
        if H.shape[-2:]!=(3,3) or not np.isfinite(H).all():raise ValueError('Expected finite3x3 incremental gradient')
        F=np.eye(3)+H;J=np.linalg.det(F)
        if np.any(J<=0):raise ValueError('Inverted or singular deformation')
        A=self.axes;HA=H@A;B=A+HA
        delta_sq=2*np.sum(A*HA,axis=-2)+np.sum(HA*HA,axis=-2)
        if np.any(delta_sq<=-1):raise ValueError('Degenerate material stretch')
        stretch=np.sqrt(1+delta_sq);e=delta_sq/(stretch+1)
        if np.any(stretch<self.minimum_stretch) or np.any(stretch>self.maximum_stretch):raise ValueError('Outside declared experimental stretch domain')
        axial=np.zeros_like(J);stress=np.zeros_like(e)
        for i in range(3):
            ut,pt,_=self.tension[i].evaluate(np.maximum(e[...,i],0));uc,pc,_=self.compression[i].evaluate(np.maximum(-e[...,i],0));axial+=ut+uc;stress[...,i]=pt-pc
        # The reference axes are orthonormal by contract; form only their
        # incremental metric change, avoiding spurious reference cross terms.
        deltaC=A.T@HA+HA.swapaxes(-2,-1)@A+HA.swapaxes(-2,-1)@HA
        off=deltaC/(stretch[..., :,None]*stretch[...,None,:])
        off=off.copy();off[...,np.arange(3),np.arange(3)]=0.
        c01,c02,c12=off[...,0,1],off[...,0,2],off[...,1,2]
        det_delta=2*c01*c02*c12-c01*c01-c02*c02-c12*c12
        if np.any(det_delta<=-1):raise ValueError('Degenerate normalized material frame')
        angular=-.5*self.G*np.log1p(det_delta)
        normalized=B/stretch[...,None,:];inverse=np.linalg.inv(np.eye(3)+off)
        angular_P=self.G*(normalized@(off@inverse)/stretch[...,None,:])@A.T
        P=(B*(stress/stretch)[...,None,:])@A.T+angular_P
        return axial+angular,P,dict(J=J,stretches=stretch,axial_energy_J_m3=axial,angular_energy_J_m3=angular)

    def evaluate(self,F):return self.evaluate_incremental(np.asarray(F,dtype=float)-np.eye(3))
