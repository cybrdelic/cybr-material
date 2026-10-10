"""One original construction-span half-width C2 bridge at each cap join."""
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parent/'r4_bounded_loop/source'))
from advance_loop import Field as OriginalField,UNIT

class LocalC2Bridge:
    def __init__(self,original,joins,half_width):
        self.original=original;self.parts=[]
        for join in joins:
            lo=join-half_width;hi=join+half_width;h=hi-lo;p0=original(lo);p1=original(hi);v0=original(lo,1)*h;v1=original(hi,1)*h;a0=original(lo,2)*h*h;a1=original(hi,2)*h*h
            c0=p0;c1=v0;c2=.5*a0;D=p1-c0-c1-c2;V=v1-c1-2*c2;A=a1-2*c2
            coefficients=np.stack((c0,c1,c2,10*D-4*V+.5*A,-15*D+7*V-A,6*D-3*V+.5*A))
            self.parts.append((lo,hi,coefficients))
    def __call__(self,s,nu=0):
        scalar=np.ndim(s)==0;values=np.atleast_1d(s).astype('f8');out=np.asarray(self.original(values,nu)).copy()
        for lo,hi,coefficients in self.parts:
            use=(values>=lo)&(values<=hi)
            if not use.any():continue
            c=coefficients.copy()
            for _ in range(nu):c=c[1:]*np.arange(1,len(c))[:,None]
            parameter=(values[use]-lo)/(hi-lo);out[use]=(parameter[:,None]**np.arange(len(c))[None])@c/(hi-lo)**nu
        return out[0] if scalar else out

class BridgeField(OriginalField):
    def __init__(self,index,data):
        super().__init__(index,data);self.original_spline=self.spline;arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(self.guide,axis=0),axis=1))];self.joins=np.array([arc[512],arc[1536]]);self.half_width=self.length/96;self.spline=LocalC2Bridge(self.original_spline,self.joins,self.half_width)
    def physical_twist_rate(self,s):return self.rate/np.linalg.norm(self.spline(s,1),axis=-1)
    def component_arrays(self,s):
        s=np.asarray(s);c=self.spline(s);v=self.spline(s,1);t=v/np.linalg.norm(v,axis=-1)[:,None];inside=np.cross(self.normal,t);inside/=np.linalg.norm(inside,axis=-1)[:,None];a=inside*np.cos(self.alpha)+self.normal*np.sin(self.alpha);b=np.cross(t,a);phase=self.phase0+self.rate*s;e1=np.cos(phase)[:,None]*a+np.sin(phase)[:,None]*b;e2=-np.sin(phase)[:,None]*a+np.cos(phase)[:,None]*b
        return c,e1,e2
    def wrapper_array(self,s):
        c,e1,e2=self.component_arrays(s);angle=np.arange(6)*np.pi/3
        return c[None]+.99*self.R*(np.cos(angle)[:,None,None]*e1[None]+np.sin(angle)[:,None,None]*e2[None])
