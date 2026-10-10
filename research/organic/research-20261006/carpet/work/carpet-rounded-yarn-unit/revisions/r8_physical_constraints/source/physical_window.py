"""Actual swept contacts, physical containment and unchanged movement disks."""
from pathlib import Path
import sys,math
import numpy as np
from scipy.sparse import coo_matrix
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parent/'r6_c2_bridge/source'))
from continue_bridge_loop import Window as ApproximateWindow,Field,UNIT,witness_search

class PhysicalSection:
    def __init__(self,reference,movement,radius):
        self.reference=reference;self.movement=movement;self.cap=1-radius;self.n=len(reference);self.counts=(self.n,self.n)
    def residual_jacobian(self,flat,jacobian=True):
        x=flat.reshape(self.n,2);norm=np.linalg.norm(x,axis=1);delta=x-self.reference;move=np.linalg.norm(delta,axis=1);raw=np.r_[norm-self.cap,move-self.movement];value=np.maximum(raw,0)
        if not jacobian:return value
        ids=np.arange(self.n);row=[];col=[];entry=[]
        for offset,gradient in ((0,x/np.maximum(norm[:,None],1e-30)),(self.n,delta/np.maximum(move[:,None],1e-30))):
            for axis in range(2):row.append(offset+ids);col.append(2*ids+axis);entry.append(gradient[:,axis]*(raw[offset:offset+self.n]>0))
        matrix=coo_matrix((np.concatenate(entry),(np.concatenate(row),np.concatenate(col))),shape=(2*self.n,2*self.n)).tocsr();matrix.eliminate_zeros();return value,matrix
    def assess(self,flat):
        x=flat.reshape(self.n,2);return {'maximum_center_radius_normalized':float(np.linalg.norm(x,axis=1).max()),'physical_center_radius_limit_normalized':self.cap,'maximum_movement_normalized':float(np.linalg.norm(x-self.reference,axis=1).max()),'movement_limit_normalized':self.movement}

class PhysicalWindow(ApproximateWindow):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.metric_diagnostics=self.parts
        self.parts=[PhysicalSection(self.refs[k],self.move,self.field.r/self.field.R) for k in range(2)]

def span_envelope_upper(field,s0,s1,p0,p1):
    # This minimal reproduction lies wholly inside the new quintic bridge.
    lo,hi,coefficients=next(part for part in field.spline.parts if part[0]<=s0 and s1<=part[1])
    ta=(s0-lo)/(hi-lo);h=(s1-s0)/(hi-lo);power=np.zeros_like(coefficients)
    for k in range(6):
        for j in range(k,6):power[k]+=coefficients[j]*math.comb(j,k)*ta**(j-k)*h**k
    difference=np.broadcast_to(-power[None],(len(p0),6,3)).copy();difference[:,0]+=p0;difference[:,1]+=p1-p0
    bezier=np.zeros_like(difference)
    for i in range(6):
        for k in range(i+1):bezier[:,i]+=difference[:,k]*math.comb(i,k)/math.comb(5,k)
    upper=np.linalg.norm(bezier,axis=2).max(1)+field.r+1e-12
    return {'maximum_continuous_surface_envelope_upper_m':float(upper.max()),'nominal_bundle_radius_m':field.R,'physical_containment_pass':bool(upper.max()<=field.R),'method':'Quintic Bernstein convex-hull bound for each straight body segment minus the exact bridge at the same parameter, plus physical fibre radius and1pm numerical allowance.'}
