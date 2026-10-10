"""Complete swept endpoint constraints in dimensionless geometry."""
from pathlib import Path
import sys,time
import numpy as np
from scipy.sparse import coo_matrix,vstack
from scipy.optimize import least_squares
ROOT=Path(__file__).resolve().parents[1];R2=ROOT.parent/'r2_section_feasibility_frozen'
sys.path.insert(0,str(R2/'source'))
from section_feasibility import SectionConstraints

def segment_closest(a,b,c,d):
    u=b-a;v=d-c;w=a-c;aa=np.sum(u*u,axis=-1);bb=np.sum(u*v,axis=-1);cc=np.sum(v*v,axis=-1);dd=np.sum(u*w,axis=-1);ee=np.sum(v*w,axis=-1)
    det=aa*cc-bb*bb;product=aa*cc;relative=np.divide(det,product,out=np.zeros_like(det),where=product>0);epsilon=64*np.finfo(float).eps
    ss=[];tt=[]
    for s in (0.,1.):ss.append(np.full_like(aa,s));tt.append(np.clip((ee+s*bb)/np.maximum(cc,np.finfo(float).tiny),0,1))
    for t in (0.,1.):ss.append(np.clip((t*bb-dd)/np.maximum(aa,np.finfo(float).tiny),0,1));tt.append(np.full_like(cc,t))
    regular=relative>epsilon
    si=np.divide(bb*ee-cc*dd,det,out=np.zeros_like(det),where=regular);ti=np.divide(aa*ee-bb*dd,det,out=np.zeros_like(det),where=regular);ss.append(si);tt.append(ti)
    s=np.stack(ss);t=np.stack(tt);difference=w[None]+s[...,None]*u[None]-t[...,None]*v[None];square=np.sum(difference*difference,axis=-1)
    square[4]=np.where(regular&(si>=0)&(si<=1)&(ti>=0)&(ti<=1),square[4],np.inf)
    which=np.argmin(square,axis=0);flat=np.arange(aa.size);shape=aa.shape
    selected_s=s.reshape(5,-1)[which.ravel(),flat].reshape(shape);selected_t=t.reshape(5,-1)[which.ravel(),flat].reshape(shape);delta=difference.reshape(5,-1,3)[which.ravel(),flat].reshape(shape+(3,));distance=np.linalg.norm(delta,axis=-1)
    # For an unresolved near-parallel interior, subtract a conservative angular
    # distance bound rather than declaring an endpoint estimate exact.
    angular=np.where(regular,0.,np.minimum(np.sqrt(aa),np.sqrt(cc))*np.sqrt(np.maximum(relative,0)+epsilon))
    roundoff=epsilon*(np.sqrt(aa)+np.sqrt(cc)+np.linalg.norm(w,axis=-1));error=angular+roundoff
    normal=delta/np.maximum(distance[...,None],1e-30)
    return distance,selected_s,selected_t,normal,error,which

def gradient_checks():
    tests=[]
    for x in (-.1,0.,.5,1.,1.1):
        a=np.array([0.,0,0]);b=np.array([1.,0,0]);c=np.array([x,-1,.2]);d=np.array([x,1,.2]);z=np.r_[a,b,c,d];dist,s,t,n,err,regime=segment_closest(a,b,c,d)
        analytic=np.r_[(1-s)*n,s*n,-(1-t)*n,-t*n];numeric=[];h=1e-6
        for i in range(12):
            plus=z.copy();minus=z.copy();plus[i]+=h;minus[i]-=h
            numeric.append((segment_closest(*plus.reshape(4,3))[0]-segment_closest(*minus.reshape(4,3))[0])/(2*h))
        error=float(np.max(abs(analytic-numeric)));assert error<5e-6
        tests.append({'crossing_parameter':x,'closest_s':float(s),'closest_t':float(t),'gradient_error':error,'selected_feature':int(regime)})
    a=np.array([0.,0,0]);b=np.array([1.,0,0]);c=np.array([.2,.3,1e-9]);d=np.array([.9,.3,2e-9]);scales=[]
    for scale in (1e-6,1.,1e6):
        dist,s,t,n,err,regime=segment_closest(a*scale,b*scale,c*scale,d*scale);scales.append({'scale':scale,'distance_rescaled':float(dist/scale),'error_bound_rescaled':float(err/scale),'feature':int(regime)})
    assert np.ptp([x['distance_rescaled'] for x in scales])<1e-12
    return {'endpoint_interior_switch_tests':tests,'scale_invariance':scales,'relative_determinant_threshold':64*np.finfo(float).eps,'passed':True,'limits':'At coincident or nonunique parallel minima the distance need not have a unique ordinary gradient; the selected valid branch is used, and the final conservative distance bound remains authoritative.'}

class SweptConstraints(SectionConstraints):
    def __init__(self,reference,k,w,r,wrap_r,margin,movement,previous,center,basis,obstacle_a,obstacle_b,obstacle_radius):
        super().__init__(reference,k,w,r,wrap_r,margin,movement=movement)
        self.previous=previous;self.center=center;self.basis=basis
        self.oa=obstacle_a;self.ob=obstacle_b;self.oradius=obstacle_radius
        self.fi=np.repeat(np.arange(self.n),len(obstacle_a));self.oi=np.tile(np.arange(len(obstacle_a)),self.n)
    def residual_jacobian(self,flat,jacobian=True):
        base=super().residual_jacobian(flat,jacobian);value,matrix=base if jacobian else (base,None)
        x=flat.reshape(self.n,2);end=self.center+x@self.basis.T
        dist,s,t,normal,error,_=segment_closest(self.previous[self.i],end[self.i],self.previous[self.j],end[self.j]);raw=2*self.r+self.margin-(dist-error);res=np.maximum(raw,0)
        distance,parameter,_,direction,bound,_=segment_closest(self.previous[self.fi],end[self.fi],self.oa[self.oi],self.ob[self.oi]);raw_fixed=self.r+self.oradius[self.oi]+self.margin-(distance-bound);fixed=np.maximum(raw_fixed,0)
        if not jacobian:return np.r_[value,res,fixed]
        rows=[];columns=[];entries=[]
        for ids,g,offset,active in ((self.i,-s[:,None]*(normal@self.basis),0,raw>0),(self.j,t[:,None]*(normal@self.basis),0,raw>0),(self.fi,-parameter[:,None]*(direction@self.basis),len(res),raw_fixed>0)):
            for dim in range(2):rows.append(np.arange(len(ids))+offset);columns.append(ids*2+dim);entries.append(g[:,dim]*active)
        extra=coo_matrix((np.concatenate(entries),(np.concatenate(rows),np.concatenate(columns))),shape=(len(res)+len(fixed),2*self.n)).tocsr();extra.eliminate_zeros()
        return np.r_[value,res,fixed],vstack((matrix,extra),format='csr')
    def assess_sweeps(self,flat):
        x=flat.reshape(self.n,2);end=self.center+x@self.basis.T
        dist,_,_,_,error,_=segment_closest(self.previous[self.i],end[self.i],self.previous[self.j],end[self.j]);distance,_,_,_,bound,_=segment_closest(self.previous[self.fi],end[self.fi],self.oa[self.oi],self.ob[self.oi])
        return {'body_swept_gap_lower_bound_normalized':float((dist-error-2*self.r).min()),'wrap_swept_gap_lower_bound_normalized':float((distance-bound-self.r-self.oradius[self.oi]).min()),'maximum_conditioning_error_normalized':float(max(error.max(),bound.max())),'body_pairs':len(self.i),'fixed_obstacle_pairs':len(self.fi),'all_pairs_complete':True}

class Stop(Exception):pass
def solve(system,reference,tolerance,seconds=30,max_evaluations=160):
    start=time.monotonic();best=reference.ravel().copy();score=np.inf;count=0;cached=None;cached_pair=None;reason='unstarted'
    def evaluate(x):
        nonlocal best,score,count,cached,cached_pair,reason
        if cached is not None and np.array_equal(x,cached):return cached_pair
        if time.monotonic()-start>seconds:reason='wall_limit';raise Stop
        values,jac=system.residual_jacobian(x);count+=1;maximum=float(values.max())
        if maximum<score:score=maximum;best=x.copy()
        cached=x.copy();cached_pair=(values,jac)
        if maximum<=tolerance:reason='feasible';raise Stop
        return values,jac
    try:
        answer=least_squares(lambda x:evaluate(x)[0],reference.ravel(),jac=lambda x:evaluate(x)[1],bounds=((reference-system.movement).ravel(),(reference+system.movement).ravel()),method='trf',tr_solver='lsmr',tr_options={'atol':1e-10,'btol':1e-10,'maxiter':300},x_scale='jac',max_nfev=max_evaluations,ftol=1e-12,xtol=1e-12,gtol=1e-12)
        reason='optimizer_stopped_without_feasibility'
    except Stop:pass
    final=system.residual_jacobian(best,False);return best.reshape(reference.shape),{'status':reason,'accepted':bool(final.max()<=tolerance),'maximum_distance_residual_normalized':float(final.max()),'section':system.assess(best),'sweeps':system.assess_sweeps(best),'seconds':time.monotonic()-start,'evaluations':count,'limits':{'wall_seconds':seconds,'evaluations':max_evaluations,'distance_residual_tolerance_normalized':tolerance}}
