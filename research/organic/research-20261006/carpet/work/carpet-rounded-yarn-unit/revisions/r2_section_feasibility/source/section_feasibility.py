"""Residual-driven dimensionless geometric feasibility; no mechanics claims."""
import time
import numpy as np
from scipy.optimize import least_squares
from scipy.sparse import coo_matrix

def metric_gradient(x,y,k,w):
    d=x-y;m=(x+y)*.5;v=w*np.stack((-m[...,1],m[...,0]),axis=-1)
    axial=1-np.sum(k*m,axis=-1);den=axial*axial+np.sum(v*v,axis=-1);u=np.sum(v*d,axis=-1)
    squared=np.sum(d*d,axis=-1)-u*u/den;distance=np.sqrt(np.maximum(squared,1e-28))
    gd=2*d-2*(u/den)[...,None]*v
    jt_d=np.stack((d[...,1],-d[...,0]),axis=-1);jt_v=np.stack((v[...,1],-v[...,0]),axis=-1)
    gm=-2*(u*w/den)[...,None]*jt_d+(u*u/(den*den))[...,None]*(-2*axial[...,None]*k+2*w*jt_v)
    return distance,(gd+.5*gm)/(2*distance[...,None]),(-gd+.5*gm)/(2*distance[...,None])

class SectionConstraints:
    def __init__(self,reference,k,w,r,wrap_r,margin,cap=.91,movement=.1):
        self.reference=np.asarray(reference,dtype='f8');self.k=np.asarray(k);self.w=w;self.r=r;self.wrap_r=wrap_r;self.margin=margin;self.cap=cap;self.movement=movement;self.n=len(reference)
        self.i,self.j=np.triu_indices(self.n,1);self.wi=np.repeat(np.arange(self.n),6);self.wj=np.tile(np.arange(6),self.n)
        angle=np.arange(6)*np.pi/3;self.wraps=.99*np.stack((np.cos(angle),np.sin(angle)),axis=-1)
        self.counts=(len(self.i),len(self.wi),self.n,self.n)

    def residual_jacobian(self,flat,jacobian=True):
        x=flat.reshape(self.n,2);body,gi,gj=metric_gradient(x[self.i],x[self.j],self.k,self.w)
        wrap,wg,_=metric_gradient(x[self.wi],self.wraps[self.wj],self.k,self.w)
        norm=np.linalg.norm(x,axis=1);move=x-self.reference;distance=np.linalg.norm(move,axis=1)
        values=[2*self.r+self.margin-body,self.r+self.wrap_r+self.margin-wrap,norm-self.cap,distance-self.movement]
        residual=np.maximum(np.concatenate(values),0)
        if not jacobian:return residual
        rows=[];columns=[];entries=[];offset=0
        def add(ids,gradient,active):
            for coordinate in range(2):
                rows.append(np.arange(len(ids))+offset);columns.append(2*ids+coordinate);entries.append(gradient[:,coordinate]*active)
        active=values[0]>0;add(self.i,-gi,active);add(self.j,-gj,active);offset+=len(self.i)
        add(self.wi,-wg,values[1]>0);offset+=len(self.wi)
        add(np.arange(self.n),x/np.maximum(norm[:,None],1e-30),values[2]>0);offset+=self.n
        add(np.arange(self.n),move/np.maximum(distance[:,None],1e-30),values[3]>0)
        matrix=coo_matrix((np.concatenate(entries),(np.concatenate(rows),np.concatenate(columns))),shape=(len(residual),2*self.n)).tocsr();matrix.eliminate_zeros()
        return residual,matrix

    def assess(self,flat):
        x=flat.reshape(self.n,2);body,_,_=metric_gradient(x[self.i],x[self.j],self.k,self.w);wrap,_,_=metric_gradient(x[self.wi],self.wraps[self.wj],self.k,self.w)
        euclidean=np.linalg.norm(x[self.i]-x[self.j],axis=1);norm=np.linalg.norm(x,axis=1);move=np.linalg.norm(x-self.reference,axis=1)
        violation={'body':max(float((2*self.r+self.margin-body).max()),0),'wrap':max(float((self.r+self.wrap_r+self.margin-wrap).max()),0),'envelope':max(float(norm.max()-self.cap),0),'movement':max(float(move.max()-self.movement),0)}
        return {'maximum_constraint_violation_normalized':max(violation.values()),'constraint_violations_normalized':violation,'minimum_body_metric_gap_normalized':float(body.min()-2*self.r),'minimum_body_euclidean_gap_normalized':float(euclidean.min()-2*self.r),'minimum_wrap_metric_gap_normalized':float(wrap.min()-self.r-self.wrap_r),'maximum_radius_normalized':float(norm.max()),'maximum_movement_normalized':float(move.max()),'minimum_guide_axial_tangent':float((1-x@self.k).min())}

def validate_gradients():
    k=np.array([.58,-.26]);w=1.49;x=np.array([.28,.14]);y=np.array([.205,.105]);value,gx,gy=metric_gradient(x,y,k,w);analytic=np.r_[gx,gy];z=np.r_[x,y];h=1e-6;numeric=[]
    for i in range(4):
        plus=z.copy();minus=z.copy();plus[i]+=h;minus[i]-=h
        numeric.append((metric_gradient(plus[:2],plus[2:],k,w)[0]-metric_gradient(minus[:2],minus[2:],k,w)[0])/(2*h))
    metric_error=float(np.max(abs(analytic-numeric)));assert metric_error<1e-7
    reference=np.array([[.7,.2],[.8,.2]]);system=SectionConstraints(reference,k,w,.045,.05,.006,cap=.91,movement=.09);z=np.array([[.90,.30],[.88,.28]]).ravel();_,jac=system.residual_jacobian(z);numerical=np.empty(jac.shape)
    for i in range(4):
        plus=z.copy();minus=z.copy();plus[i]+=h;minus[i]-=h
        numerical[:,i]=(system.residual_jacobian(plus,False)-system.residual_jacobian(minus,False))/(2*h)
    residual_error=float(np.max(abs(jac.toarray()-numerical)));assert residual_error<1e-7
    return {'metric_midpoint_curvature_twist_max_gradient_error':metric_error,'full_active_residual_max_gradient_error':residual_error,'finite_difference_step_normalized':h,'passed':True}

class BoundedStop(Exception):pass

def solve_section(reference,k,w,r,wrap_r,margin,movement,seconds=30,max_evaluations=160,tolerance=1e-6):
    system=SectionConstraints(reference,k,w,r,wrap_r,margin,movement=movement);start=time.monotonic();best=reference.ravel().copy();best_score=np.inf;evaluations=0;status='unstarted';cache_x=None;cache_value=None;cache_jac=None
    def evaluate(x):
        nonlocal best,best_score,evaluations,status,cache_x,cache_value,cache_jac
        if cache_x is not None and np.array_equal(x,cache_x):return cache_value,cache_jac
        if time.monotonic()-start>seconds:status='wall_limit';raise BoundedStop
        value,jac=system.residual_jacobian(x);evaluations+=1;score=float(value.max())
        if score<best_score:best_score=score;best=x.copy()
        cache_x=x.copy();cache_value=value;cache_jac=jac
        if score<=tolerance:status='feasible';raise BoundedStop
        return value,jac
    try:
        result=least_squares(lambda x:evaluate(x)[0],reference.ravel(),jac=lambda x:evaluate(x)[1],bounds=((reference-movement).ravel(),(reference+movement).ravel()),method='trf',tr_solver='lsmr',tr_options={'atol':1e-10,'btol':1e-10,'maxiter':300},max_nfev=max_evaluations,ftol=1e-12,xtol=1e-12,gtol=1e-12,x_scale='jac')
        status='optimizer_stopped_without_feasibility';message=result.message
    except BoundedStop:message=status
    assessment=system.assess(best);accepted=assessment['maximum_constraint_violation_normalized']<=tolerance
    report={'status':status,'accepted':accepted,'optimizer_message':message,'evaluations':evaluations,'wall_seconds':time.monotonic()-start,'limits':{'wall_seconds':seconds,'evaluations':max_evaluations,'movement_normalized':movement,'acceptance_tolerance_normalized':tolerance},'initial':system.assess(reference.ravel()),'final':assessment,'authoritative_acceptance':'Body, fixed wraps, circular envelope and movement all meet residual tolerance. No clamp follows the solve. Failure returns rejected state and must not advance or interpolate.'}
    return best.reshape(reference.shape),report
