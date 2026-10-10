"""Error-bounded representation of frozen uniform Catmull-Rom fibre curves.

No material construction or crown repair occurs here. Key counts adapt per fibre.
Continuous correspondence bounds come from restricted cubic Bezier differences.
A derivative Bernstein bound also limits the displacement of swept radius circles.
"""
import numpy as np
from functools import lru_cache

MENU=(24,32,40,48,56,64,80,96,112,128,144)
POSITION_LIMIT=0.75e-6
SURFACE_LIMIT=1.0e-6

def bezier(p):
    prev=np.concatenate((p[:,:1],p[:,:-1]),axis=1)
    nex=np.concatenate((p[:,1:],p[:,-1:]),axis=1)
    tangent=(nex-prev)*.5
    return np.stack((p[:,:-1],p[:,:-1]+tangent[:,:-1]/3,
                     p[:,1:]-tangent[:,1:]/3,p[:,1:]),axis=2)

def basis(t):
    t=np.asarray(t);u=1-t
    return np.stack((u**3,3*u*u*t,3*u*t*t,t**3),axis=-1)

def deriv_basis(t):
    t=np.asarray(t);u=1-t
    return np.stack((-3*u*u,3*u*u-6*u*t,6*u*t-3*t*t,3*t*t),axis=-1)

def second_basis(t):
    t=np.asarray(t)
    return np.stack((6-6*t,-12+18*t,6-18*t,6*t),axis=-1)

def restrict_matrix(a,b):
    aa=basis(a);bb=basis(b);d=(b-a)[...,None]/3
    return np.stack((aa,aa+d*deriv_basis(a),bb-d*deriv_basis(b),bb),axis=-2)

@lru_cache(None)
def geometric_layout(n,m,subdiv=4):
    knots=np.asarray(m,dtype='f8') if isinstance(m,tuple) else np.linspace(0,n-1,m)
    count=len(knots);steps=np.diff(knots)
    merged=np.sort(np.r_[np.arange(n,dtype='f8'),knots]);merged=merged[np.r_[True,np.diff(merged)>1e-9]]
    ts=(merged[:-1,None]+np.arange(subdiv)[None,:]/subdiv*np.diff(merged)[:,None]).ravel();ts=np.r_[ts,n-1]
    oi=np.minimum(np.floor(ts).astype('i4'),n-2);ci=np.clip(np.searchsorted(knots,ts,side='right')-1,0,count-2)
    mid=(ts[:-1]+ts[1:])*.5;oin=np.minimum(np.floor(mid).astype('i4'),n-2);cin=np.clip(np.searchsorted(knots,mid,side='right')-1,0,count-2)
    v=(ts-knots[ci])/steps[ci]
    return ts,oi,ci,oin,cin,restrict_matrix(ts[:-1]-oin,ts[1:]-oin),v

@lru_cache(None)
def layout(n,m):
    knots=np.linspace(0,n-1,m);step=(n-1)/(m-1)
    qi=np.minimum(np.floor(knots).astype('i4'),n-2);qt=knots-qi
    merged=np.sort(np.r_[np.arange(n,dtype='f8'),knots])
    merged=merged[np.r_[True,np.diff(merged)>1e-9]]
    a=merged[:-1];b=merged[1:];mid=(a+b)*.5
    oi=np.minimum(np.floor(mid).astype('i4'),n-2)
    ci=np.minimum(np.floor(mid/step).astype('i4'),m-2)
    return qi,basis(qt),oi,ci,restrict_matrix(a-oi,b-oi),restrict_matrix(a/step-ci,b/step-ci)

def candidate_and_bound(p,radius,m,details=False,q_override=None):
    count,n,_=p.shape
    if m==n or (isinstance(m,tuple) and len(m)==n):
        values=(p.astype('f4'),np.zeros(count),np.zeros(count),np.zeros(count),np.zeros(count,dtype=bool))
        return values+({'failed_candidate_spans':[]},) if details else values
    src=bezier(p.astype('f8'))
    if q_override is not None:q=q_override.astype('f4')
    elif isinstance(m,tuple):q=p[:,list(m)].astype('f4')
    else:
        qi,qbasis,oi,ci,to,tc=layout(n,m)
        q=np.einsum('ka,ckad->ckd',qbasis,src[:,qi],optimize=True).astype('f4')
    assert np.array_equal(q[:,0],p[:,0]) and np.array_equal(q[:,-1],p[:,-1])
    reduced=bezier(q.astype('f8'))
    # Match geometry rather than assigning equal parameter speed. Within each
    # candidate span, closest-parameter endpoints define a continuous monotone
    # piecewise-affine correspondence. Its cubic restrictions retain a full bound.
    ts,oi,ci,oin,cin,to,initial_v=geometric_layout(n,m)
    target=np.einsum('ka,ckad->ckd',basis(ts-oi),src[:,oi],optimize=True)
    selected=reduced[:,ci];v=np.broadcast_to(initial_v,(count,len(ts))).copy()
    for _ in range(8):
        point=np.einsum('cka,ckad->ckd',basis(v),selected,optimize=True)
        d1=np.einsum('cka,ckad->ckd',deriv_basis(v),selected,optimize=True)
        d2=np.einsum('cka,ckad->ckd',second_basis(v),selected,optimize=True)
        denom=np.sum(d1*d1+(point-target)*d2,axis=-1)
        change=np.sum((point-target)*d1,axis=-1)/np.maximum(denom,1e-24)
        change[np.linalg.norm(point-target,axis=-1)<1e-12]=0
        v=np.clip(v-np.clip(change,-.2,.2),0,1)
    boundary=(np.abs(initial_v)<1e-9)|(np.abs(initial_v-1)<1e-9)
    v[:,boundary]=initial_v[boundary]
    absolute=v+ci[None]
    va=absolute[:,:-1]-cin[None];vb=absolute[:,1:]-cin[None]
    monotone=(vb>=va-1e-10)&(va>=-1e-9)&(vb<=1+1e-9)
    a=np.einsum('kab,ckbd->ckad',to,src[:,oin],optimize=True)
    b=np.einsum('ckab,ckbd->ckad',restrict_matrix(va,vb),reduced[:,cin],optimize=True)
    position=np.linalg.norm(a-b,axis=-1).max(axis=-1)
    exact=np.all(src[:,oin]==reduced[:,cin],axis=(2,3))
    # Derivative vectors are quadratic Bernstein polynomials in the same
    # restricted interval coordinate. Their cross/dot products are quartics.
    da=3*np.diff(a,axis=2);db=3*np.diff(b,axis=2)
    amid=(da[:,:,0]+2*da[:,:,1]+da[:,:,2])*.25
    bmid=(db[:,:,0]+2*db[:,:,1]+db[:,:,2])*.25
    aref=amid/np.maximum(np.linalg.norm(amid,axis=-1)[:,:,None],1e-30)
    bref=bmid/np.maximum(np.linalg.norm(bmid,axis=-1)[:,:,None],1e-30)
    lower_a=np.sum(da*aref[:,:,None],axis=-1).min(axis=-1)
    lower_b=np.sum(db*bref[:,:,None],axis=-1).min(axis=-1)
    max_cross=np.zeros_like(lower_a);min_dot=np.full_like(lower_a,np.inf)
    comb2=(1,2,1);comb4=(1,4,6,4,1)
    for k in range(5):
        cross=np.zeros_like(amid);dot=np.zeros_like(lower_a)
        for i in range(3):
            j=k-i
            if 0<=j<3:
                weight=comb2[i]*comb2[j]/comb4[k]
                cross+=weight*np.cross(da[:,:,i],db[:,:,j])
                dot+=weight*np.sum(da[:,:,i]*db[:,:,j],axis=-1)
        max_cross=np.maximum(max_cross,np.linalg.norm(cross,axis=-1));min_dot=np.minimum(min_dot,dot)
    degenerate=(lower_a<=1e-14)|(lower_b<=1e-14)|(min_dot<=0)|(~monotone)
    sine=np.minimum(1,max_cross/np.maximum(lower_a*lower_b,1e-30))
    angle=np.arcsin(sine);angle[degenerate]=np.pi
    position[exact]=0;angle[exact]=0;degenerate[exact]=False
    # Minimal rotation maps each original circle to the corresponding new one.
    # This is a conservative correspondence bound for the swept circle surface.
    surface=position+2*radius[:,None]*np.sin(angle*.5)
    values=(q,position.max(1),surface.max(1),angle.max(1),degenerate.any(1))
    if details:
        bad=((position>POSITION_LIMIT)|(surface>SURFACE_LIMIT)|degenerate).any(0)
        return values+({'failed_candidate_spans':np.unique(cin[bad]).tolist(),'failed_source_spans':np.unique(oin[bad]).tolist(),'fit_target':target,'fit_segments':ci,'fit_parameters':v},)
    return values

def fit_keys(p,radii,plan,protected,steps=3):
    """Small geometric least-squares fit; fixed endpoint/crown keys stay exact."""
    q=p[:,list(plan)].copy();history=[];f,m,_=q.shape
    fixed=np.array([j for j,k in enumerate(plan) if k in protected],dtype='i4')
    free=np.array([j for j in range(m) if j not in fixed],dtype='i4')
    for iteration in range(steps+1):
        q,pos,surf,ang,deg,info=candidate_and_bound(p,radii,tuple(plan),True,q)
        history.append({'iteration':iteration,'position_bound_um':float(pos.max()*1e6),'surface_bound_um':float(surf.max()*1e6),'passing_fibres':int(((pos<=POSITION_LIMIT)&(surf<=SURFACE_LIMIT)&(~deg)).sum())})
        if iteration==steps:break
        t=info['fit_parameters'];segment=info['fit_segments'];target=info['fit_target'];k=t.shape[1]
        weights=np.stack((-.5*t+t*t-.5*t**3,1-2.5*t*t+1.5*t**3,.5*t+2*t*t-1.5*t**3,-.5*t*t+.5*t**3),axis=-1)
        design=np.zeros((f,k,m),dtype='f8');ii=np.arange(f)[:,None];kk=np.arange(k)[None,:]
        for j in range(4):
            idx=np.clip(segment+j-1,0,m-1);design[ii,kk,idx[None]]+=weights[:,:,j]
        rhs=target-np.matmul(design[:,:,fixed],q[:,fixed].astype('f8'))
        mat=design[:,:,free];ata=np.matmul(mat.transpose(0,2,1),mat);atb=np.matmul(mat.transpose(0,2,1),rhs)
        lam=1e-10;ata+=np.eye(len(free))[None]*lam;atb+=q[:,free]*lam
        new=np.linalg.solve(ata,atb).astype('f4');q[:,free]=new
    assert np.array_equal(q[:,0],p[:,0]) and np.array_equal(q[:,-1],p[:,-1])
    return q,pos,surf,ang,deg,history

def reduce_batch(points,radii,menu=MENU):
    count,n,_=points.shape
    counts=np.full(count,n,dtype='i4');result=[None]*count
    pos=np.zeros(count);surf=np.zeros(count);angle=np.zeros(count);degen_attempt=np.zeros(count,dtype=bool)
    remaining=np.arange(count)
    for m in menu:
        if not len(remaining):break
        q,pp,ss,aa,dd=candidate_and_bound(points[remaining],radii[remaining],m)
        passed=(pp<=POSITION_LIMIT)&(ss<=SURFACE_LIMIT)
        degen_attempt[remaining]|=dd
        for local in np.flatnonzero(passed):
            idx=remaining[local];counts[idx]=m;result[idx]=q[local];pos[idx]=pp[local];surf[idx]=ss[local];angle[idx]=aa[local]
        remaining=remaining[~passed]
    assert len(remaining)==0 and all(x is not None for x in result)
    return result,counts,pos,surf,angle,degen_attempt

def sample_dense(points,samples=16):
    b=bezier(points.astype('f8'));t=np.linspace(0,1,samples,endpoint=False)
    out=np.einsum('sa,cfad->cfsd',basis(t),b,optimize=True).reshape(len(points),-1,3)
    return np.concatenate((out,points[:,-1:]),axis=1)
