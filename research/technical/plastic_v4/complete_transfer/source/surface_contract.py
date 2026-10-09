"""Unchanged original form + retained native field, with analytic derivatives.
This is a reconstruction contract. It does not recover a historical scene binary.
"""
from pathlib import Path
import hashlib, sys
import numpy as np
from PIL import Image
ROOT=Path(__file__).resolve().parents[1]
EXPERIMENTS=ROOT.parent
MACRO=EXPERIMENTS/'plastic_curved_transfer/source/analytic_surface.py'
MAP=EXPERIMENTS/'plastic_native_relief_transfer/maps/production_retained_height_native4096.png'
MAP_SHA='78ad8b37255340889593dcf47180913d611ce0124e756f3dcbaac187d9eb2aee'
MACRO_SHA='31280c3ec450d251fa5da8a77f4bf4485fdccfbe0fe805571f8052e67e29698f'
assert hashlib.sha256(MACRO.read_bytes()).hexdigest()==MACRO_SHA
assert hashlib.sha256(MAP.read_bytes()).hexdigest()==MAP_SHA
sys.path.insert(0,str(MACRO.parent))
from analytic_surface import base_and_direction
HEIGHT=(np.flipud(np.asarray(Image.open(MAP))).astype('f8')/65535-.5)*80e-6
PITCH=.08/4096
TANGENTS=np.array([[[1,0,0],[0,1,0]],[[0,1,0],[1,0,0]],[[0,1,0],[0,0,1]],[[0,0,1],[0,1,0]],[[1,0,0],[0,0,1]],[[0,0,1],[1,0,0]]],dtype='f8')

def native_field(xy):
    raw=(np.asarray(xy,dtype='f8')/.08+.5)*4096-.5
    ij=np.clip(raw,0,4095);lo=np.floor(ij).astype('i4');hi=np.minimum(lo+1,4095);f=ij-lo
    h00=HEIGHT[lo[...,1],lo[...,0]];h10=HEIGHT[lo[...,1],hi[...,0]]
    h01=HEIGHT[hi[...,1],lo[...,0]];h11=HEIGHT[hi[...,1],hi[...,0]]
    x,y=f[...,0],f[...,1]
    h=(1-x)*(1-y)*h00+x*(1-y)*h10+(1-x)*y*h01+x*y*h11
    hx=((1-y)*(h10-h00)+y*(h11-h01))/PITCH
    hy=((1-x)*(h01-h00)+x*(h11-h10))/PITCH
    hx=np.where((raw[...,0]<0)|(raw[...,0]>=4095),0,hx)
    hy=np.where((raw[...,1]<0)|(raw[...,1]>=4095),0,hy)
    return h,np.stack((hx,hy),-1)

def macro_differential(q,t):
    """Directional derivative of retained b,d; knot values use interior clip branch."""
    q=np.asarray(q,dtype='f8');t=np.broadcast_to(np.asarray(t,dtype='f8'),q.shape)
    half=np.array([.037,.037,.001]);c=np.clip(q,-half,half)
    dc=t*(np.abs(q)<half);v=q-c;dv=t-dc;r=np.linalg.norm(v,axis=-1,keepdims=True)
    n=v/r;dn=(dv-n*np.sum(n*dv,axis=-1,keepdims=True))/r
    s=c+.003*n+[0,0,.004];ds=dc+.003*dn;z=s[...,2]-.004
    draft=np.tan(np.deg2rad(.8))/.04;a=1-draft*np.abs(z);da=-draft*np.sign(z)*ds[...,2]
    e=np.exp(-.5*(z/.000027)**2);fac=np.maximum(0,1-n[...,2]**2)
    seam=.000022*e*fac;dseam=.000022*e*((-z/.000027**2)*ds[...,2]*fac-2*n[...,2]*dn[...,2])
    db=ds.copy();db[...,:2]=ds[...,:2]*a[...,None]+s[...,:2]*da[...,None]+dn[...,:2]*seam[...,None]+n[...,:2]*dseam[...,None]
    p=s.copy();p[...,:2]=s[...,:2]*a[...,None]+n[...,:2]*seam[...,None]
    mask=n[...,2]<-.999
    for cx in [-.023,.023]:
        for cy in [-.023,.023]:
            delta=p[...,:2]-[cx,cy];dist=np.linalg.norm(delta,axis=-1)
            dr=np.divide(np.sum(delta*db[...,:2],axis=-1),dist,out=np.zeros_like(dist),where=dist>0)
            ring=.000006*np.exp(-.5*((dist-.0022)/.000095)**2)
            exponent=np.minimum(50,(dist-.0022)/.00006);ex=np.exp(exponent)
            derivative=-ring*(dist-.0022)/.000095**2-.000002*ex/(1+ex)**2/.00006*((dist-.0022)/.00006<50)
            db[...,2]+=mask*derivative*dr
    nz=np.maximum(n[...,2],0);dd=dn*nz[...,None]**6+n*(6*nz**5*dn[...,2])[...,None]
    return db,dd

def surface(q):
    b,d=base_and_direction(q);h,_=native_field(b[...,:2]);return b+d*h[...,None]

def normal(q,chart):
    q=np.asarray(q,dtype='f8');chart=np.broadcast_to(np.asarray(chart,dtype='i4'),q.shape[:-1]);b,d=base_and_direction(q);h,g=native_field(b[...,:2]);tt=[]
    for i in [0,1]:
        db,dd=macro_differential(q,TANGENTS[chart,i]);tt.append(db+dd*h[...,None]+d*np.sum(g*db[...,:2],axis=-1)[...,None])
    n=np.cross(tt[0],tt[1]);return n/np.linalg.norm(n,axis=-1,keepdims=True)
