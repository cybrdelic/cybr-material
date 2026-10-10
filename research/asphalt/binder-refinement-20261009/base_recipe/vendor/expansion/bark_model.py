"""Authored pine-type plated periderm: finite fractured flakes, not Voronoi cells.

Growth-aligned strips carry staggered, partially overlapping outer flakes.
Uneven shard edges and lifted distal lips are actual closed geometry. No species
measurement claim; this models old plated pine bark at a 180 mm inspection scale.
"""
import math,random
import numpy as np

def flakes(seed,extent):
    rng=random.Random(seed+911);out=[];x=-extent/2
    while x<extent/2:
        width=rng.uniform(.009,.024);cx=x+width/2;y=-extent/2-rng.uniform(0,.035)
        while y<extent/2:
            length=rng.uniform(.018,.048)
            if y+length>extent/2:length=extent/2-y
            if length<.004:break
            out.append(dict(cx=cx,cy=y+length/2,width=width*rng.uniform(.87,.97),length=length,
                            angle=rng.uniform(-.13,.13),lift=rng.uniform(.00015,.00125),
                            depth=rng.uniform(.00045,.0011),z=rng.uniform(.0072,.0098),
                            left=[rng.uniform(-.16,.10) for _ in range(9)],
                            right=[rng.uniform(-.10,.16) for _ in range(9)],
                            phases=[rng.uniform(0,math.tau) for _ in range(5)],tone=rng.uniform(-.025,.035)))
            y+=length-rng.uniform(.0004,.0015)
        x+=width
    return out

def sample(f,u,v):
    """Fractured piece boundary and growth-correlated relief; works on arrays."""
    u=np.asarray(u);v=np.asarray(v);phase=f['phases'];width=f['width'];length=f['length']
    left=np.interp(v.ravel(),np.linspace(0,1,9),f['left']).reshape(v.shape);right=np.interp(v.ravel(),np.linspace(0,1,9),f['right']).reshape(v.shape)
    xx=width*((u-.5)+left*(1-u)+right*u)
    # Ragged top/bottom ends terminate at unequal positions, as fractured flakes.
    end=.0006*np.sin(u*math.pi*7+phase[0])+.00035*np.sin(u*math.pi*15+phase[1])
    yy=(v-.5)*length+end*(2*v-1)
    center=.00035*np.sin(v*4.1+phase[2])+.0002*np.sin(v*11.3+phase[3])
    axial=(xx-center)/width
    # Subordinate fibers differ per flake, with broken shallow depressions.
    ridge=.00018*np.sin(axial*(38+phase[0])*math.pi+phase[1])
    ridge+=.00011*np.sin(axial*97+v*3.4+phase[2])*(.5+.5*np.sin(v*7+phase[4]))
    fissure=.00038*np.exp(-((axial-.12*np.sin(v*4+phase[0]))/.035)**2)*(np.clip(np.sin(v*8+phase[1]),0,1))
    roof=.0008*(1-(u*2-1)**2)*(.45+.55*np.sin(v*math.pi))
    lip=f['lift']*(v**6)*(1-.35*np.cos(u*9+phase[3]))
    z=f['z']+roof+ridge-fissure+lip+.0003*np.sin(v*4+u*3+phase[4])
    ca=math.cos(f['angle']);sa=math.sin(f['angle'])
    return f['cx']+xx*ca-yy*sa,f['cy']+xx*sa+yy*ca,z

def projected_maps(s,n):
    """Analytic upper envelope sampled from the same flake parameterization.

Maps approximate the upper surface only; undercuts and lifted lips require mesh.
"""
    size=s['tile_m'];x=(np.arange(n,dtype=np.float32)[None,:]+.5)/n*size-size/2
    y=(np.arange(n,dtype=np.float32)[:,None]+.5)/n*size-size/2
    # Inner periderm at fissure floors has intrinsic different reddish-brown color.
    height=np.full((n,n),.0065,np.float32);tone=np.zeros((n,n),np.float32);coverage=np.zeros((n,n),np.float32)
    for f in flakes(s['seed'],size):
        ca=math.cos(f['angle']);sa=math.sin(f['angle']);dx=x-f['cx'];dy=y-f['cy']
        xx=dx*ca+dy*sa;yy=-dx*sa+dy*ca
        v=yy/f['length']+.5;u=xx/f['width']+.5
        left=np.interp(v,np.linspace(0,1,9),f['left']);right=np.interp(v,np.linspace(0,1,9),f['right'])
        # Invert the varying broken-width shape, neglecting the small ragged ends.
        u=(u-left)/(1+right-left);mask=(u>=0)&(u<=1)&(v>=0)&(v<=1)
        _,_,z=sample(f,u,v);top=mask&(z>height)
        height=np.where(top,z,height);tone=np.where(top,f['tone'],tone);coverage=np.where(top,1,coverage)
    h=.5+(height-.008)/s['height_scale_m']
    grain=.009*np.sin(x*3900+y*21)+.005*np.sin(x*8500+y*81)
    rgb=np.stack(np.broadcast_arrays(.245+tone+grain,.12+tone*.7+grain*.65,.048+tone*.4+grain*.35),axis=-1)
    inner=np.array((.19,.082,.026),np.float32)
    rgb=rgb*coverage[...,None]+inner*(1-coverage[...,None])
    rough=.86+.05*(1-coverage)
    return rgb.astype(np.float32),np.clip(h,0,1).astype(np.float32),rough.astype(np.float32),coverage
