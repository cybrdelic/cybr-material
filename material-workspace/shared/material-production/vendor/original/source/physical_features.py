"""Finite structures with named causes, physical sizes, and toroidal boundaries.

Random numbers position and vary shapes. They never become pixel color or bump.
"""
import math
import numpy as np
from PIL import Image, ImageDraw
from scipy.ndimage import gaussian_filter
from surface_math import clamp, smooth, TAU


def deposits(n, seed, count, radii, aspect=(.65, 1.4), angle=0., spread=math.pi,
             location_mask=None, angularity=.22, floor=.5):
    """Torn mineral grains / cavities with different depths and irregular rims.

    Radii are fractions of a tile. Sizes follow a log distribution; subpixel
    shapes retain area coverage instead of becoming oversized single pixels.
    """
    rng=np.random.default_rng(seed)
    coverage=np.zeros((n,n),np.float32);depth=coverage.copy();tone=coverage.copy()
    for _ in range(count):
        cx,cy=rng.uniform(0,n,2)
        if location_mask is not None and rng.random()>location_mask[int(cy)%n,int(cx)%n]:
            continue
        r=np.exp(rng.uniform(np.log(radii[0]),np.log(radii[1])))*n
        elong=np.exp(rng.uniform(np.log(aspect[0]),np.log(aspect[1])))
        rx,ry=r*math.sqrt(elong),r/math.sqrt(elong)
        theta=angle+rng.uniform(-spread,spread);ct,st=math.cos(theta),math.sin(theta)
        extent=max(rx,ry)*1.45+2
        ix=np.arange(int(cx-extent),int(cx+extent+2));iy=np.arange(int(cy-extent),int(cy+extent+2))
        px=ix[None,:]-cx;py=iy[:,None]-cy
        u=(ct*px+st*py)/max(rx,.60);v=(-st*px+ct*py)/max(ry,.60)
        # A variable polygonal perimeter avoids the shared three-lobed
        # flower shape produced by a strong fixed angular sine harmonic.
        corners=int(rng.integers(7,15))
        angles=np.arange(corners)*TAU/corners+rng.uniform(-.20,.20,corners)*TAU/corners
        sizes=rng.uniform(1-angularity,1+angularity,corners)
        angles=np.r_[angles[-1]-TAU,angles,angles[0]+TAU]
        sizes=np.r_[sizes[-1],sizes,sizes[0]]
        a=np.mod(np.arctan2(v,u),TAU)
        contour=np.interp(a,angles,sizes)
        radial=np.sqrt(u*u+v*v)/contour
        aa=min(1.,rx/.60)*min(1.,ry/.60)
        mask=(1-smooth(radial,.75,1.12))*aa
        # Torn cavities have a shallowly irregular floor and a steeper rim,
        # rather than several triangular ramps meeting at a central point.
        bottom=.91+.06*np.sin(u*rng.uniform(2,5)+rng.uniform(0,TAU))+.03*np.sin(v*4)
        profile=(1-smooth(radial,.63+.26*floor,1.03))*bottom*rng.uniform(.30,1.0)*aa
        target=np.ix_(iy%n,ix%n)
        coverage[target]=np.maximum(coverage[target],mask)
        depth[target]=np.maximum(depth[target],profile)
        tone[target]=np.where(mask>coverage[target]*.95,rng.uniform(-1,1)*mask,tone[target])
    return coverage,depth,tone


def fracture_network(n,seed,primaries=18,branches=650,width=.00018):
    """Growing curved shrinkage cracks, arrested at existing cracks (T joints).

    This is deliberately not a Voronoi boundary texture. A low-resolution
    topology grid decides contacts; native-resolution strokes retain detail.
    """
    rng=np.random.default_rng(seed);g=768
    occupied=np.zeros((g,g),bool)
    im=Image.new('L',(n,n));draw=ImageDraw.Draw(im)
    def grow(start,theta,length):
        points=[np.asarray(start,dtype=float)];angle=theta;curvature=rng.normal(0,.006)
        for k in range(max(2,int(length*g/2))):
            curvature=.92*curvature+rng.normal(0,.007)
            angle+=curvature
            pos=points[-1]+2*np.array([math.cos(angle),math.sin(angle)])
            xx,yy=(np.rint(pos).astype(int)%g)
            points.append(pos)
            if k>3 and occupied[yy,xx]:break
        return points
    for i in range(primaries+branches):
        start=rng.uniform(0,g,2);xx,yy=np.rint(start).astype(int)%g
        if occupied[yy,xx]:continue
        theta=rng.uniform(0,TAU)
        length=rng.uniform(.35,.8) if i<primaries else rng.uniform(.018,.13)
        left=grow(start,theta,length);right=grow(start,theta+math.pi,length)
        points=list(reversed(left))+right[1:]
        for p in points:
            px,py=np.rint(p).astype(int)%g
            for dx,dy in ((0,0),(1,0),(-1,0),(0,1),(0,-1)):
                occupied[(py+dy)%g,(px+dx)%g]=True
        pts=np.asarray(points)*n/g
        w=width*n*rng.uniform(.55,1.35)*(1 if i<primaries else .62)
        tone=int(255*min(1,w)*rng.uniform(.6,1))
        for ox in (-n,0,n):
            for oy in (-n,0,n):draw.line(list(zip(pts[:,0]+ox,pts[:,1]+oy)),fill=tone,width=max(1,round(w)))
    return gaussian_filter(np.asarray(im,dtype=np.float32)/255,.40,mode='wrap')


def trowel_layers(n,seed,count=36):
    """Overlapping tapered lime applications with a lip on the leading side."""
    rng=np.random.default_rng(seed)
    layers=np.zeros((n,n),np.float32);lips=layers.copy();polish=layers.copy()
    for _ in range(count):
        cx,cy=rng.uniform(0,n,2);length=rng.uniform(.18,.52)*n;width=rng.uniform(.035,.11)*n
        angle=rng.uniform(-.42,.48);ct,st=math.cos(angle),math.sin(angle)
        extent=(length+width)*.65+4
        ix=np.arange(int(cx-extent),int(cx+extent+2));iy=np.arange(int(cy-extent),int(cy+extent+2))
        xx=ix[None,:]-cx;yy=iy[:,None]-cy
        u=(ct*xx+st*yy)/length;v=(-st*xx+ct*yy)/width
        v+=.11*np.sin(u*4+rng.uniform(0,TAU));edge=.5-.10*u
        ends=(1-smooth(abs(u),.40,.52))
        mask=(1-smooth(abs(v),edge-.075,edge+.03))*ends
        lip=np.exp(-((v-edge)/.023)**2)*ends*(.7+.3*np.cos(u*7))
        target=np.ix_(iy%n,ix%n);amount=rng.uniform(.15,1)
        layers[target]=np.maximum(layers[target],mask*amount)
        lips[target]=np.maximum(lips[target],lip*amount)
        polish[target]=np.maximum(polish[target],mask*(.5+.5*amount))
    return layers,lips,polish
