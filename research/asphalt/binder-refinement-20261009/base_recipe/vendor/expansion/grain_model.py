"""Graded individual mineral grains with finite crowns, not polygon-cell fills."""
import math
import numpy as np
from surface_math import smooth

def grains(n,seed,count,radii,angularity=.18):
    rng=np.random.default_rng(seed);mask=np.zeros((n,n),np.float32);height=mask.copy();tone=mask.copy()
    for i in range(count):
        cx,cy=rng.uniform(0,n,2);r=math.exp(rng.uniform(math.log(radii[0]),math.log(radii[1])))*n
        aspect=math.exp(rng.uniform(-.4,.4));rx=r*math.sqrt(aspect);ry=r/math.sqrt(aspect)
        extent=max(rx,ry)*1.35+1;ix=np.arange(int(cx-extent),int(cx+extent+2));iy=np.arange(int(cy-extent),int(cy+extent+2))
        dx=ix[None,:]-cx;dy=iy[:,None]-cy;angle=rng.uniform(0,math.tau);ct=math.cos(angle);st=math.sin(angle)
        u=(ct*dx+st*dy)/max(.55,rx);v=(-st*dx+ct*dy)/max(.55,ry)
        corners=rng.integers(6,11);angles=np.linspace(0,math.tau,corners,endpoint=False);outline=rng.uniform(1-angularity,1+angularity,corners)
        phi=np.mod(np.arctan2(v,u),math.tau);contour=np.interp(phi,np.r_[angles[-1]-math.tau,angles,angles[0]+math.tau],np.r_[outline[-1],outline,outline[0]])
        rr=np.sqrt(u*u+v*v)/contour;aa=min(1,rx/.55)*min(1,ry/.55)
        coverage=(1-smooth(rr,.84,1.08))*aa
        crown=np.sqrt(np.clip(1-(rr/1.08)**2,0,1))*coverage*rng.uniform(.45,1)
        tone_value=rng.uniform(-1,1);target=np.ix_(iy%n,ix%n);wins=crown>height[target]
        tone[target]=np.where(wins,tone_value,tone[target]);height[target]=np.maximum(height[target],crown)
        mask[target]=np.maximum(mask[target],coverage)
    return mask,height,tone
