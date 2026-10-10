"""Mature Douglas-fir-inspired ridge/furrow hierarchy, no staggered strips.

Shape reference only: NC State Extension Pseudotsuga menziesii description and
OSU Landscape Plants bark views. No reference pixels enter this generator.
Dimensions are authored inspection-scale estimates, not botanical measurements.
"""
import math,random
import numpy as np
from PIL import Image,ImageDraw
from scipy.ndimage import gaussian_filter,distance_transform_edt
from surface_math import field,paths,clamp,color,blend,smooth

def structure(s,n):
    rng=random.Random(s['seed']+3101);extent=s['tile_m']
    major=Image.new('L',(n,n));secondary=Image.new('L',(n,n))
    draw=ImageDraw.Draw(major);sub=ImageDraw.Draw(secondary)
    # Irregular renewal positions rather than equal columns or a staggered lattice.
    positions=[];x=-.025
    while x<1.03:
        x+=rng.uniform(.060,.170);positions.append(x)
    tracks=[]
    for i,x in enumerate(positions):
        controls=[];xx=x;velocity=rng.uniform(-.07,.07)
        for j in range(25):
            velocity=.55*velocity+rng.uniform(-.022,.022)
            xx+=velocity*.16;controls.append((xx,j/24))
        tracks.append(controls)
        for j,(a,b) in enumerate(zip(controls,controls[1:])):
            width=rng.uniform(.0011,.0047)/extent*n
            draw.line([(a[0]*n,a[1]*n),(b[0]*n,b[1]*n)],fill=rng.randint(175,255),width=max(1,round(width)))
    # Offshoots start on a major fissure and arrest inside a neighboring ridge.
    for i,track in enumerate(tracks):
        for k in range(rng.randint(4,10)):
            j=rng.randrange(2,22);cx,cy=track[j];length=rng.uniform(.025,.16);direction=rng.choice((-1,1))
            points=[(cx*n,cy*n)];x=cx;y=cy
            for step in range(12):
                x+=direction*length/12*rng.uniform(.5,1.2);y+=rng.uniform(-.012,.018)
                points.append((x*n,y*n))
            sub.line(points,fill=rng.randint(140,245),width=max(1,round(rng.uniform(.0003,.0011)/extent*n)))
    a=np.asarray(major,np.float32)/255;b=np.asarray(secondary,np.float32)/255
    # Broken shoulders and variable depth, not flat flake surfaces with repeat ribs.
    primary=gaussian_filter(a,max(.7,n*.00038/extent))
    small=gaussian_filter(b,max(.35,n*.00015/extent))
    d=distance_transform_edt(a<.15)*extent/n
    shoulders=np.exp(-(d/(.0015+.0005*field(n,(19,11),s['seed']+3110)))**1.25)
    ridge=field(n,(31,63),s['seed']+3111,2)
    broad=field(n,(6,12),s['seed']+3112,2)
    grain=field(n,(130,330),s['seed']+3113)
    interrupted=paths(n,s['seed']+3114,2400,width=(.00010,.00035),length=(.001,.045),angle=math.pi/2,jitter=.31,bend=.12)
    h=.67+.067*ridge+.05*broad+.012*grain-.46*primary-.17*shoulders-.20*small-.038*interrupted
    # Outer cork gray/brown vs dark inner periderm; no bright smooth ochre gaps.
    rgb=color((.245,.208,.166),.023*broad+.016*ridge+.005*grain)
    rgb=blend(rgb,(.093,.070,.048),clamp(primary*.6+small*.22))
    rough=.87+.04*primary+.018*ridge
    return clamp(rgb).astype(np.float32),clamp(h).astype(np.float32),clamp(rough).astype(np.float32),clamp(1-primary*.36-small*.18),primary
