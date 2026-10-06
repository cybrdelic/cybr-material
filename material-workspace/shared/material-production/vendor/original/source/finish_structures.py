"""Finish-specific structures in physical scale. No photographic inputs.

Quiet finish maps describe fabrication, not arbitrary widespread damage.
Mesh-specific bends, seams and exposed edges are authored by the scene builder.
"""
import numpy as np
from scipy.ndimage import gaussian_filter
from surface_math import field, smooth, clamp, color, blend, paths, TAU, voronoi
from physical_features import deposits, trowel_layers


def marble(x,y,n):
    # One 2.8 m quarry slab: sparse primary shear seams, subsidiary seams,
    # mineral intergrowth and feathered fracture ends. No periodic ribbons.
    geo=field(n,(12,9),911,3);grain=field(n,180,912,2)
    vein=np.zeros((n,n),np.float32);gold=vein.copy();mist=vein.copy()
    for i,(origin,slope,width) in enumerate(((.13,.67,.0018),(.62,-.31,.0026),(.91,.42,.0010))):
        track=origin+slope*(y-.5)+.028*np.sin(y*9+i*2)+.011*np.sin(y*31+i)
        distance=abs(x-track-.0045*geo-.0012*grain)
        w=width*(.45+1.2*smooth(field(n,(41,23),916+i,2),-.5,.4))
        vein=np.maximum(vein,1-smooth(distance,w*.15,w))
        gold=np.maximum(gold,np.exp(-((distance-w*.8)/(w*.21+.00006))**2))
        mist=np.maximum(mist,np.exp(-(distance/(w*5+.002))**2)*.13)
        for k in range(7):
            start=.09+k*.13;direction=(-1 if k%2 else 1)
            branch=track+direction*(y-start)*(.18+.08*np.sin(k))
            active=smooth(y,start,start+.018)*(1-smooth(y,start+.07,start+.24))
            vein=np.maximum(vein,(1-smooth(abs(x-branch-.003*geo),.0002,.0008))*active)
    cloudy=field(n,(21,15),930,3)
    rgb=color((.82,.809,.776),.018*cloudy+.003*grain)
    rgb=blend(rgb,(.54,.55,.535),mist)
    rgb=blend(rgb,(.36,.39,.385),vein*(.44+.18*smooth(grain,-.6,.5)))
    ochre=gold*smooth(field(n,13,931),.12,.55)
    rgb=blend(rgb,(.55,.397,.22),ochre*.40)
    h=.5+.003*grain-.010*vein
    rough=.255+.034*vein+.012*cloudy
    return rgb,h,rough,np.ones((n,n),np.float32),np.zeros((n,n),np.float32),ochre*.1


def travertine(x,y,n):
    flow=.007*field(n,(7,12),941,2)
    phase=(y+flow)*24
    strata=field(n,(56,5),942,2)+.14*np.sin(phase*TAU)
    beds=.18+.82*smooth(strata,-.1,.65)
    pores,depth,_=deposits(n,943,4200,(.0007,.0048),aspect=(.6,2.9),angle=.12,spread=.48,location_mask=beds,angularity=.18,floor=.10)
    small,sd,_=deposits(n,944,20000,(.00010,.0007),aspect=(.7,1.9),location_mask=.25+.75*beds,angularity=.12)
    rgb=color((.731,.664,.548),.028*strata+.008*field(n,(130,26),945,2))
    # Color does not outline the cavities. Relief and correlated roughness do.
    h=.60+.007*strata-.46*depth-.035*sd
    rough=.54+.12*pores+.025*small+.018*strata
    return rgb,h,rough,clamp(1-.48*depth-.12*sd),np.zeros((n,n),np.float32),pores*.08


def metal(x,y,n,steel=False):
    seed=960 if steel else 950
    tooling=field(n,(8,1300),seed,2)
    fine=field(n,(26,1900),seed+1)
    # Scores are rare and narrow. Oxide is an intentionally intact finish.
    scratch=paths(n,seed+2,42,width=(.00005,.00012),length=(.005,.04),angle=np.pi/2,jitter=.06,bend=.002)
    rgb=color((.132,.147,.151) if steel else (.77,.646,.405),.0035*tooling+.0015*fine)
    h=.50+.008*tooling+.002*fine-.015*scratch
    rough=(.38 if steel else .30)+.016*tooling+.022*scratch
    metallic=np.full((n,n),.68 if steel else 1,np.float32)
    return rgb,h,rough,np.ones((n,n),np.float32),metallic,scratch*.18


def porcelain(x,y,n):
    peel=field(n,180,971,2);flow=field(n,(6,24),972,2)
    rgb=color((.835,.817,.761),.003*flow)
    h=.5+.005*peel+.002*flow
    rough=.18+.014*peel+.006*flow
    return rgb,h,rough,np.ones((n,n),np.float32),np.zeros((n,n),np.float32),np.zeros((n,n),np.float32)


def leather(x,y,n):
    # 0.4–1.0 mm irregular grain clusters; height range under 55 microns.
    gap,dist,tone=voronoi(n,260,981)
    body=smooth(gap,.003,.11)
    fine=field(n,630,982,2)
    follicles,fd,_=deposits(n,983,18000,(.00012,.00024),aspect=(.7,1.4),angularity=.06)
    dye=field(n,(13,17),984,2)
    rgb=color((.378,.205,.104),.009*dye+.003*tone)
    rgb=blend(rgb,(.28,.135,.064),follicles*.035)
    h=.49+.044*body+.007*fine-.018*fd
    rough=.50+.025*(1-body)+.009*dye
    return rgb,h,rough,clamp(1-.035*fd),np.zeros((n,n),np.float32),np.zeros((n,n),np.float32)


def plaster(x,y,n):
    layers,lips,polish=trowel_layers(n,991,19)
    skim=field(n,(9,11),992,2);fine=field(n,380,993,2)
    rgb=color((.774,.755,.713),.017*layers+.009*skim+.0015*fine)
    h=.46+.16*layers+.043*lips+.005*fine
    rough=.75-.08*polish+.018*skim
    return rgb,h,rough,np.ones((n,n),np.float32),np.zeros((n,n),np.float32),lips*.1
