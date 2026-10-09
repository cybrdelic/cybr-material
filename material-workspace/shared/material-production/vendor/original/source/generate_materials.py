#!/usr/bin/env python3
"""CYBR MATERIAL 3: physically sized structures and cause-specific wear.

Original deterministic procedural authorship. No photographs, baked lights,
independent pixel noise, image upscaling, or claim of measured scan fidelity.
"""
from __future__ import annotations
import argparse, gc, json
from pathlib import Path
import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter
from surface_math import clamp,smooth,field,color,blend,paths,warp_map,TAU
from physical_features import deposits,fracture_network,trowel_layers

ROOT=Path(__file__).resolve().parents[1]
SPECS=[
 dict(id='01_calacatta_oro',name='Calacatta Oro / Weathered',family='Stone',tile_m=.8,height_scale_m=.00055,roughness=.34,ior=1.52,metallic=0,coat=.015,coat_roughness=.38,preview_diameter_m=.24,description='Ivory marble with a hierarchy of ragged mineral seams, fractured inclusions, crystalline intergrowth and shallow abrasion.'),
 dict(id='02_roman_travertine',name='Roman Travertine / Open Pore',family='Stone',tile_m=.65,height_scale_m=.006,roughness=.62,ior=1.5,metallic=0,preview_diameter_m=.24,description='Vein-cut travertine with bed-clustered torn cavities, varied cavity floors, mineral laminae and fine calcite pores.'),
 dict(id='03_american_walnut',name='American Walnut / Timeworn',family='Wood',tile_m=.55,tile_y_m=.55,height_scale_m=.00042,roughness=.47,ior=1.5,metallic=0,coat=.025,coat_roughness=.52,anisotropy=.12,preview_diameter_m=.18,tiling='finite-cut',description='Procedural flat-sawn walnut sampled from a radial growth volume, with nonuniform annual rings, diffuse lumen windows, fine rays and shallow worn oil finish.'),
 dict(id='04_fumed_oak',name='Fumed Oak / Open Grain',family='Wood',tile_m=.45,tile_y_m=.45,height_scale_m=.00045,roughness=.55,ior=1.5,metallic=0,coat=.015,coat_roughness=.58,anisotropy=.12,preview_diameter_m=.16,tiling='finite-cut',description='Procedural rift-cut fumed oak with growth-volume earlywood vessels, finite medullary ray plates, shallow lumen relief and worn low-gloss finish.'),
 dict(id='05_champagne_brass',name='Champagne Brass / Aged',family='Metal',tile_m=.2,height_scale_m=.00020,roughness=.37,ior=1.5,metallic=1,anisotropy=.19,preview_diameter_m=.10,description='Aged brass with subdued tooling, localized abrasion bundles, granular tarnish fronts and small dielectric verdigris deposits.'),
 dict(id='06_blackened_steel',name='Blackened Steel / Worn',family='Metal',tile_m=.25,height_scale_m=.00055,roughness=.49,ior=1.5,metallic=0,anisotropy=.14,preview_diameter_m=.12,description='Black oxide steel with directionally rubbed exposure, terraced rust crusts, corrosion pits and short abrasive score groups.'),
 dict(id='07_bone_porcelain',name='Bone Porcelain / Crazed',family='Ceramic',tile_m=.24,height_scale_m=.00045,roughness=.24,ior=1.49,metallic=0,coat=.12,coat_roughness=.24,preview_diameter_m=.10,description='Ivory porcelain with a growing arrested crack network, subdued glaze ripples, pinholes and torn chips exposing the ceramic body.'),
 dict(id='08_saddle_leather',name='Saddle Leather / Used',family='Leather',tile_m=.18,height_scale_m=.00075,roughness=.50,ior=1.48,metallic=0,coat=.015,coat_roughness=.48,preview_diameter_m=.065,description='Cognac hide with soft stretched skin creases, fine follicles, compressed folds, burnished areas and grouped finish scuffs.'),
 dict(id='09_lime_plaster',name='Lime Plaster / Weathered',family='Plaster',tile_m=.65,height_scale_m=.004,roughness=.76,ior=1.48,metallic=0,preview_diameter_m=.22,description='Layered lime with overlapping trowel passes, directional leading lips, exposed angular aggregate, delamination and fine pinholes.'),
 dict(id='10_natural_linen',name='Natural Linen / Frayed',family='Textile',tile_m=.085,height_scale_m=.00065,roughness=.80,ior=1.46,metallic=0,sheen=.32,sheen_roughness=.66,preview_diameter_m=.028,description='Fine plain flax weave with 120 variable yarns per tile, continuous over-under bending, fibrils, attached slubs and real transparent inter-yarn openings.'),
]
_wood_profiles=json.loads((ROOT/'source/wood_profiles.json').read_text())
for _spec,_profile in zip(SPECS[2:4],(_wood_profiles['american_walnut'],_wood_profiles['fumed_oak'])):
    _spec.update(tile_m=_profile['extent_m'],tile_y_m=_profile['extent_m'],height_scale_m=_profile['height_scale_m'],roughness=_profile['roughness'])


def marble(x,y,n):
    geology=field(n,(7,5),112,3);breccia=field(n,58,113,2)
    vein=np.zeros((n,n),np.float32);rim=vein.copy();halo=vein.copy()
    # Integer winding and periodic oscillations close every geological seam.
    for i,(base,winding,width) in enumerate(((-.10,1,.008),(.50,-1,.013),(.82,1,.0045))):
        track=base+winding*y+.073*np.sin(TAU*y+i)+.034*np.sin(TAU*y*3+i*2)
        d=np.abs((x-track+.013*geology+.003*breccia+.5)%1-.5)
        w=width*(.65+1.75*smooth(field(n,(17,11),121+i,2),-.55,.55))
        fractured=w*(1+.36*breccia)
        vein=np.maximum(vein,1-smooth(d,fractured*.32,fractured))
        rim=np.maximum(rim,np.exp(-((d-fractured*.90)/(w*.19+.00035))**2))
        halo=np.maximum(halo,(1-smooth(d,w*1.0,w*3.6))*smooth(breccia,-.40,.27))
    branches=paths(n,144,390,width=(.00015,.0010),length=(.016,.18),angle=.82,jitter=.65,bend=.024)
    branches*=.15+.85*smooth(gaussian_filter(vein,n*.009,mode='wrap'),.015,.13)
    crystal,_,ct=deposits(n,145,11000,(.00018,.0022),aspect=(.4,2.7),location_mask=.12+.88*vein,angularity=.35)
    fragments,_,ft=deposits(n,153,12000,(.0007,.0038),aspect=(.4,2.3),location_mask=halo,angularity=.42)
    calcite,_,_=deposits(n,154,36000,(.00016,.0011),aspect=(.6,1.8),angularity=.36)
    fracture=fracture_network(n,146,primaries=9,branches=280,width=.00012)*smooth(vein,.08,.55)
    mineral=smooth(field(n,81,147,2),-.15,.4)*vein
    rgb=color((.814,.800,.759),.027*geology+.010*field(n,73,148))
    rgb=blend(rgb,(.485,.491,.470),vein*(.55+.30*smooth(breccia,-.4,.5)))
    rgb=blend(rgb,(.738,.711,.643),mineral*.47+crystal*.24)
    rgb=blend(rgb,(.585,.585,.540),fragments*.36)
    rgb=blend(rgb,(.842,.828,.789),calcite*.075)
    rgb=blend(rgb,(.388,.397,.380),branches*.43+fracture*.23)
    ochre=rim*smooth(field(n,9,149),.04,.40)
    rgb=blend(rgb,(.561,.410,.221),ochre*.61)
    pits,pitdepth,_=deposits(n,150,3100,(.00011,.00065),angularity=.3)
    abrasion=paths(n,151,58,width=(.00012,.00038),length=(.005,.035),angle=.18,jitter=.4,bend=.006)
    h=.55+.009*field(n,93,152)-.16*pitdepth-.028*fracture-.018*branches-.018*crystal-.025*abrasion+.004*calcite
    rough=.31+.11*vein+.052*crystal+.032*fragments+.014*calcite+.10*pits+.11*abrasion+.016*geology
    return rgb,h,rough,1-.28*pitdepth,np.zeros((n,n),np.float32),clamp(pits*.55+abrasion*.7+fracture*.3)


def travertine(x,y,n):
    flow=.016*field(n,(9,14),201,2)+.009*np.sin(TAU*x*2)
    bedding=29*(y+flow);frac=bedding%1;idx=np.floor(bedding).astype(int)%29
    rng=np.random.default_rng(202);tones=rng.uniform(-1,1,29).astype(np.float32)
    bands=tones[idx]*(1-smooth(frac,.64,1))+tones[(idx+1)%29]*smooth(frac,.64,1)
    porous=.12+.88*smooth(field(n,(48,16),203,2)+.58*bands,-.30,.60)
    cavities,depth,ct=deposits(n,204,2600,(.0011,.010),aspect=(2,12),angle=0,spread=.055,location_mask=porous,angularity=.36,floor=.26)
    tiny,tinydepth,_=deposits(n,205,19000,(.00013,.0012),aspect=(.5,2.4),location_mask=.30+.70*porous,angularity=.28)
    mineral,_,mt=deposits(n,206,4900,(.0002,.002),aspect=(.6,2.7),angularity=.4)
    lamina=np.exp(-((frac-.13)/.042)**2)*(1-cavities)
    rgb=color((.735,.667,.548),.045*bands+.016*field(n,(100,35),207,2))
    # Cavity shadows come from geometry. Color only describes mineral residues.
    residue=cavities*smooth(ct,-.08,.62)
    rgb=blend(rgb,(.536,.435,.283),residue*.44+tiny*.055)
    rgb=blend(rgb,(.817,.761,.645),mineral*.24+lamina*.11)
    rgb=blend(rgb,(.648,.558,.409),smooth(mt,.1,.65)*mineral*.22)
    h=.66+.021*bands+.014*lamina-.57*depth-.083*tinydepth+.008*mineral
    rough=.58+.18*cavities+.065*tiny+.022*bands+.025*mineral
    return rgb,h,rough,clamp(1-.66*depth-.17*tinydepth),np.zeros((n,n),np.float32),clamp(cavities+tiny*.23)


def wood(x,y,n,oak=False):
    from wood_anatomy import wood as procedural_wood
    return procedural_wood(x,y,n,oak)


def metal(x,y,n,steel=False):
    seed=460 if steel else 420
    # Micrometer tooling relief, rather than high-contrast corrugated stripes.
    tooling=field(n,(28,780),seed+1,2)
    brush=paths(n,seed+2,1600,width=(.00010,.00028),length=(.015,.20),angle=np.pi/2,jitter=.012,bend=.003)
    contact=smooth(field(n,(6,9),seed+3,2),.02,.51)
    scores=paths(n,seed+4,245,width=(.00010,.00042),length=(.008,.09),angle=.32,jitter=.28,bend=.004)*contact
    scuffs=paths(n,seed+5,460,width=(.00013,.00058),length=(.003,.032),angle=.31,jitter=.15,bend=.002)*contact
    geography=field(n,8,seed+6,3)+.24*field(n,53,seed+7,2)
    fronts=smooth(geography,.08,.40)
    particles,pd,pt=deposits(n,seed+8,46000,(.00018,.0022),aspect=(.4,2.7),location_mask=.07+.93*fronts,angularity=.40)
    grain=smooth(field(n,195,seed+9),-.20,.43)
    crust=fronts*(.37+.40*particles+.23*grain)
    pits,pitdepth,_=deposits(n,seed+10,3800,(.0002,.0024),location_mask=.10+.90*fronts,angularity=.37)
    dents,dd,_=deposits(n,seed+11,42,(.0016,.006),aspect=(.8,1.5),angularity=.13)
    if steel:
        rubbed=contact*(.13+.87*smooth(field(n,(24,110),seed+12,2),-.01,.40))
        exposed=clamp(rubbed*.56+scores*.42+scuffs*.53)*(1-fronts*.6)
        rust=smooth(geography,.28,.54)*(.45+.40*particles+.15*grain)*(1-exposed)
        rgb=color((.159,.178,.183),.006*tooling+.009*field(n,18,seed+13))
        rgb=blend(rgb,(.458,.489,.499),exposed*.86)
        rgb=blend(rgb,(.325,.161,.072),rust*.90)
        rgb=blend(rgb,(.469,.282,.143),rust*particles*smooth(pt,-.10,.65)*.57)
        rgb=blend(rgb,(.169,.101,.056),rust*pits*.41)
        metallic=exposed*(1-rust)
        rough=.47+.23*rust-.19*exposed+.045*particles*rust+.012*tooling
        h=.50+.009*tooling-.012*brush-.09*scores-.07*scuffs-.18*dd-.14*pitdepth+.28*rust*(.34+.66*pd)
        wear=clamp(exposed+rust+dents*.4)
    else:
        patina=smooth(geography,.42,.64)*(.15+.85*particles)*(1-contact*.5)
        rgb=color((.790,.682,.463),.005*tooling+.008*field(n,18,seed+13))
        rgb=blend(rgb,(.425,.299,.134),crust*.81)
        rgb=blend(rgb,(.245,.363,.288),patina*.86)
        rgb=blend(rgb,(.647,.464,.251),fronts*particles*smooth(pt,.1,.60)*.20)
        rgb=blend(rgb,(.825,.724,.504),scores*.16+scuffs*.13)
        metallic=clamp(1-smooth(crust,.12,.47)*.98-patina*.70)
        rough=.34+.30*crust+.12*patina-.062*contact+.085*scuffs+.012*tooling
        h=.52+.011*tooling-.018*brush-.13*scores-.09*scuffs-.25*dd-.16*pitdepth+.19*crust*pd+.06*patina
        wear=clamp(crust+scores*.55+scuffs*.45+dents)
    return rgb,h,rough,clamp(1-.26*dd-.32*pitdepth),metallic,wear


def porcelain(x,y,n):
    cracks=fracture_network(n,501,primaries=22,branches=1650,width=.00016)
    pinholes,pd,_=deposits(n,502,1700,(.00016,.00062),angularity=.22)
    chips,cd,ct=deposits(n,503,49,(.001,.0063),aspect=(.6,2.2),angularity=.40,floor=.18)
    burnish=smooth(field(n,8,504,2),.05,.48)
    scuff=paths(n,505,180,width=(.00010,.00032),length=(.004,.020),angle=.6,jitter=.16,bend=.007)*burnish
    orange_peel=field(n,145,506,2)
    rgb=color((.825,.809,.756),.009*field(n,16,507)+.003*orange_peel)
    rgb=blend(rgb,(.510,.444,.338),cracks*.43+pinholes*.06)
    rgb=blend(rgb,(.679,.592,.451),chips*.81)
    rgb=blend(rgb,(.857,.840,.781),scuff*.11)
    h=.55+.007*orange_peel-.038*cracks-.15*pd-.58*cd-.012*scuff
    rough=.22+.13*cracks+.26*chips-.065*burnish+.14*scuff+.018*pinholes+.013*orange_peel
    return rgb,h,rough,clamp(1-.16*pd-.37*cd),np.zeros((n,n),np.float32),clamp(chips+scuff+burnish*.35)


def leather(x,y,n):
    stress=field(n,(7,9),601,3)
    u=x+.008*field(n,18,602,2)+.003*stress;v=y+.006*field(n,17,603,2)
    skin=warp_map(field(n,(213,167),604,2),u,v)
    fine=warp_map(field(n,387,605),u,v)
    # Smooth stretched skin contours replace straight Voronoi cell borders.
    crease=np.exp(-((skin+.12*fine)/.105)**2)
    softskin=gaussian_filter(np.tanh(skin*2),max(.45,n*.00019),mode='wrap')
    folds=paths(n,606,19,width=(.0007,.0024),length=(.10,.45),angle=.72,jitter=.30,bend=.048)
    folds=gaussian_filter(folds,max(.6,n*.0005),mode='wrap')
    follicles,fd,_=deposits(n,607,15200,(.00015,.00038),aspect=(.65,1.75),angle=.8,spread=.2,angularity=.08)
    burnish=smooth(field(n,(9,7),608,2),.02,.5)*(1-crease*.27)
    rub=smooth(field(n,(13,9),609,2),.18,.54)
    scuffs=paths(n,610,470,width=(.00013,.00049),length=(.002,.016),angle=.19,jitter=.23,bend=.005)*rub
    dry=clamp(scuffs*.82+rub*.12)
    rgb=color((.441,.256,.145),.025*stress+.007*softskin)
    rgb=blend(rgb,(.311,.141,.061),crease*.10+folds*.21+follicles*.11)
    rgb=blend(rgb,(.471,.281,.151),burnish*.32)
    rgb=blend(rgb,(.610,.417,.249),dry*.65)
    h=.53+.10*softskin-.068*crease-.027*fd-.16*folds-.031*scuffs+.017*stress
    rough=.52+.045*crease+.040*folds-.115*burnish+.18*dry+.010*fine
    return rgb,h,rough,clamp(1-.105*crease-.21*folds-.12*fd),np.zeros((n,n),np.float32),clamp(burnish*.4+dry+folds*.2)


def plaster(x,y,n):
    layers,lips,polish=trowel_layers(n,701,47)
    aggregate,ad,at=deposits(n,702,85000,(.00014,.0014),aspect=(.5,2),angularity=.40)
    skim=field(n,11,703,2)
    exposed=aggregate*(.15+.85*smooth(skim,-.20,.47))
    pinholes,pd,_=deposits(n,704,15000,(.00014,.00085),location_mask=.35+.65*(1-polish),angularity=.25)
    chips,cd,ct=deposits(n,705,97,(.0018,.020),aspect=(.5,2.5),location_mask=.12+.88*smooth(layers,.3,.65),angularity=.39,floor=.23)
    fissure=fracture_network(n,706,primaries=6,branches=95,width=.00014)*smooth(skim,.04,.5)
    sweep=paths(n,707,320,width=(.00013,.00042),length=(.012,.070),angle=.16,jitter=.12,bend=.011)*polish
    matrix=field(n,560,708,2)
    rgb=color((.762,.742,.693),.025*skim+.018*layers+.005*matrix)
    rgb=blend(rgb,(.678,.631,.536),exposed*smooth(at,.15,.70)*.46+chips*.31)
    rgb=blend(rgb,(.853,.838,.791),exposed*(1-smooth(at,-.7,-.1))*.28+lips*.08)
    rgb=blend(rgb,(.684,.656,.592),pinholes*.055+fissure*.12)
    h=.47+.150*layers+.070*lips+.085*ad*exposed+.008*matrix*(1-.5*polish)-.13*pd-.30*cd-.035*fissure+.005*sweep
    rough=.77-.075*polish+.043*exposed+.065*chips+.025*pinholes+.018*fissure
    return rgb,h,rough,clamp(1-.23*pd-.33*cd-.08*fissure),np.zeros((n,n),np.float32),clamp(chips+exposed*.18+fissure*.35)


def linen(x,y,n):
    threads=120;rng=np.random.default_rng(801)
    gx=x*threads+.16*np.sin(TAU*y*5)+.12*field(n,(12,70),802)
    gy=y*threads+.13*np.sin(TAU*x*7)+.11*field(n,(70,12),803)
    ix=np.floor(gx).astype(int)%threads;iy=np.floor(gy).astype(int)%threads
    ux=gx%1-.5;uy=gy%1-.5
    widths=rng.uniform(.36,.49,(2,threads)).astype(np.float32)
    tones=rng.uniform(-.029,.029,(2,threads)).astype(np.float32)
    tx=TAU*(y*91+ix*.174);ty=TAU*(x*96+iy*.192)
    bx=clamp(1-(ux/(widths[0,ix]*(1+.035*np.sin(tx))))**2)**.72
    by=clamp(1-(uy/(widths[1,iy]*(1+.040*np.sin(ty))))**2)**.72
    a=bx*(.62+.14*np.cos(np.pi*(gy+ix)));b=by*(.62-.14*np.cos(np.pi*(gx+iy)))
    top=np.maximum(a,b);fx=np.zeros((n,n),np.float32);fy=fx.copy()
    for strand in range(16):
        phase=strand*TAU/16
        sx=-.395+strand*.0526+.014*np.sin(tx+phase)
        sy=-.395+strand*.0526+.014*np.sin(ty+phase)
        fx+=np.exp(-((ux-sx)/.011)**2)*np.cos(tx+phase)*.08
        fy+=np.exp(-((uy-sy)/.011)**2)*np.cos(ty+phase)*.08
    fibrils=np.where(a>b,fx,fy)*top
    slubx,_,_=deposits(n,804,650,(.0005,.0014),aspect=(.07,.20),angle=0,spread=.035,angularity=.17)
    sluby,_,_=deposits(n,805,590,(.0005,.0014),aspect=(5,14),angle=0,spread=.035,angularity=.17)
    slubs=slubx*bx+sluby*by
    fibers=paths(n,806,920,width=(.00008,.00018),length=(.001,.010),angle=.4,jitter=2.8,bend=.13)
    worn=smooth(field(n,9,807,2),.23,.54)
    tonal=np.where(a>b,tones[0,ix],tones[1,iy])
    rgb=color((.713,.676,.582),tonal+.012*fibrils+.010*slubs)
    rgb=blend(rgb,(.797,.764,.669),fibers*.22+worn*.07)
    h=.16+.51*top+.052*fibrils+.030*slubs+.025*fibers-.022*worn*top
    rough=.79+.023*(1-top)+.045*worn+.014*slubs
    ao=1-.22*clamp(1-bx*by)
    opacity=clamp(smooth(np.maximum(bx,by),.015,.16)+fibers*.58)
    return rgb,h,rough,ao,np.zeros((n,n),np.float32),clamp(worn*.55+fibers*.3),opacity


from finish_structures import marble,travertine,metal,porcelain,leather,plaster
SPECS[0].update(tile_m=2.8,height_scale_m=.00012,name='Calacatta Oro / Honed Quarry Slab',roughness=.255,description='One procedural 2.8 m slab with narrow shear seams, branching fractures and mineral intergrowth; honed intact finish.')
SPECS[1].update(height_scale_m=.0032,name='Roman Travertine / Honed Vein Cut',description='Fine bed-correlated irregular open pores without painted cavity outlines; shallow honed sedimentary surface.')
SPECS[4].update(name='Champagne Brass / Satin Brushed',description='Intact satin-brushed brass; fine directional tooling and restrained isolated scores.',anisotropy=.35)
SPECS[5].update(name='Blackened Steel / Satin Oxide',description='Intact black oxide on satin machined steel; exposed edges are assigned on fabricated meshes.')
SPECS[6].update(name='Bone Porcelain / Clear Glaze',coat=.32,coat_roughness=.12,description='Clear glaze over warm ivory ceramic; fine firing flow and shallow glaze peel, without random chips.')
SPECS[7].update(name='Saddle Leather / Full Grain',description='Fine irregular submillimeter hide grain and follicles; restrained dye variation. Mesh seams and bends describe construction.')
SPECS[8].update(name='Lime Plaster / Troweled Skim',height_scale_m=.002,description='Overlapping broad trowel passes and polished leading lips with restrained fine aggregate.')
GENERATORS=[marble,travertine,lambda x,y,n:wood(x,y,n),lambda x,y,n:wood(x,y,n,True),metal,lambda x,y,n:metal(x,y,n,True),porcelain,leather,plaster,linen]

def save_rgb(path,a):Image.fromarray(np.rint(clamp(a)*255).astype(np.uint8)).save(path,compress_level=6)
def save_gray(path,a,bit16=False):Image.fromarray(np.rint(clamp(a)*(65535 if bit16 else 255)).astype(np.uint16 if bit16 else np.uint8)).save(path,compress_level=6)
def normal_from_height(h,s,n):
    if s.get('tiling')=='finite-cut':
        dy,dx=np.gradient(h)
        dx*=n*s['height_scale_m']/s['tile_m']
        dy*=n*s['height_scale_m']/s.get('tile_y_m',s['tile_m'])
    else:
        dx=(np.roll(h,-1,1)-np.roll(h,1,1))*(n/2)*s['height_scale_m']/s['tile_m']
        dy=(np.roll(h,-1,0)-np.roll(h,1,0))*(n/2)*s['height_scale_m']/s['tile_m']
    v=np.stack((-dx,dy,np.ones_like(dx)),axis=-1)
    v/=np.linalg.norm(v,axis=-1,keepdims=True)
    return .5+.5*v
def seam(a):
    edge=np.concatenate((abs(a[0]-a[-1]).ravel(),abs(a[:,0]-a[:,-1]).ravel()))
    inside=np.concatenate((abs(a[1]-a[0]).ravel(),abs(a[:,1]-a[:,0]).ravel()))
    return dict(wrap_mean=float(edge.mean()),neighbor_mean=float(inside.mean()),wrap_max=float(edge.max()),neighbor_max=float(inside.max()))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--resolution',type=int,default=4096)
    ap.add_argument('--only',nargs='*',default=[]);ap.add_argument('--draft',action='store_true');args=ap.parse_args()
    requested_resolution=args.resolution
    folder=ROOT/('draft_materials' if args.draft else 'materials');folder.mkdir(exist_ok=True)
    checks=json.loads((folder/'validation.json').read_text()) if args.only and (folder/'validation.json').exists() else {}
    manifest=dict(name='CYBR MATERIAL / Fabrication and Physical Studies',version='3.2-study',resolution=requested_resolution,workflow='metallic-roughness',normal='OpenGL +Y and DirectX -Y',base_color_space='sRGB',data_color_space='linear / Non-Color',height_bit_depth=16,height_midlevel=.5,units='meters',provenance='Ten original deterministic procedural materials. No source images, image generation, white-noise bump or upscaling. Not scans.',materials=SPECS)
    map_types={'BaseColor.png':'sRGB RGB8; intrinsic color, no lighting','Normal_OpenGL.png':'Non-Color tangent +Y RGB8','Normal_DirectX.png':'Non-Color tangent -Y RGB8','Normal_Micro_OpenGL.png':'Non-Color +Y RGB8; residual height only, paired with Height_Macro','Height.png':'Non-Color linear grayscale16; displacement = (height - .5) * height_scale_m','Height_Macro.png':'Non-Color linear grayscale16; band-limited geometry, same physical scale','Roughness.png':'Non-Color linear grayscale8','Metallic.png':'Non-Color grayscale8; exposed metal versus dielectric oxides','AO.png':'Non-Color grayscale8; local cavity approximation, do not bake into color','WearMask.png':'Non-Color grayscale8; material-specific abrasion or weathering','ORM.png':'Non-Color RGB8; AO / roughness / metallic','Opacity.png':'Non-Color grayscale8; woven openings for linen, opaque for other materials'}
    for s,gen in zip(SPECS,GENERATORS):
        if args.only and s['id'] not in args.only:continue
        n=s.get('resolution',requested_resolution)
        x=((np.arange(n,dtype=np.float32)+.5)/n)[None,:];y=((np.arange(n,dtype=np.float32)+.5)/n)[:,None]
        print('GENERATING '+s['id']+' '+str(n)+'px',flush=True)
        out=folder/s['id'];out.mkdir(exist_ok=True)
        values=gen(x,y,n)
        rgb,h,rough,ao,metal,wear=[clamp(v).astype(np.float32) for v in values[:6]]
        opacity=values[6] if len(values)==7 else np.ones((n,n),np.float32)
        mode='nearest' if s.get('tiling')=='finite-cut' else 'wrap'
        h=gaussian_filter(h,.35,mode=mode)
        sigma=max(.5,n/1024)
        if s.get('tiling')=='finite-cut':sigma=max(sigma,.62*s['preview_diameter_m']*1.8/384*n/s['tile_m'])
        macro=gaussian_filter(h,sigma,mode=mode)
        save_rgb(out/'BaseColor.png',rgb);save_gray(out/'Height.png',h,True);save_gray(out/'Height_Macro.png',macro,True)
        save_rgb(out/'Normal_Micro_OpenGL.png',normal_from_height(h-macro,s,n))
        for name,a in [('Roughness',rough),('AO',ao),('Metallic',metal),('WearMask',wear),('Opacity',opacity)]:save_gray(out/(name+'.png'),a)
        save_rgb(out/'ORM.png',np.stack((ao,rough,metal),axis=-1))
        normal=normal_from_height(h,s,n);encoded=np.rint(clamp(normal)*255).astype(np.uint8)
        Image.fromarray(encoded).save(out/'Normal_OpenGL.png',compress_level=6)
        encoded[...,1]=255-encoded[...,1];Image.fromarray(encoded).save(out/'Normal_DirectX.png',compress_level=6)
        checks[s['id']]=dict(height_min=float(h.min()),height_max=float(h.max()),relief_range_m=float((h.max()-h.min())*s['height_scale_m']),roughness_min=float(rough.min()),roughness_max=float(rough.max()),metallic_min=float(metal.min()),metallic_max=float(metal.max()),base_color_seam=seam(rgb),height_seam=seam(h),finite=bool(all(np.isfinite(v).all() for v in values)))
        (out/'material.json').write_text(json.dumps(dict(s,resolution=n,maps=map_types),indent=2)+'\n')
        Image.fromarray(np.rint(rgb*255).astype(np.uint8)).resize((768,768),Image.Resampling.LANCZOS).save(out/'swatch.jpg',quality=95)
        print('  relief %.3f mm; roughness %.2f to %.2f'%((h.max()-h.min())*s['height_scale_m']*1000,rough.min(),rough.max()),flush=True)
        del values,rgb,h,rough,ao,metal,wear,opacity,macro,normal,encoded;gc.collect()
    (folder/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (folder/'validation.json').write_text(json.dumps(checks,indent=2)+'\n')
    print('MAPS_COMPLETE',flush=True)

if __name__=='__main__':main()
