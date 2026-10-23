"""Live deterministic recipes extracted from pinned MATERIAL 3.2 source.

Six earlier shadowed definitions are intentionally absent. The immutable
vendor/original snapshot remains the historical comparison oracle.
"""
from pathlib import Path
import json
import numpy as np
from .surface_math import clamp,smooth,field,color,blend,paths,warp_map,TAU
from .physical_features import deposits


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

_wood_profiles=json.loads(Path(__file__).with_name('wood_profiles.json').read_text())

for _spec,_profile in zip(SPECS[2:4],(_wood_profiles['american_walnut'],_wood_profiles['fumed_oak'])):
    _spec.update(tile_m=_profile['extent_m'],tile_y_m=_profile['extent_m'],height_scale_m=_profile['height_scale_m'],roughness=_profile['roughness'])

def wood(x,y,n,oak=False):
    from .wood_anatomy import wood as procedural_wood
    return procedural_wood(x,y,n,oak)

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

from .finish_structures import marble,travertine,metal,porcelain,leather,plaster

SPECS[0].update(tile_m=2.8,height_scale_m=.00012,name='Calacatta Oro / Honed Quarry Slab',roughness=.255,description='One procedural 2.8 m slab with narrow shear seams, branching fractures and mineral intergrowth; honed intact finish.')

SPECS[1].update(height_scale_m=.0032,name='Roman Travertine / Honed Vein Cut',description='Fine bed-correlated irregular open pores without painted cavity outlines; shallow honed sedimentary surface.')

SPECS[4].update(name='Champagne Brass / Satin Brushed',description='Intact satin-brushed brass; fine directional tooling and restrained isolated scores.',anisotropy=.35)

SPECS[5].update(name='Blackened Steel / Satin Oxide',description='Intact black oxide on satin machined steel; exposed edges are assigned on fabricated meshes.')

SPECS[6].update(name='Bone Porcelain / Clear Glaze',coat=.32,coat_roughness=.12,description='Clear glaze over warm ivory ceramic; fine firing flow and shallow glaze peel, without random chips.')

SPECS[7].update(name='Saddle Leather / Full Grain',description='Fine irregular submillimeter hide grain and follicles; restrained dye variation. Mesh seams and bends describe construction.')

SPECS[8].update(name='Lime Plaster / Troweled Skim',height_scale_m=.002,description='Overlapping broad trowel passes and polished leading lips with restrained fine aggregate.')

GENERATORS=[marble,travertine,lambda x,y,n:wood(x,y,n),lambda x,y,n:wood(x,y,n,True),metal,lambda x,y,n:metal(x,y,n,True),porcelain,leather,plaster,linen]
