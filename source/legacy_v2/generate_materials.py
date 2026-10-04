#!/usr/bin/env python3
"""AUREL II: seamless, original material structures and localized wear.

Feature randomness varies named shapes at physical scales; no independent
pixel color, white-noise bump, baked light, third-party images or upscaling.
"""
from __future__ import annotations
import argparse, gc, json
from pathlib import Path
import numpy as np
from PIL import Image, ImageDraw
from scipy.ndimage import gaussian_filter, map_coordinates, zoom
from scipy.interpolate import CubicSpline
ROOT=Path(__file__).resolve().parents[1]
TAU=np.float32(2*np.pi)
SPECS=[
 dict(id='01_calacatta_oro',name='Calacatta Oro / Weathered',family='Stone',tile_m=.8,height_scale_m=.0012,roughness=.31,ior=1.52,metallic=0,coat=.025,coat_roughness=.3,preview_diameter_m=.24,description='Brecciated ivory marble: branching graphite and ochre seams, crystalline inclusions, mineral fissures and abraded pits.'),
 dict(id='02_roman_travertine',name='Roman Travertine / Open Pore',family='Stone',tile_m=.65,height_scale_m=.005,roughness=.59,ior=1.5,metallic=0,preview_diameter_m=.24,description='Weathered vein-cut travertine with varied sediment beds, torn elongated cavities, fossil-like mineral fragments and unfilled pores.'),
 dict(id='03_american_walnut',name='American Walnut / Timeworn',family='Wood',tile_m=.55,height_scale_m=.001,roughness=.44,ior=1.5,metallic=0,coat=.045,coat_roughness=.37,anisotropy=.16,preview_diameter_m=.18,description='Rich walnut with knot-distorted growth, broken longitudinal fibers, open vessels, dry checks and rubbed oil finish.'),
 dict(id='04_fumed_oak',name='Fumed Oak / Open Grain',family='Wood',tile_m=.45,height_scale_m=.0014,roughness=.51,ior=1.5,metallic=0,coat=.02,coat_roughness=.45,anisotropy=.13,preview_diameter_m=.16,description='Dark oak with irregular earlywood vessels, medullary ribbons, splits, flattened grain peaks and pale worn fibers.'),
 dict(id='05_champagne_brass',name='Champagne Brass / Aged',family='Metal',tile_m=.2,height_scale_m=.00022,roughness=.32,ior=1.5,metallic=1,anisotropy=.43,preview_diameter_m=.10,description='Brushed brass with discrete crossing scratches, shallow dents, warm tarnish islands and small verdigris deposits. Exposed metal and dielectric oxide have separate metalness.'),
 dict(id='06_blackened_steel',name='Blackened Steel / Worn',family='Metal',tile_m=.25,height_scale_m=.00045,roughness=.44,ior=1.5,metallic=0,anisotropy=.22,preview_diameter_m=.12,description='Black oxide steel with exposed silver wear, directionally scored tooling, impact pits and localized rust. Bare steel is metallic; oxide and rust are dielectric.'),
 dict(id='07_bone_porcelain',name='Bone Porcelain / Crazed',family='Ceramic',tile_m=.24,height_scale_m=.00055,roughness=.21,ior=1.49,metallic=0,coat=.24,coat_roughness=.17,preview_diameter_m=.10,description='Ivory glaze with fine irregular crazing, pinholes, isolated chips exposing ceramic body, and rubbed glaze patches.'),
 dict(id='08_saddle_leather',name='Saddle Leather / Used',family='Leather',tile_m=.18,height_scale_m=.0011,roughness=.46,ior=1.48,metallic=0,coat=.025,coat_roughness=.38,preview_diameter_m=.065,description='Cognac leather with fine stretched grain, follicle pits, larger compressed folds, pale abrasion and short surface scuffs.'),
 dict(id='09_lime_plaster',name='Lime Plaster / Weathered',family='Plaster',tile_m=.65,height_scale_m=.003,roughness=.73,ior=1.48,metallic=0,preview_diameter_m=.22,description='Layered lime plaster with trowel lips, mineral aggregate, pinholes, chipped laminae and fine shrinkage checks.'),
 dict(id='10_natural_linen',name='Natural Linen / Frayed',family='Textile',tile_m=.085,height_scale_m=.0008,roughness=.78,ior=1.46,metallic=0,sheen=.34,sheen_roughness=.65,preview_diameter_m=.028,description='Irregular flax weave with individual twisted fibrils, variable yarn widths, attached slubs, loose fibers and small frayed wear patches.'),
]

def clamp(a,lo=0,hi=1):return np.clip(a,lo,hi)
def smooth(a,lo,hi):
 t=clamp((a-lo)/(hi-lo));return t*t*(3-2*t)

def field(n,cells,seed,octaves=1):
 """Periodic process variation, cubic interpolated at finite feature scales."""
 rng=np.random.default_rng(seed);out=np.zeros((n,n),np.float32);norm=0
 shape=(cells,cells) if isinstance(cells,int) else cells
 for i in range(octaves):
  sy,sx=[int(v*2**i) for v in shape]
  grid=rng.uniform(-1,1,(sy,sx)).astype(np.float32)
  a=zoom(grid,(n/sy,n/sx),order=3,mode='grid-wrap',grid_mode=True)
  weight=0.5**i;out+=a*weight;norm+=weight
 return out/np.float32(norm)

def color(base,mod):return np.stack([v+mod for v in base],axis=-1).astype(np.float32)
def blend(rgb,tint,mask):
 for c,v in enumerate(tint):rgb[...,c]=rgb[...,c]*(1-mask)+v*mask
 return rgb

def specks(n,seed,count,rx_range,ry_range,power=.75):
 """Ragged compact inclusions/cavities with varied, finite-support outlines."""
 rng=np.random.default_rng(seed);out=np.zeros((n,n),np.float32)
 for _ in range(count):
  cx,cy=rng.uniform(0,n,2);rx,ry=rng.uniform(*rx_range)*n,rng.uniform(*ry_range)*n
  ix=np.arange(int(cx-rx*1.3-1),int(cx+rx*1.3+2));iy=np.arange(int(cy-ry*1.3-1),int(cy+ry*1.3+2))
  u=(ix[None,:]-cx)/max(rx,.6);v=(iy[:,None]-cy)/max(ry,.6)
  angle=np.arctan2(v,u);phase=rng.uniform(0,TAU)
  r=(u*u+v*v)*(1+.17*np.sin(3*angle+phase)+.09*np.sin(7*angle-phase))
  patch=clamp(1-r).astype(np.float32)**power*rng.uniform(.55,1)
  target=np.ix_(iy%n,ix%n);out[target]=np.maximum(out[target],patch)
 return out

def paths(n,seed,count,width=(.0003,.001),length=(.02,.15),angle=np.pi/2,jitter=.2,bend=.08):
 """Finite scratches, checks or fibers; curved paths wrapped around the tile."""
 rng=np.random.default_rng(seed);im=Image.new('L',(n,n));draw=ImageDraw.Draw(im)
 for _ in range(count):
  center=rng.uniform(0,n,2);theta=angle+rng.uniform(-jitter,jitter);size=rng.uniform(*length)*n
  t=np.linspace(-.5,.5,7);control=rng.uniform(-1,1,7)*size*bend;control[0]*=.3;control[-1]*=.3
  fine=np.linspace(-.5,.5,64);off=CubicSpline(t,control)(fine)
  xx=center[0]+np.cos(theta)*fine*size-np.sin(theta)*off
  yy=center[1]+np.sin(theta)*fine*size+np.cos(theta)*off
  line_width=max(1,round(rng.uniform(*width)*n));tone=int(rng.uniform(120,255))
  for oy in (-n,0,n):
   for ox in (-n,0,n):draw.line(list(zip(xx+ox,yy+oy)),fill=tone,width=line_width)
 return gaussian_filter(np.asarray(im,dtype=np.float32)/255,.42,mode='wrap')

def warp_map(a,u,v):
 n=a.shape[0];u,v=np.broadcast_arrays(u,v)
 return map_coordinates(a,np.stack([(v%1)*n-.5,(u%1)*n-.5]),order=1,mode='grid-wrap',prefilter=False)

def voronoi(n,cells,seed,u=None,v=None):
 """Distance gap, seed distance and cell variation on an irregular tiling."""
 rng=np.random.default_rng(seed);jx=rng.uniform(-.46,.46,(cells,cells)).astype(np.float32);jy=rng.uniform(-.46,.46,(cells,cells)).astype(np.float32)
 values=rng.uniform(-1,1,(cells,cells)).astype(np.float32)
 if u is None:u=(np.arange(n,dtype=np.float32)[None,:]+.5)/n
 if v is None:v=(np.arange(n,dtype=np.float32)[:,None]+.5)/n
 xx=np.broadcast_to(u*cells,(n,n));yy=np.broadcast_to(v*cells,(n,n));ix=np.floor(xx).astype(np.int16);iy=np.floor(yy).astype(np.int16)
 first=np.full((n,n),100,np.float32);second=first.copy();best=np.zeros((n,n),np.float32)
 for dy in (-1,0,1):
  for dx in (-1,0,1):
   cx=(ix+dx)%cells;cy=(iy+dy)%cells
   ds=(xx-(ix+dx+.5+jx[cy,cx]))**2+(yy-(iy+dy+.5+jy[cy,cx]))**2
   low=ds<first;second=np.where(low,first,np.minimum(second,ds));best=np.where(low,values[cy,cx],best);first=np.minimum(first,ds)
 return np.sqrt(second)-np.sqrt(first),np.sqrt(first),best

def marble(x,y,n):
 drift=field(n,4,113,3);fracture=field(n,24,23,2);coord=2*x+y+.35*drift+.032*fracture
 distance=np.abs(np.sin(np.pi*coord));width=.055+.12*smooth(field(n,9,25),-.4,.5)
 vein=np.exp(-(distance/width)**2);shoulder=np.exp(-(distance/(width*2.8))**2)*.24
 branchfield=field(n,9,54,3)+.18*drift
 branch=np.exp(-((branchfield-.09)/.022)**2)
 hairline=np.exp(-((branchfield+.18)/.009)**2)*smooth(field(n,13,953),-.15,.4)
 gap,grain,stone=voronoi(n,257,811);mineral=field(n,48,833,2);broken=clamp(.56+.44*mineral)
 rgb=color((.82,.796,.723),.045*field(n,12,931,2)+.023*stone)
 rgb=blend(rgb,(.33,.349,.346),clamp(vein*(.51+.34*broken)+shoulder+branch*.38))
 crystalline=smooth(mineral,.08,.48)*vein;rgb=blend(rgb,(.64,.637,.574),crystalline*.65)
 crystals=specks(n,933,3400,(.0007,.0038),(.0006,.0045),power=.7)*(.25+.75*smooth(vein+branch,.05,.6))
 rgb=blend(rgb,(.63,.615,.541),crystals*.43)
 rgb=blend(rgb,(.417,.403,.352),hairline*.47)
 ochre=np.exp(-((distance-width*.91)/(.02+width*.15))**2)*smooth(field(n,5,335),-.05,.42)
 ochre=clamp(ochre+branch*.36*smooth(field(n,8,955),.1,.45));rgb=blend(rgb,(.49,.316,.126),ochre*.71)
 fissure=np.exp(-(gap/.027)**2)*smooth(vein+branch,.28,.7)
 pits=specks(n,910,4600,(.00022,.0013),(.0002,.001),power=.65)
 abrasion=paths(n,223,85,width=(.00022,.0007),length=(.004,.055),angle=.22,jitter=1.2,bend=.018)
 wear=clamp(pits*.75+abrasion*.42+fissure*.6);rgb=blend(rgb,(.365,.35,.291),fissure*.48+pits*.23);rgb=blend(rgb,(.85,.835,.778),abrasion*.28)
 h=.52+.023*field(n,80,129,2)+.012*stone-.26*pits-.08*fissure-.035*abrasion-.026*hairline-.018*crystals
 rough=.26+.12*vein+.055*branch+.04*crystals+.13*wear+.04*field(n,12,428)
 return rgb,h,rough,1-.3*pits-.11*fissure,np.zeros((n,n),np.float32),wear

def travertine(x,y,n):
 flow=.017*field(n,(7,18),164,2);coord=38*(y+flow+.008*np.sin(TAU*x*2));rng=np.random.default_rng(234);tones=rng.uniform(-1,1,38).astype(np.float32)
 idx=np.floor(coord).astype(np.int16)%38;frac=coord%1;bands=tones[idx]*(1-smooth(frac,.65,1))+tones[(idx+1)%38]*smooth(frac,.65,1)
 fine=.022*field(n,(55,8),245,2);pits=specks(n,164,4600,(.00065,.018),(.0003,.0048),power=.5);micro=specks(n,1194,13500,(.00015,.0013),(.0002,.001),power=.65)
 pits*=.35+.65*smooth(bands,-.35,.65);gaps=np.exp(-((frac-.06)/.038)**2)
 checks=paths(n,1641,160,width=(.00022,.00085),length=(.008,.09),angle=0,jitter=.13,bend=.055)
 rgb=color((.708,.625,.475),.085*bands+fine+.013*field(n,190,612));rgb=blend(rgb,(.334,.246,.127),clamp(pits*.79+micro*.41+checks*.27));rgb=blend(rgb,(.803,.739,.612),gaps*.22)
 wear=clamp(pits+checks*.65);h=.66+.026*bands+.012*field(n,230,745)-.58*pits-.15*micro-.055*checks
 rough=.53+.18*pits+.06*micro+.05*bands+.05*checks
 return rgb,h,rough,clamp(1-.64*pits-.22*micro-.12*checks),np.zeros((n,n),np.float32),wear

def wood(x,y,n,oak=False):
 seed=89 if oak else 21;count=68 if oak else 43;rng=np.random.default_rng(seed)
 drift=.011*field(n,(4,7),seed+1,2)+.007*np.sin(TAU*y);u=np.broadcast_to(x+drift,(n,n)).copy();v=np.broadcast_to(y,(n,n))
 knots=np.zeros((n,n),np.float32);knot_pattern=knots.copy();heart=knots.copy()
 for cx,cy,rx,ry in ((.29,.33,.065,.11),(.76,.80,.046,.065),(.62,.12,.026,.042)) if not oak else ((.71,.39,.023,.04),(.19,.87,.032,.053)):
  dx=(x-cx+.5)%1-.5;dy=(y-cy+.5)%1-.5;radius=np.sqrt((dx/rx)**2+(dy/ry)**2)
  u+=rx*.77*np.tanh(dx/(rx*.38))*np.exp(-((dx/(rx*2))**2+(dy/(ry*2.8))**2))
  km=1-smooth(radius,.78,1.9);knots=np.maximum(knots,km)
  theta=np.arctan2(dy/ry,dx/rx)
  distorted=radius*(1+.10*np.sin(theta*3+.8)+.07*np.sin(theta*7+cy*11))+.11*field(n,29,seed+85,2)
  pattern=(.5+.5*np.sin(TAU*(distorted*(5.8+cx*4)+.14*field(n,41,seed+86))))*(.75+.25*field(n,34,seed+87))
  knot_pattern=np.maximum(knot_pattern,pattern*km);heart=np.maximum(heart,np.exp(-(distorted/.34)**2))
 breaks=field(n,(42,270),seed+13,2);coord=u*count+.13*field(n,(6,36),seed+14)
 idx=np.floor(coord).astype(np.int16)%count;frac=coord%1;tone=rng.uniform(-.8,1,count).astype(np.float32);width=rng.uniform(.065,.20,count).astype(np.float32)
 band=tone[idx]*.65+.35*np.cos(frac*np.pi);late=np.exp(-((frac-.14)/width[idx])**2)*(.72+.28*breaks)
 fiber_grid=field(n,(28,512),seed+37,2);fiber=warp_map(fiber_grid,u,v)
 vessel=specks(n,seed+9,25000 if oak else 20500,(.00016,.00054),(.00065,.014),power=.7);vessel=warp_map(vessel,u,v)*(.48+.52*smooth(late,.05,.7))
 checks=warp_map(paths(n,seed+39,155,width=(.00025,.0011),length=(.006,.14),angle=np.pi/2,jitter=.024,bend=.012),u,v)
 worn=smooth(field(n,(6,10),seed+119,2),.13,.55)*(.42+.58*smooth(fiber,-.3,.5))
 scratches=paths(n,seed+70,90,width=(.0002,.0008),length=(.008,.068),angle=1.2,jitter=.7,bend=.018)
 base=(.338,.258,.165) if oak else (.459,.296,.166);dark=(.125,.079,.032) if oak else (.198,.091,.031)
 rgb=color(base,.069*band+.031*fiber+.017*breaks);rgb=blend(rgb,dark,clamp(late*.51+vessel*.52+checks*.64))
 rgb=blend(rgb,(.255,.137,.063) if not oak else (.252,.181,.093),knots*(.15+.45*knot_pattern));rgb=blend(rgb,(.078,.041,.017),heart*.85)
 if oak:rays=warp_map(specks(n,seed+21,950,(.004,.034),(.0003,.0012),power=.8),u,v);rgb=blend(rgb,(.566,.425,.234),rays*.45)
 else:rays=np.zeros((n,n),np.float32)
 exposed=(.551,.435,.266) if oak else (.612,.438,.255);rgb=blend(rgb,exposed,clamp(worn*.28+scratches*.52))
 wear=clamp(worn*.6+checks*.75+scratches);h=.56+.082*band+.032*fiber-.30*vessel-.35*checks-.09*scratches+.018*rays;h+=knots*.04*knot_pattern-.18*heart
 rough=(.49 if oak else .39)+.13*vessel+.11*checks+.17*worn+.12*scratches-.015*late
 return rgb,h,rough,clamp(1-.30*vessel-.34*checks),np.zeros((n,n),np.float32),wear

def metal(x,y,n,steel=False):
 seed=619 if steel else 610;rng=np.random.default_rng(seed);profile=rng.uniform(-1,1,650).astype(np.float32)
 p=zoom(profile,n/650,order=3,mode='grid-wrap',grid_mode=True)[None,:];lines=np.broadcast_to(p,(n,n))*(.6+.4*field(n,(35,75),seed+1))
 brush=paths(n,seed+2,640,width=(.00012,.00045),length=(.06,.46),angle=np.pi/2,jitter=.005,bend=.001)
 scratches=paths(n,seed+4,85,width=(.00025,.0014),length=(.016,.28),angle=.23,jitter=1.3,bend=.012);scuffs=paths(n,seed+5,95,width=(.0006,.0035),length=(.008,.065),angle=1.2,jitter=.85,bend=.006)
 dents=specks(n,seed+8,95,(.0015,.008),(.001,.006),power=.7);domains=field(n,9,seed+12,3)+.24*field(n,50,seed+13)
 oxide=smooth(domains,.08,.39);reaction=field(n,130,seed+14,2);oxide*=.63+.37*smooth(reaction,-.5,.5)
 if steel:
  exposed=clamp(scratches*.85+scuffs*.92+smooth(field(n,12,seed+20),.49,.7)*.75);rust=smooth(domains,.3,.54)*(1-exposed)
  rgb=color((.115,.134,.142),.015*field(n,18,seed+15)+.006*lines);rgb=blend(rgb,(.596,.621,.625),exposed);rgb=blend(rgb,(.395,.177,.067),rust*(.7+.3*reaction))
  metallic=exposed*(1-rust);rough=.40+.085*lines+.19*rust-.21*exposed+.09*oxide
  h=.52+.022*lines-.06*brush-.14*scratches-.24*dents+.10*rust*reaction;wear=clamp(exposed+rust+dents*.5)
 else:
  rgb=color((.771,.633,.395),.014*lines+.014*field(n,15,seed+15));rgb=blend(rgb,(.361,.228,.092),oxide*.73)
  patina=smooth(domains,.39,.63)*smooth(reaction,-.05,.5);rgb=blend(rgb,(.187,.326,.252),patina*.92);rgb=blend(rgb,(.855,.755,.516),scratches*.45+scuffs*.22)
  metallic=clamp(1-oxide*.92-patina*.8);rough=.27+.07*lines+.30*oxide+.18*patina-.11*scratches+.12*scuffs
  h=.53+.03*lines-.08*brush-.22*scratches-.22*dents+.06*oxide*reaction;wear=clamp(oxide+scratches*.7+scuffs*.5+dents)
 return rgb,h,rough,1-.28*dents,metallic,wear

def porcelain(x,y,n):
 u=x+.009*field(n,8,788,2);v=y+.009*field(n,8,789,2);gap,dist,cell=voronoi(n,37,87,u,v)
 craze=np.exp(-(gap/(.016+.008*field(n,14,875)))**2);subgap,_,_=voronoi(n,121,881,u,v);microcrack=np.exp(-(subgap/.016)**2)*smooth(field(n,12,887),.12,.45)
 pinholes=specks(n,872,2900,(.00026,.0013),(.00025,.001),power=.6);chips=specks(n,873,85,(.0018,.012),(.001,.007),power=.55)*smooth(field(n,11,810),-.1,.4)
 rub=smooth(field(n,7,813,2),.17,.53);scratches=paths(n,814,70,width=(.00012,.00045),length=(.004,.036),angle=.8,jitter=1.7,bend=.013)
 rgb=color((.811,.783,.698),.013*field(n,14,878,2)+.007*cell);rgb=blend(rgb,(.492,.425,.311),clamp(craze*.41+microcrack*.18+pinholes*.23));rgb=blend(rgb,(.584,.483,.337),chips*.82);rgb=blend(rgb,(.869,.841,.766),rub*.12+scratches*.17)
 h=.54+.022*field(n,31,879,2)-.058*craze-.026*microcrack-.25*pinholes-.62*chips-.022*scratches;wear=clamp(chips+rub*.65+scratches*.65)
 rough=.16+.12*craze+.05*microcrack+.29*chips+.22*rub+.10*scratches+.025*pinholes
 return rgb,h,rough,clamp(1-.10*craze-.24*pinholes-.38*chips),np.zeros((n,n),np.float32),wear

def leather(x,y,n):
 foldfield=field(n,(7,5),411,3)
 fold_paths=paths(n,412,28,width=(.0011,.0042),length=(.10,.60),angle=.7,jitter=.6,bend=.10)
 folds=gaussian_filter(fold_paths,max(.5,n*.0006),mode='wrap')
 folds+=.38*np.exp(-(np.abs(np.sin(np.pi*(13*x+5*y+.44*foldfield)))/.065)**2)*smooth(field(n,7,413),.05,.4)
 folds=clamp(folds)
 u=x+.007*field(n,17,415,2)+.0014*field(n,95,423)+.004*foldfield
 v=y+.005*field(n,17,416,2)+.0012*field(n,87,424)
 gap,dist,cell=voronoi(n,153,417,u,v)
 grain=np.exp(-(gap/(.09+.019*cell))**2)
 dome=gaussian_filter(1-np.exp(-(gap/.32)**1.12),max(.5,n/153*.12),mode='wrap')
 smallgap,_,smallcell=voronoi(n,393,418,u,v)
 fine=np.exp(-(smallgap/.036)**2)*smooth(field(n,27,425),-.05,.4)
 pores=specks(n,419,22500,(.00019,.00065),(.00019,.0008),power=.9);scuffs=paths(n,420,85,width=(.00022,.0012),length=(.003,.04),angle=.35,jitter=1.7,bend=.04)
 worn=smooth(field(n,(7,9),421,2),.07,.50)*(1-grain*.45);mellowness=field(n,11,422,2)
 rgb=color((.449,.226,.106),.021*mellowness+.007*cell);rgb=blend(rgb,(.257,.089,.027),grain*.26+folds*.30+pores*.13);rgb=blend(rgb,(.626,.385,.193),worn*.44+scuffs*.43)
 wear=clamp(worn*.8+scuffs+folds*.13)
 h=.40+.13*dome*(.9+.20*cell)-.013*grain-.007*fine-.028*pores-.22*folds-.055*scuffs+.025*foldfield
 rough=.47+.07*grain+.08*folds+.14*worn+.12*scuffs+.013*smallcell
 return rgb,h,rough,clamp(1-.16*grain-.19*folds-.08*pores),np.zeros((n,n),np.float32),wear

def plaster(x,y,n):
 trowel=field(n,(7,9),12,2);lips=np.exp(-((trowel+.13*field(n,17,18))/.034)**2);deposit=smooth(trowel,-.35,.4);aggregate=field(n,340,19)
 pinholes=specks(n,20,18000,(.00022,.0018),(.0002,.0016),power=.6);chips=specks(n,21,140,(.002,.025),(.0015,.016),power=.45)*smooth(field(n,14,22),-.15,.35)
 fissures=paths(n,23,195,width=(.0002,.00075),length=(.007,.14),angle=.7,jitter=2,bend=.11);sweep=paths(n,24,120,width=(.00065,.002),length=(.035,.26),angle=.24,jitter=.8,bend=.3)
 rgb=color((.729,.698,.619),.075*trowel+.019*aggregate+.012*field(n,70,26));rgb=blend(rgb,(.536,.493,.407),chips*.69+pinholes*.26+fissures*.20);rgb=blend(rgb,(.805,.787,.727),lips*.25+sweep*.10)
 h=.46+.18*deposit+.027*aggregate+.10*lips-.22*pinholes-.34*chips-.10*fissures+.025*sweep;wear=clamp(chips+fissures*.8);rough=.75-.065*deposit+.065*chips+.035*pinholes-.035*sweep
 return rgb,h,rough,clamp(1-.3*pinholes-.29*chips-.19*fissures),np.zeros((n,n),np.float32),wear

def linen(x,y,n):
 threads=82;rng=np.random.default_rng(983);gx=x*threads+.095*np.sin(TAU*y*5)+.065*field(n,(8,16),984);gy=y*threads+.09*np.sin(TAU*x*7)+.06*field(n,(16,8),985)
 ix=np.floor(gx).astype(np.int16)%threads;iy=np.floor(gy).astype(np.int16)%threads;ux=gx%1-.5;uy=gy%1-.5
 widths=rng.uniform(.345,.49,(2,threads)).astype(np.float32);tones=rng.uniform(-.065,.065,(2,threads)).astype(np.float32);twistx=TAU*(y*37+ix*.172);twisty=TAU*(x*43+iy*.183)
 bx=clamp(1-(ux/(widths[0,ix]*(1+.06*np.sin(twistx))))**2)**.58;by=clamp(1-(uy/(widths[1,iy]*(1+.07*np.sin(twisty))))**2)**.58
 over=((ix+iy)%2)==0;a=bx*(.73+.27*over);b=by*(1-.27*over);top=np.maximum(a,b);fibersx=np.zeros((n,n),np.float32);fibersy=fibersx.copy()
 for strand in range(10):
  phase=strand*TAU/10;fx=-.38+strand*.084+.017*np.sin(twistx+phase);fy=-.38+strand*.084+.019*np.sin(twisty+phase)
  fibersx+=np.exp(-((ux-fx)/.019)**2)*np.cos(twistx+phase)*.1;fibersy+=np.exp(-((uy-fy)/.018)**2)*np.cos(twisty+phase)*.1
 fibers=np.where(a>b,fibersx,fibersy)*top
 slubs=specks(n,986,4200,(.0005,.0018),(.0015,.007),power=.75)*bx;slubs+=specks(n,987,3800,(.0015,.006),(.0004,.0016),power=.75)*by
 loose=paths(n,988,1350,width=(.00008,.00023),length=(.0015,.017),angle=.6,jitter=2.5,bend=.18);worn=smooth(field(n,7,989,2),.32,.63)
 tonal=np.where(a>b,tones[0,ix],tones[1,iy]);rgb=color((.672,.626,.510),tonal+.023*fibers+.017*slubs);rgb=blend(rgb,(.786,.742,.624),loose*.36+worn*.13)
 h=.16+.54*top+.12*fibers+.055*slubs+.04*loose-.09*worn*top;rough=.755+.04*(1-top)+.055*worn+.013*fibers;ao=1-.32*clamp(1-bx*by)
 return rgb,h,rough,ao,np.zeros((n,n),np.float32),clamp(worn*.8+loose*.6)

GENERATORS=[marble,travertine,lambda x,y,n:wood(x,y,n),lambda x,y,n:wood(x,y,n,True),metal,lambda x,y,n:metal(x,y,n,True),porcelain,leather,plaster,linen]
def save_rgb(path,data):Image.fromarray(np.rint(clamp(data)*255).astype(np.uint8)).save(path,compress_level=6)
def save_gray(path,data,bit16=False):Image.fromarray(np.rint(clamp(data)*(65535 if bit16 else 255)).astype(np.uint16 if bit16 else np.uint8)).save(path,compress_level=6)
def normal_from_height(h,spec,n):
 dx=(np.roll(h,-1,1)-np.roll(h,1,1))*(n/2)*spec['height_scale_m']/spec['tile_m'];dy=(np.roll(h,-1,0)-np.roll(h,1,0))*(n/2)*spec['height_scale_m']/spec['tile_m']
 normal=np.stack((-dx,dy,np.ones_like(dx)),axis=-1);normal/=np.sqrt(np.sum(normal*normal,axis=-1,keepdims=True));return .5+.5*normal
def analyze_seam(a):
 edge=np.concatenate((np.abs(a[0]-a[-1]).ravel(),np.abs(a[:,0]-a[:,-1]).ravel()));inside=np.concatenate((np.abs(a[1]-a[0]).ravel(),np.abs(a[:,1]-a[:,0]).ravel()))
 return dict(wrap_mean=float(edge.mean()),neighbor_mean=float(inside.mean()),wrap_max=float(edge.max()),neighbor_max=float(inside.max()))
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--resolution',type=int,default=4096);ap.add_argument('--only',nargs='*',default=[]);ap.add_argument('--draft',action='store_true');args=ap.parse_args()
 n=args.resolution;x=((np.arange(n,dtype=np.float32)+.5)/n)[None,:];y=((np.arange(n,dtype=np.float32)+.5)/n)[:,None]
 base=ROOT/('draft_materials' if args.draft else 'materials');base.mkdir(exist_ok=True);checks=json.loads((base/'validation.json').read_text()) if args.only and (base/'validation.json').exists() else {}
 manifest=dict(name='AUREL II / Natural Detail + Wear',version='2.0',resolution=n,workflow='metallic-roughness',normal='OpenGL +Y; DirectX -Y also supplied',base_color_space='sRGB',data_color_space='linear / Non-Color',height_bit_depth=16,height_midlevel=.5,units='meters',provenance='Original deterministic procedural materials, not scans. No external textures or image upscaling.',materials=SPECS)
 for spec,gen in zip(SPECS,GENERATORS):
  if args.only and spec['id'] not in args.only:continue
  print(f"Generating {spec['name']} at {n} x {n}",flush=True);folder=base/spec['id'];folder.mkdir(exist_ok=True)
  rgb,h,rough,ao,metallic,wear=gen(x,y,n);rgb,h,rough,ao,metallic,wear=[clamp(a) for a in (rgb,h,rough,ao,metallic,wear)];h=gaussian_filter(h,.35,mode='wrap')
  save_rgb(folder/'BaseColor.png',rgb);save_gray(folder/'Height.png',h,True)
  macro=gaussian_filter(h,max(.5,n/1024),mode='wrap')
  save_gray(folder/'Height_Macro.png',macro,True)
  save_rgb(folder/'Normal_Micro_OpenGL.png',normal_from_height(h-macro,spec,n))
  for name,data in [('Roughness',rough),('AO',ao),('Metallic',metallic),('WearMask',wear)]:save_gray(folder/(name+'.png'),data)
  save_rgb(folder/'ORM.png',np.stack((ao,rough,metallic),axis=-1));normal=normal_from_height(h,spec,n);encoded=np.rint(clamp(normal)*255).astype(np.uint8)
  Image.fromarray(encoded).save(folder/'Normal_OpenGL.png',compress_level=6);encoded[...,1]=255-encoded[...,1];Image.fromarray(encoded).save(folder/'Normal_DirectX.png',compress_level=6)
  checks[spec['id']]=dict(height_min=float(h.min()),height_max=float(h.max()),relief_range_m=float((h.max()-h.min())*spec['height_scale_m']),roughness_min=float(rough.min()),roughness_max=float(rough.max()),metallic_min=float(metallic.min()),metallic_max=float(metallic.max()),base_color_seam=analyze_seam(rgb),height_seam=analyze_seam(h),finite=bool(all(np.isfinite(a).all() for a in (rgb,h,rough,ao,normal,metallic,wear))))
  meta=dict(spec,resolution=n,maps={'BaseColor.png':'sRGB, RGB8; intrinsic color, no illumination','Normal_OpenGL.png':'Non-Color, tangent +Y RGB8','Normal_DirectX.png':'Non-Color, tangent -Y RGB8','Roughness.png':'Non-Color, linear grayscale8','Metallic.png':'Non-Color, linear grayscale8; metal/oxide coverage','AO.png':'Non-Color, grayscale8; local cavity only','Height.png':'Non-Color, linear grayscale16; displacement = (height - .5) * height_scale_m','ORM.png':'Non-Color, RGB8: AO / roughness / metalness','WearMask.png':'Non-Color, grayscale8; abrasion, exposure, chips or patina'})
  meta['maps']['Height_Macro.png']='Non-Color, grayscale16; band-limited geometry height, same physical scale'
  meta['maps']['Normal_Micro_OpenGL.png']='Non-Color, +Y RGB8; residual relief for use with Height_Macro only'
  (folder/'material.json').write_text(json.dumps(meta,indent=2)+'\n');Image.fromarray(np.rint(rgb*255).astype(np.uint8)).resize((768,768),Image.Resampling.LANCZOS).save(folder/'swatch.jpg',quality=96)
  print(f"  relief {checks[spec['id']]['relief_range_m']*1e6:.1f} um; roughness {rough.min():.2f}–{rough.max():.2f}",flush=True)
  del rgb,h,rough,ao,normal,encoded,metallic,wear;gc.collect()
 (base/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');(base/'validation.json').write_text(json.dumps(checks,indent=2)+'\n');print('Map generation complete.',flush=True)
if __name__=='__main__':main()
