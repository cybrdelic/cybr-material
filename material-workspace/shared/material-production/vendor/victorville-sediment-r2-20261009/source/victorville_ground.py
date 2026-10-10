"""CYBR 36: procedural Victorville-area alluvial desert ground, SI units.
No raster inputs. Authored regional interpretation, not a measured site.
The packed grains, mineral identity, dust exposure and metric relief share state.
"""
from __future__ import annotations
import math, hashlib
import numpy as np
from scipy.ndimage import zoom, gaussian_filter

DEFAULT = dict(id='36_victorville_desert_ground', name='Victorville / Granitic desert ground',
 family='Desert ground', variant='Dust and sand mantle over compacted dirt and buried alluvial gravel',
 tile_m=.60, height_scale_m=.032, seed=92394, roughness=.87, metallic=0., ior=1.50,
 gravel_density=1.0, crust_strength=.32, dust_amount=.74, sediment_depth_m=.0016, sand_deposit_m=.0012, clod_density=2200., compaction_smoothing_m=.00085, loose_sand_scale=.15, 
 color=[.590,.535,.455], representation='metric displaced ground and coupled PBR fields',
 physical_values='Authored approximations; not laboratory-measured or site-calibrated',
 tiling='periodic toroidal field', region='Victorville, western Mojave, California')

# Intrinsic sRGB authoring values, not sampled from photographs.
PALETTE=np.array([[.56,.54,.49],[.66,.62,.54],[.48,.45,.40],[.72,.69,.61],
                  [.40,.41,.39],[.58,.49,.40],[.77,.74,.66]],'f4')
PROB=np.array([.28,.20,.18,.12,.08,.09,.05])

def smooth(a,lo,hi):
 t=np.clip((a-lo)/(hi-lo),0,1);return t*t*(3-2*t)

def field(n,cells,seed):
 rng=np.random.default_rng(seed);shape=(cells,cells) if isinstance(cells,int) else cells
 q=rng.uniform(-1,1,shape).astype('f4')
 return zoom(q,(n/shape[0],n/shape[1]),order=3,mode='grid-wrap',grid_mode=True).astype('f4')

def hash01(x,y,seed):
 v=x.astype('u8')*np.uint64(374761393)+y.astype('u8')*np.uint64(668265263)+np.uint64(seed)
 v=(v^(v>>np.uint64(13)))*np.uint64(1274126177);v=v^(v>>np.uint64(16))
 return (v&np.uint64(0xffffff)).astype('f4')/16777216.

def environment(x,y,tile,seed):
 """Periodic low-energy sheetwash, not a giant cracked playa surface."""
 u=np.asarray(x)/tile;v=np.asarray(y)/tile;ph=(seed%71)/71*math.tau
 wave=v+.055*np.sin(math.tau*u+ph)+.018*np.sin(3*math.tau*u-ph)
 dist=np.abs(np.sin(math.pi*(2*wave+.17)))
 wash=np.exp(-(dist/.10)**2)
 patch=.5+.25*np.sin(math.tau*(u+v)+ph)+.18*np.sin(math.tau*(2*u-v)-ph)
 return wash,np.clip(patch,0,1)

def sand_field(n,tile,seed):
 """Sparse explicit submillimetre angular grains seated within compacted fines."""
 cells=max(3,round(tile/.00082));pitch=tile/cells
 a=(np.arange(n,dtype='f4')+.5)/n*tile;idx=np.floor(a/pitch).astype('i4')
 height=np.zeros((n,n),'f4');cover=height.copy();tone=height.copy()
 # Compact grains fit their own jittered cells; 9-neighbour raster evaluation.
 for oy in (-1,0,1):
  cy=idx[:,None]+oy
  for ox in (-1,0,1):
   cx=idx[None,:]+ox;r=hash01(cx%cells,cy%cells,seed);q=hash01(cx%cells,cy%cells,seed+19)
   k=hash01(cx%cells,cy%cells,seed+49)
   dx=a[None,:]-(cx+.5+(r-.5)*.90)*pitch;dy=a[:,None]-(cy+.5+(q-.5)*.90)*pitch
   co=np.cos(k*math.tau);si=np.sin(k*math.tau)
   u=dx*co+dy*si;v=(-dx*si+dy*co)/(.62+.34*r)
   radius=.00012+.00028*q
   rho=np.maximum(np.abs(u),np.maximum(np.abs(.48*u+.876*v),np.abs(.48*u-.876*v)))
   d=radius-rho;m=smooth(d,-tile/n*.45,tile/n*.45)*(k>.42)
   h=np.maximum(0,np.minimum(d*.65,.00012+.00012*r))*m
   win=h>height;cover=np.maximum(cover,m);height=np.maximum(height,h);tone=np.where(win,k,tone)
 return height,cover,tone

def crust_field(n,tile,seed):
 """Weak, incomplete platy crust boundaries from a periodic irregular cell field."""
 cells=max(3,round(tile/.035));a=(np.arange(n,dtype='f4')+.5)/n*cells
 idx=np.floor(a).astype('i4');first=np.full((n,n),1e5,'f4');second=first.copy();label=np.zeros((n,n),'f4')
 for oy in (-1,0,1):
  cy=idx[:,None]+oy
  for ox in (-1,0,1):
   cx=idx[None,:]+ox;r=hash01(cx%cells,cy%cells,seed);q=hash01(cx%cells,cy%cells,seed+11)
   dx=a[None,:]-(cx+.5+(r-.5)*.85);dy=a[:,None]-(cy+.5+(q-.5)*.85)
   d=dx*dx+dy*dy;win=d<first;second=np.where(win,first,np.minimum(second,d));first=np.minimum(first,d);label=np.where(win,r,label)
 gap=(np.sqrt(second)-np.sqrt(first))*tile/cells
 crack=np.exp(-(gap/.00048)**2)*smooth(label,.24,.62)
 return crack,label

def scatter_grains(n,s):
 """Physical-size, toroidally packed angular gravel. No repeated cell mosaics.
 Conservative bounding-disk exclusion is approximate packing, not DEM.
 The exact accepted primitives are independent of raster resolution.
 """
 tile=s['tile_m'];rng=np.random.default_rng(s['seed']+401);pixel=tile/n
 m=np.zeros((n,n),'f4');h=m.copy();rough=m.copy();dust=m.copy();micro_h=m.copy();ids=np.zeros((n,n),'i4');rgb=np.zeros((n,n,3),'f4')
 # Number densities per square metre; coarse gravel islands plus granules.
 bands=[(.014,.038,480),(.005,.014,3600),(.0016,.005,22000)]
 bin_size=.024;nb=math.ceil(tile/bin_size);bin_size=tile/nb;bins={};records=[];primitives=[];counts=[];serial=0
 for low,high,density in bands:
  desired=round(density*tile*tile*s['gravel_density']);accepted=0;trials=0
  for trial in range(desired*18):
   if accepted>=desired:break
   trials+=1;cx,cy=rng.uniform(0,tile,2);radius=math.exp(rng.uniform(math.log(low*.5),math.log(high*.5)))
   wash,patch=environment(cx,cy,tile,s['seed']);prob=(.45+.55*patch)*(1-.70*wash)
   if rng.random()>prob:continue
   bx,by=int(cx/bin_size),int(cy/bin_size);reach=math.ceil((radius+.022)/bin_size);overlap=False
   for ix in range(bx-reach,bx+reach+1):
    for iy in range(by-reach,by+reach+1):
     for x0,y0,r0 in bins.get((ix%nb,iy%nb),[]):
      dx=abs(cx-x0);dy=abs(cy-y0);dx=min(dx,tile-dx);dy=min(dy,tile-dy)
      if dx*dx+dy*dy < (.80*(radius+r0))**2:overlap=True;break
     if overlap:break
    if overlap:break
   if overlap:continue
   bins.setdefault((bx,by),[]).append((cx,cy,radius));accepted+=1;serial+=1
   k=int(rng.integers(5,9));phase=rng.uniform(0,math.tau);ang=np.arange(k)*math.tau/k+phase+rng.uniform(-.12,.12,k)
   radii=radius*rng.uniform(.79,1.10,k);aspect=rng.uniform(.58,1.0);rot=rng.uniform(0,math.tau)
   vx=np.cos(ang)*radii;vy=np.sin(ang)*radii*aspect;c,ss=math.cos(rot),math.sin(rot);vx,vy=vx*c-vy*ss,vx*ss+vy*c
   ex=radius*1.25;xi=np.arange(math.floor((cx-ex)/pixel),math.ceil((cx+ex)/pixel)+1);yi=np.arange(math.floor((cy-ex)/pixel),math.ceil((cy+ex)/pixel)+1)
   dx=(xi[None,:]+.5)*pixel-cx;dy=(yi[:,None]+.5)*pixel-cy;d=np.full((len(yi),len(xi)),1,'f4')
   for j in range(k):
    px=vx[(j+1)%k]-vx[j];py=vy[(j+1)%k]-vy[j]
    d=np.minimum(d,(px*(dy-vy[j])-py*(dx-vx[j]))/math.hypot(px,py))
   # Ragged chipped outlines plus a small number of independent fracture planes.
   d+=radius*.026*np.sin(dx/radius*23+phase)*np.sin(dy/radius*27-phase)
   coverage=smooth(d,-pixel*.45,pixel*.45)
   slopes=rng.uniform(-.5,.5,(4,2))
   facets=np.minimum.reduce([np.broadcast_to(1+sa*dx/radius+sb*dy/radius,d.shape) for sa,sb in slopes])
   peak=radius*rng.uniform(.19,.43);burial=rng.uniform(.14,.38)
   body=np.minimum(np.maximum(d/(radius*.33),0),np.clip(facets,.20,1.2))
   height=np.maximum(0,body-burial)*peak*coverage
   species=int(rng.choice(len(PALETTE),p=PROB));col=PALETTE[species]*rng.uniform(.94,1.06)
   # Dust accumulates at the seated edge; weathered exposed faces warm slightly.
   coat=np.clip(s['dust_amount']*(.40+.50*np.exp(-np.maximum(d,0)/(radius*.18)))+.10*wash,0,.85)
   mat=np.array(s['color'],'f4');patchrgb=col[None,None,:]*(1-coat[...,None])+mat[None,None,:]*coat[...,None]
   # Mineral zoning is anisotropic and shared with a tiny fractured surface term.
   grain_pitch=max(.00045,radius*.075)
   gx=np.floor((dx+radius*2)/grain_pitch).astype('i4');gy=np.floor((dy+radius*2)/grain_pitch).astype('i4')
   tx=(dx+radius*2)/grain_pitch-gx;ty=(dy+radius*2)/grain_pitch-gy;tx=tx*tx*(3-2*tx);ty=ty*ty*(3-2*ty)
   crystal=(hash01(gx,gy,s['seed']+serial*7)*(1-tx)+hash01(gx+1,gy,s['seed']+serial*7)*tx)*(1-ty)+(hash01(gx,gy+1,s['seed']+serial*7)*(1-tx)+hash01(gx+1,gy+1,s['seed']+serial*7)*tx)*ty
   small=smooth(crystal,.66,.83);dark=1-smooth(crystal,.13,.26)
   patchrgb=patchrgb*(1-dark[...,None]*.26)+small[...,None]*.055
   # Sub-grain fracture is irregular, never a sinusoidal checkerboard.
   mineral_relief=(small-dark)*.000030*coverage
   height=np.maximum(0,height+mineral_relief)
   target=np.ix_(yi%n,xi%n);win=(coverage>.001)&(height>=h[target])
   h[target]=np.where(win,height,h[target]);m[target]=np.maximum(m[target],coverage)
   rough[target]=np.where(win,np.clip(rng.uniform(.66,.82)+coat*.11-.065*dark+.025*small,.55,.95),rough[target]);micro_h[target]=np.where(win,mineral_relief,micro_h[target]);dust[target]=np.where(win,coat,dust[target])
   ids[target]=np.where(win,serial,ids[target])
   for ch in range(3):rgb[...,ch][target]=np.where(win,patchrgb[...,ch],rgb[...,ch][target])
   records.append([float(cx),float(cy),float(radius),float(peak),species])
   primitives.append(dict(center=[float(cx),float(cy)],radius=float(radius),outline=np.stack((vx,vy),axis=-1).tolist(),roof_slopes=slopes.tolist(),peak=float(peak),burial=float(burial)))
  counts.append(dict(diameter_m=[low,high],target_count=desired,count=accepted,trials=trials))
 return m,h,rgb,rough,dust,ids,micro_h,dict(bands=counts,primitive_sha256=hashlib.sha256(np.asarray(records,'f8').tobytes()).hexdigest(),stone_count=serial,geometry_primitives=primitives)

def dirt_clods(n,tile,seed,density):
 """Finite soil agglomerates; irregular low relief, never a second stone palette."""
 rng=np.random.default_rng(seed);height=np.zeros((n,n),'f4');mask=height.copy();pixel=tile/n
 count=round(tile*tile*density)
 for _ in range(count):
  cx,cy=rng.uniform(0,tile,2);radius=rng.uniform(.0018,.0048);aspect=rng.uniform(.55,1);phase=rng.uniform(0,math.tau)
  x=np.arange(math.floor((cx-radius*1.3)/pixel),math.ceil((cx+radius*1.3)/pixel)+1);y=np.arange(math.floor((cy-radius*1.3)/pixel),math.ceil((cy+radius*1.3)/pixel)+1)
  dx=((x[None,:]+.5)*pixel-cx)/radius;dy=((y[:,None]+.5)*pixel-cy)/radius/aspect
  angle=np.arctan2(dy,dx);rho=(dx*dx+dy*dy)*(1+.18*np.sin(angle*3+phase)+.09*np.sin(angle*7-phase))
  z=np.maximum(0,1-rho)**1.2*rng.uniform(.00016,.00055);m=smooth(1-rho,0,.35)
  ix=np.ix_(y%n,x%n);height[ix]=np.maximum(height[ix],z);mask[ix]=np.maximum(mask[ix],m)
 return height,mask,count

def evaluate(spec=None,n=4096):
 s=DEFAULT.copy();s.update(spec or {})
 if not isinstance(n,int) or not 16<=n<=8192:raise ValueError('Resolution must be an integer 16..8192')
 if not .10<=s['tile_m']<=2:raise ValueError('tile_m must be 0.1..2 metres')
 if not .0<s['gravel_density']<=2:raise ValueError('gravel_density must be >0..2')
 if not .016<=s['height_scale_m']<=.1:raise ValueError('height_scale_m must be .016..0.1')
 if not 0<=s['sediment_depth_m']<=.010 or not 0<=s['sand_deposit_m']<=.006:raise ValueError('Sediment depths outside supported physical range')
 if not 0<=s['clod_density']<=10000:raise ValueError('Unsupported clod density')
 tile=s['tile_m'];seed=s['seed'];a=(np.arange(n,dtype='f4')+.5)/n*tile
 wash,deposit=environment(a[None,:],a[:,None],tile,seed)
 crack,crust=crust_field(n,tile,seed+97)
 crack*=smooth(deposit,.48,.80)
 meso=field(n,max(4,round(tile/.006)),seed+31);micro=field(n,max(4,round(tile/.00075)),seed+32)
 sand_h,sand_m,sand_t=sand_field(n,tile,seed+109)
 # Sheetwash lowers relief and deposits paler fines; crust relief and pigment share state.
 compact=.00025*meso+.000045*micro-.00125*wash+.00028*(crust-.5)
 compact-=crack*.00045*s['crust_strength']*(1-.7*wash)
 ground_h=compact+sand_h*(1-.7*wash)*(.32+.56*deposit)
 soil=np.empty((n,n,3),'f4');tint=np.array(s['color'],'f4')
 for ch in range(3):soil[...,ch]=tint[ch]+.025*deposit+.020*wash-.014*crack*s['crust_strength']+.004*meso
 sand_rgb=np.stack([.55+.10*sand_t,.51+.09*sand_t,.44+.08*sand_t],axis=-1)
 sand_visibility=sand_m*(.24+.32*sand_t)*(1-.4*wash)
 soil=soil*(1-sand_visibility[...,None])+sand_rgb*sand_visibility[...,None]
 m,h,rgb,r,dust,ids,micro_h,packing=scatter_grains(n,s)
 # A real sediment surface overtops existing stones instead of lifting them.
 # The foundation anchors gravel; soil-only geometry uses the deposited surface.
 support=ground_h.copy()
 sand_patch=smooth(field(n,max(3,round(tile/.075)),seed+811),-.30,.40)
 shelter=gaussian_filter((h>.0015).astype('f4'),max(.5,n/tile*.0035),mode='wrap')
 shelter=np.roll(shelter,round(n/tile*.006),axis=1)
 depth=s['sediment_depth_m']*(.47+.68*deposit)+s['sand_deposit_m']*sand_patch+.0008*wash+.0018*shelter
 clod_h,clod_m,clod_count=dirt_clods(n,tile,seed+922,s['clod_density'])
 # Loose sand is concentrated into deposits; compacted areas have subdued microrelief.
 deposited_base=gaussian_filter(support,max(.5,n/tile*s['compaction_smoothing_m']),mode='wrap')
 deposited_h=deposited_base+depth+clod_h*(.35+.65*(1-sand_patch))+sand_h*s['loose_sand_scale']*(.13+.87*sand_patch)
 above=h-(deposited_h-support)
 exposed=m*smooth(above,-.00008,.00020)
 full_h=np.maximum(deposited_h,support+h)
 soil_tint=np.array(s['color'],'f4');deposited_rgb=np.empty((n,n,3),'f4')
 for ch in range(3):deposited_rgb[...,ch]=soil_tint[ch]+.025*sand_patch+.018*deposit-.025*clod_m*(1-.6*sand_patch)+.0025*meso
 # Same buried/exposed state controls dust film, color and optical roughness.
 coating=np.clip(s['dust_amount']+.22*np.exp(-np.maximum(above,0)/.0013)+.08*sand_patch,.2,.985)
 coated_rock=rgb*(1-coating[...,None])+deposited_rgb*coating[...,None]
 base=deposited_rgb*(1-exposed[...,None])+coated_rock*exposed[...,None]
 rough=(.92+.025*clod_m-.025*sand_patch)*(1-exposed)+(r*(1-coating)+.92*coating)*exposed
 sediment_report={'initial_gravel_coverage':float((m>.5).mean()),'exposed_gravel_coverage':float((exposed>.5).mean()),'soil_sand_dust_coverage':float((exposed<=.5).mean()),'sediment_depth_m':[float(depth.min()),float(depth.mean()),float(depth.max())],'maximum_gravel_exposure_m':float(np.maximum(above,0).max()),'visible_gravel_bodies':int(len(np.unique(ids[exposed>.5]))),'dirt_clod_count':clod_count,'foundation_fixed':True}
 # Preserve primitive geometry and numerical ownership; expose only visible optics.
 ground_h=deposited_h;micro_h*=exposed*(1-coating)
 m=exposed

 normalized=.5+full_h/s['height_scale_m']
 if normalized.min()<0 or normalized.max()>1:raise ValueError('Physical height range exceeded; increase height_scale_m')
 fields={'BaseColor':np.clip(base,0,1).astype('f4'),'Height':normalized.astype('f4'),
  'Roughness':np.clip(rough,0,1).astype('f4'),'Metallic':np.zeros((n,n),'f4'),
  'Opacity':np.ones((n,n),'f4'),'AO':np.ones((n,n),'f4'),
  'SupportHeight':(.5+support/s['height_scale_m']).astype('f4'),'SedimentDepth':np.clip(depth/.01,0,1).astype('f4'),'DustCoating':coating.astype('f4'),'DirtClodMask':clod_m,'SandDepositMask':sand_patch,'SoilHeight':(.5+ground_h/s['height_scale_m']).astype('f4'),'GravelMicroHeight':(.5+micro_h/.00012).astype('f4'),'AggregateMask':m,'CrustMask':crack*(1-m),'WashMask':wash,'SandMask':sand_m*(1-m),'GrainID':ids}
 meta={**s,'resolution':n,'native_evaluation':True,'upscaled':False,'packing':packing,
  'sediment':sediment_report,'gravel_coverage':float((m>.5).mean()),'metric_height_minmax_m':[float(full_h.min()),float(full_h.max())],
  'height_encoding':'height_m = (Height - 0.5) * height_scale_m',
  'normal_contract':'Full metric object-space normal for static planar ground; not additive residual normal',
  'provenance':'Entirely numerical geometry/fields. No photographs, scans or generated-image inputs.',
  'limits':'Regional alluvial desert interpretation, not exact site, soil-series identification or calibrated optics. Weak patchy crust; no living biocrust claim.'}
 return fields,meta
