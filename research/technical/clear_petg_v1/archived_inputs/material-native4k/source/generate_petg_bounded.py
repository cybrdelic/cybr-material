"""Re-evaluate the selected seeded PETG surface recipe at native resolution, channel-sequential."""
import sys,json,hashlib,gc,argparse,resource
from pathlib import Path
import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter
SRC=Path('/workspace/shared/material-technical-rebuild/source');sys.path.insert(0,str(SRC))
from surface_math import field,clamp,smooth,color
from expansion.catalog import get
R=Path(__file__).resolve().parents[1]
def generate(cid,n):
 s=get(cid);out=(R/'sets' if n==4096 else R/'tests'/('reference_'+str(n)))/cid;out.mkdir(parents=True,exist_ok=True);maps={}
 def save(name,a,sixteen=False,space='Non-Color'):
  assert a.shape[:2]==(n,n) and np.isfinite(a).all()
  limits=[float(a.min()),float(a.max())];assert limits[0]>=0 and limits[1]<=1,(name,limits)
  p=out/(name+'.png');enc=np.empty(a.shape,np.uint16 if sixteen else np.uint8)
  for j in range(0,n,64):enc[j:j+64]=np.rint(a[j:j+64]*(65535 if sixteen else 255)).astype(enc.dtype)
  Image.fromarray(enc).save(p,compress_level=4)
  maps[name]={'file':p.name,'resolution':[n,n],'bit_depth':16 if sixteen else 8,'color_space':space,'range':limits,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()};print('MAP',cid,name,flush=True)
 fine=field(n,64,s['seed']+1);broad=field(n,5,s['seed']+2,2)
 rgb=clamp(color(s['color'],.005*broad)).astype(np.float32);save('BaseColor',rgb,space='sRGB');del rgb,broad;gc.collect()
 y=(np.arange(n,dtype=np.float32)[:,None]+.5)/n;phase=(y*s['tile_m']/s['layer_height_m'])%1;bead=np.sqrt(clamp(1-((phase-.5)/.52)**2));boundary=np.broadcast_to(1-smooth(np.minimum(phase,1-phase),.015,.12),(n,n)).copy();save('PrintBoundary',boundary)
 rough=clamp(s['roughness']+.14*boundary+.006*fine).astype(np.float32);save('Roughness',rough)
 orm=np.empty((n,n,3),np.float32);orm[:,:,0]=1;orm[:,:,1]=rough;orm[:,:,2]=s['metallic'];save('ORM',orm);del rough,orm,boundary;gc.collect()
 h=clamp(.22+.55*bead+.008*fine).astype(np.float32);del fine;save('Height',h,True);macro=gaussian_filter(h,max(.5,n/512),mode='wrap');save('Height_Macro',macro,True)
 def normal(name,height,invert_green=False):
  # Derivatives evaluated from full native scalar state, encoded in bounded row strips.
  enc=np.empty((n,n,3),np.uint8);scale=n*.5*s['height_scale_m']/s['tile_m']
  for y0 in range(0,n,64):
   y1=min(n,y0+64);rows=np.arange(y0,y1);a=height[y0:y1];dx=(np.roll(a,-1,1)-np.roll(a,1,1))*scale;dy=(height[(rows+1)%n]-height[(rows-1)%n])*scale
   v=np.stack((-dx,(-dy if invert_green else dy),np.ones_like(dx)),axis=-1);v/=np.linalg.norm(v,axis=-1,keepdims=True);enc[y0:y1]=np.rint((.5+.5*v)*255).astype(np.uint8)
  p=out/(name+'.png');Image.fromarray(enc).save(p,compress_level=4);maps[name]={'file':p.name,'resolution':[n,n],'bit_depth':8,'color_space':'Non-Color','range':[float(enc.min()/255),float(enc.max()/255)],'sha256':hashlib.sha256(p.read_bytes()).hexdigest()};print('MAP',cid,name,flush=True)
 normal('Normal_OpenGL',h);normal('Normal_DirectX',h,True);h-=macro;del macro;normal('Normal_Micro_OpenGL',h);del h;gc.collect()
 for name,value,sixteen in [('Metallic',s['metallic'],False),('AO',1,False),('Opacity',1,False),('Transmission',s.get('transmission',0),False),('Thickness',1,True)]:
  a=np.full((n,n),value,np.float32);save(name,a,sixteen);del a
 meta={**s,'resolution':[n,n],'provenance':'Native procedural re-evaluation of the selected legacy PETG generator; no image input and no upsampling','formation_claim':'Analytic fused-bead surface approximation with seeded finite-scale color/roughness variation; no thermal or extrusion solver in this production branch','source_generator':str(SRC/'expansion/maps.py'),'source_generator_sha256':hashlib.sha256((SRC/'expansion/maps.py').read_bytes()).hexdigest(),'export_generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'maps':maps,'height_midlevel':.5,'height_formula_m':'(Height - 0.5) * 0.00016','thickness_formula_m':'Thickness * 0.0012','normal_convention':'OpenGL +Y tangent normal; DirectX green-inverted alternative','basecolor_encoding':'sRGB bytes from original display-space authored palette; no baked lighting','packing':'ORM = AO,Roughness,Metallic','geometry_dependency':'Real wall/rim/seam thickness and bead profile remain in supplied geometry. Do not connect Height/Normal to existing resolved bead geometry; those channels are alternative flat-surface approximations. Constant AO is neutral, not a geometry-occlusion bake. Clear PETG requires true thickness, IOR1.57 and transmission.','usage':'Either use real bead geometry with BaseColor/Roughness/Metallic, or a flat-surface normal approximation. Do not double-count bead relief.','peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'material.json').write_text(json.dumps(meta,indent=2));print('COMPLETE',cid,len(maps),meta['peak_rss_kib'],flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('id');p.add_argument('--resolution',type=int,default=4096);a=p.parse_args();generate(a.id,a.resolution)
