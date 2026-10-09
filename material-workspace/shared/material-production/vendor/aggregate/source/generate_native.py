"""Native4096 evaluation of the selected procedural recipe, not image enlargement."""
from pathlib import Path
import sys,json,time,gc,hashlib,struct,zlib,resource,shutil,math
import numpy as np
from PIL import Image
from scipy.ndimage import map_coordinates
ROOT=Path('/workspace/shared/material-aggregate-rebuild');OUT=ROOT/'native4k';sys.path.insert(0,str(ROOT/'source'));from mineral_construction_bounded import construction
from expansion.maps import normals
# Exact pointwise algebra in row blocks avoids allocating four simultaneous
# full-resolution float64 temporaries in board-form thresholding.
_original_smooth=construction.__globals__['smooth']
def bounded_smooth(a,lo,hi):
 if not isinstance(a,np.ndarray) or a.ndim!=2 or a.dtype!=np.float64 or a.size<1048576:return _original_smooth(a,lo,hi)
 out=np.empty_like(a);low=np.broadcast_to(lo,a.shape) if np.ndim(lo) else None;high=np.broadcast_to(hi,a.shape) if np.ndim(hi) else None
 for row in range(0,a.shape[0],64):out[row:row+64]=_original_smooth(a[row:row+64],low[row:row+64] if low is not None else lo,high[row:row+64] if high is not None else hi)
 return out
construction.__globals__['smooth']=bounded_smooth
while (OUT/'GENERATION_HOLD').exists():
 print('GENERATION_HELD_FOR_SHARED_MEMORY',flush=True);time.sleep(10)
recipe=sys.argv[1];NORMAL_BITS=int(sys.argv[2]) if len(sys.argv)>2 else 8;N=4096;src=ROOT/'expansion/maps'/recipe;dest=OUT/'maps'/recipe;dest.mkdir(parents=True,exist_ok=True);s=json.loads((src/'material.json').read_text());start=time.time();checks={}
# Preserve the exact generator snapshot and source recipe references for this output.
snapshot=OUT/'source'/recipe;snapshot.mkdir(parents=True,exist_ok=True)
for p in [ROOT/'source/expansion/mineral_construction.py',ROOT/'source/expansion/grain_model.py',ROOT/'source/surface_math.py',OUT/'source/mineral_construction_bounded.py']:
 shutil.copy2(p,snapshot/p.name)
shutil.copy2(src/'material.json',snapshot/'selected_recipe_1536.json')

def png16rgb(path,a):
 h,w,c=a.shape;assert c==3
 def chunk(f,tag,data):f.write(struct.pack('>I',len(data)));f.write(tag);f.write(data);f.write(struct.pack('>I',zlib.crc32(tag+data)&0xffffffff))
 with path.open('wb') as f:
  f.write(b'\x89PNG\r\n\x1a\n');chunk(f,b'IHDR',struct.pack('>IIBBBBB',w,h,16,2,0,0,0));compress=zlib.compressobj(4)
  for row in a:
   raw=np.ascontiguousarray(row.astype('>u2')).view('u1').ravel();filtered=np.empty_like(raw);filtered[:6]=raw[:6];np.subtract(raw[6:],raw[:-6],out=filtered[6:]);data=compress.compress(b'\x01'+filtered.tobytes())
   if data:chunk(f,b'IDAT',data)
  data=compress.flush()
  if data:chunk(f,b'IDAT',data)
  chunk(f,b'IEND',b'')
def save(name,a,bits=8,meaning=None):
 assert a.shape[:2]==(N,N) and np.isfinite(a).all(),name
 lo,hi=float(a.min()),float(a.max());assert lo>=-1e-5 and hi<=1.00001,(name,lo,hi)
 path=dest/(name+'.png');encoded=np.rint(np.clip(a,0,1)*(65535 if bits==16 else 255)).astype('uint16' if bits==16 else 'uint8')
 if bits==16 and a.ndim==3:png16rgb(path,encoded)
 else:Image.fromarray(encoded).save(path,compress_level=4)
 raw=path.read_bytes();assert raw[:8]==b'\x89PNG\r\n\x1a\n';width,height=struct.unpack('>II',raw[16:24]);assert (width,height,raw[24])==(N,N,bits)
 checks[name]={'file':path.name,'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw),'resolution':[N,N],'bits_per_channel':bits,'colorspace':'sRGB' if name=='BaseColor' else 'Non-Color / linear data','range':[lo,hi],'meaning':meaning};print('SAVED',recipe,name,bits,len(raw),flush=True);del encoded,raw;gc.collect()

print('NATIVE_GENERATION_START',recipe,N,flush=True);values,construction_report=construction(s,N);print('ANALYTIC_FIELDS_READY',recipe,time.time()-start,flush=True)
height=values.pop('Height');save('Height',height,16,'Full source surface height. displacement_m=(value-0.5)*height_scale_m')
for name in ['BaseColor','Roughness','AggregateMask']:
 a=values.pop(name);save(name,a,16 if name=='Roughness' else 8);del a;gc.collect()
# Constant physical channels are included explicitly rather than invented detail.
for name,constant in [('Metallic',0.),('Opacity',1.)]:
 a=values.pop(name);save(name,a,8,'Constant '+str(constant)+'; nonmetal/opaque material');del a;gc.collect()
values.clear();gc.collect()
normal=normals(height,s);save('Normal_OpenGL',normal,NORMAL_BITS,'Full tangent-space normal, OpenGL +Y. Use without macro displacement or without a second full normal.');normal[...,1]=1-normal[...,1];save('Normal_DirectX',normal,NORMAL_BITS,'Full tangent-space normal, DirectX -Y. Green channel inverted from OpenGL.');del normal;gc.collect()
# Native bake from the retained selected mesh, not enlargement of a1536 bitmap.
mesh_height=np.load(src/'GeometryHeight.npy') if (src/'GeometryHeight.npy').exists() else None
if mesh_height is not None:
 shutil.copy2(src/'GeometryHeight.npy',dest/'GeometryHeight.npy');grid=mesh_height.shape[0]-1;macro=np.empty((N,N),'f4');px=(np.arange(N,dtype='f4')+.5)*grid/N
 for row0 in range(0,N,128):
  py=(np.arange(row0,min(row0+128,N),dtype='f4')+.5)*grid/N;yy,xx=np.meshgrid(py,px,indexing='ij');macro[row0:row0+len(py)]=map_coordinates(mesh_height,np.stack((yy,xx)),order=1,mode='nearest')
else:
 # Terrazzo selected flat tiles do not use a geometric height field.
 macro=np.full((N,N),.5,'f4')
save('Height_Macro',macro,16,'Native rasterization of retained selected GeometryHeight mesh, for matched displacement/normal split');residual=height-macro;micro=normals(residual,s);save('Normal_Micro_OpenGL',micro,NORMAL_BITS,'Residual tangent normal against retained mesh height; selected slab uses this to avoid double relief');del macro,residual,micro;gc.collect()
# Optional diffuse horizon visibility: no light direction, color, or baked shadow.
h=(height-.5)*s['height_scale_m'];ao=np.zeros_like(height);steps_m=[.00025,.0007,.0017,.0045,.010]
for angle in np.arange(8)*math.tau/8:
 horizon=np.zeros_like(height)
 for dist in steps_m:
  sx=int(round(math.cos(angle)*dist/s['tile_m']*N));sy=int(round(math.sin(angle)*dist/s['tile_m']*N))
  if sx==sy==0:continue
  distance=math.hypot(sx,sy)*s['tile_m']/N;shifted=np.roll(np.roll(h,sy,0),sx,1);slope=(shifted-h)/distance
  if s.get('tiling')=='finite specimen':
   if sy>0:slope[:sy]=0
   if sy<0:slope[sy:]=0
   if sx>0:slope[:,:sx]=0
   if sx<0:slope[:,sx:]=0
  np.maximum(horizon,slope,out=horizon)
 ao+=1/(1+horizon*horizon)
 del horizon,shifted,slope;gc.collect()
ao/=8;save('AO',ao,8,'Optional8-direction diffuse horizon visibility from full physical height; not used to darken selected albedo');del ao,h,height;gc.collect()
metadata={'id':recipe,'selected_source_metadata':str(src/'material.json'),'resolution':[N,N],'native_evaluation':True,'bounded_memory_math':'Large float64 smooth uses exact row blocks; dead full-resolution intermediates released early. Physical recipe and arithmetic unchanged','upscaled_from_bitmap':False,'normal_bits':NORMAL_BITS,'exporter_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'generator':'Original selected procedural recipe evaluated directly at4096; retained mesh rasterized directly for residual split','generator_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in snapshot.glob('*.py')},'tile_m':s['tile_m'],'height_scale_m':s['height_scale_m'],'height_midlevel':.5,'displacement_equation':'(Height-0.5)*height_scale_m','normal_convention':'OpenGL +Y; DirectX -Y optional','selected_geometry_pair':'retained GeometryHeight.npy with Normal_Micro_OpenGL; full Normal_OpenGL for flat/no-displacement surfaces','metallic':0,'opacity':1,'construction':construction_report,'maps':checks,'quality':'Native resolution upgrade of current selected authored material. Does not claim causal simulation or cure rejected courtyard wall appearance.','seconds':time.time()-start,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
(dest/'material_4k.json').write_text(json.dumps(metadata,indent=2)+'\n');(OUT/'receipts'/(recipe+'_native4k.json')).write_text(json.dumps(metadata,indent=2)+'\n');print('NATIVE_4K_COMPLETE',recipe,metadata['seconds'],metadata['peak_rss_kib'],flush=True)
