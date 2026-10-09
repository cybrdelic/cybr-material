from pathlib import Path
import sys,importlib.util,json,hashlib,time,struct,zlib,gc,resource
import numpy as np
from scipy.ndimage import gaussian_filter,map_coordinates
from PIL import Image
R=Path(__file__).resolve().parents[1];root=R.parents[1];B=R/'base_recipe';sys.path[:0]=[str(R/'source'),str(B/'vendor'),str(B/'source')]
from fractured_construction import construction
from height_normals import height_normals,encode16
s=json.loads((B/'material.json').read_text());N=4096;out=R/'native4096';out.mkdir(exist_ok=True);start=time.monotonic();v,report=construction(s,N)
files={}
for name in ['BaseColor','Roughness','Metallic']:
 p=out/(name+'.png');Image.fromarray(np.rint(np.clip(v[name],0,1)*255).astype('u1')).save(p,compress_level=4);files[name]=hashlib.sha256(p.read_bytes()).hexdigest()
h=v['Height'].copy();np.save(out/'HeightState.npy',h);del v;gc.collect();macro=gaussian_filter(h,N/512*1.15,mode='wrap');coord=np.linspace(-.5,N-.5,513,dtype='f4');y,x=np.meshgrid(coord,coord,indexing='ij');g=map_coordinates(macro,np.stack((y,x)),order=1,mode='grid-wrap');np.save(out/'GeometryHeight.npy',g.astype('f4'));del macro,x,y;gc.collect()
n=height_normals(h,tile_m=s['tile_m'],height_scale_m=s['height_scale_m'],periodic=True);a=encode16(n);del n;gc.collect()
def chunk(t,b):return struct.pack('>I',len(b))+t+b+struct.pack('>I',zlib.crc32(t+b)&0xffffffff)
p=out/'Normal_Object_RGB16.png';comp=zlib.compressobj(4);blocks=[]
for row in a:blocks.append(comp.compress(b'\0'+row.astype('>u2').tobytes()))
blocks.append(comp.flush());p.write_bytes(b'\x89PNG\r\n\x1a\n'+chunk(b'IHDR',struct.pack('>IIBBBBB',N,N,16,2,0,0,0))+chunk(b'IDAT',b''.join(blocks))+chunk(b'IEND',b''));files['Normal_Object']=hashlib.sha256(p.read_bytes()).hexdigest()
r={'native_resolution':[N,N],'upscaled':False,'source_seed':s['seed'],'construction':report,'fracture_model':s['fracture_model'],'files':files,'seconds':time.monotonic()-start,'peak_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'physical_qualification':False,'selected':False};(out/'receipt.json').write_text(json.dumps(r,indent=2));print(json.dumps(r))
