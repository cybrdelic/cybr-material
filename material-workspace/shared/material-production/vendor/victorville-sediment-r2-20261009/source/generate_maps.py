"""Native procedural evaluation and metric PBR export; no upsampling."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2');os.environ.setdefault('OMP_NUM_THREADS','2')
import argparse,gc,json,time,hashlib,struct,zlib,resource
from pathlib import Path
import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter,map_coordinates
from victorville_ground import evaluate,DEFAULT

def png16rgb(path,a):
 def chunk(t,b):return struct.pack('>I',len(b))+t+b+struct.pack('>I',zlib.crc32(t+b)&0xffffffff)
 n,m,_=a.shape;comp=zlib.compressobj(4);blocks=[]
 for row in a:blocks.append(comp.compress(b'\0'+row.astype('>u2').tobytes()))
 blocks.append(comp.flush());path.write_bytes(b'\x89PNG\r\n\x1a\n'+chunk(b'IHDR',struct.pack('>IIBBBBB',m,n,16,2,0,0,0))+chunk(b'IDAT',b''.join(blocks))+chunk(b'IEND',b''))

def export(out,n,s=None):
 out=Path(out);out.mkdir(parents=True,exist_ok=True);start=time.monotonic();v,meta=evaluate(s,n);print('Evaluated',n,'procedural state',flush=True)
 h=v.pop('Height');support=v.pop('SupportHeight');np.save(out/'SupportHeightState.npy',support);del support;soil=v.pop('SoilHeight');primitives=meta['packing'].pop('geometry_primitives');(out/'GravelPrimitives.json').write_text(json.dumps(primitives));np.save(out/'SoilHeightState.npy',soil);ids=v.pop('GrainID');np.save(out/'GrainID.npy',ids);del ids
 for name,a in v.items():Image.fromarray(np.rint(np.clip(a,0,1)*255).astype('u1')).save(out/(name+'.png'),compress_level=4)
 del v;gc.collect();Image.fromarray(np.rint(h*65535).astype('u2')).save(out/'Height.png');np.save(out/'HeightState.npy',h)
 grid=min(512,n//2);macro=gaussian_filter(h,max(.5,n/grid*.8),mode='wrap')
 coord=np.linspace(-.5,n-.5,grid+1,dtype='f4');y,x=np.meshgrid(coord,coord,indexing='ij')
 g=map_coordinates(macro,np.stack((y,x)),order=1,mode='grid-wrap').astype('f4');np.save(out/'GeometryHeight.npy',g);soil_macro=gaussian_filter(soil,max(.5,n/grid*.8),mode='wrap');gs=map_coordinates(soil_macro,np.stack((y,x)),order=1,mode='grid-wrap').astype('f4');np.save(out/'SoilGeometryHeight.npy',gs);del macro,soil,soil_macro,gs,x,y,g
 # PNG top row maps to Blender UV y=1. Green uses world +Y (negative array row).
 step=meta['tile_m']/n;scale=meta['height_scale_m'];dx=(np.roll(h,-1,1)-np.roll(h,1,1))*scale/(2*step);dy=-(np.roll(h,-1,0)-np.roll(h,1,0))*scale/(2*step)
 inv=1/np.sqrt(1+dx*dx+dy*dy);normal=np.stack((-dx*inv,-dy*inv,inv),axis=-1);a=np.rint((normal*.5+.5)*65535).astype('u2');png16rgb(out/'Normal_Object_RGB16.png',a);del a,normal,dx,dy,inv
 meta.update(geometry_grid=grid,map_contract={'BaseColor':'sRGB reflectance','Height':'linear UNORM16, midpoint0.5','Roughness':'linear','Metallic':'linear,zero','Normal_Object_RGB16':'linear RGB16 full object-space SI height normal'},seconds=time.monotonic()-start,peak_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
 meta['files']={p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.iterdir() if p.is_file()}
 (out/'material.json').write_text(json.dumps(meta,indent=2));print(json.dumps({k:meta[k] for k in ['gravel_coverage','metric_height_minmax_m','seconds','peak_rss_mib']}),flush=True)
 return meta
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--resolution',type=int,default=4096);p.add_argument('--out',required=True);p.add_argument('--parameters');args=p.parse_args();s=json.loads(Path(args.parameters).read_text()) if args.parameters else None;export(args.out,args.resolution,s)
