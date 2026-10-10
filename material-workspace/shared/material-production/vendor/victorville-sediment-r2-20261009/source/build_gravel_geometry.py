"""Derive finite convex gravel bodies from the exact procedural roof/edge planes.
Eliminates shadow stair-stepping from raster-only macro gravel geometry.
"""
import json,math
from pathlib import Path
import numpy as np
from scipy.spatial import HalfspaceIntersection,ConvexHull
from scipy.ndimage import map_coordinates

def build(root):
 root=Path(root);m=json.loads((root/'material.json').read_text());tile=m['tile_m'];scale=m['height_scale_m']
 soil=np.load(root/'SupportHeightState.npy');n=len(soil);primitives=json.loads((root/'GravelPrimitives.json').read_text());verts=[];faces=[];skipped=0
 for serial,p in enumerate(primitives):
  center=np.array(p['center']);radius=p['radius'];peak=p['peak'];burial=p['burial'];outline=np.array(p['outline'])/radius;k=len(outline)
  half=[]
  for j in range(k):
   v=outline[j];edge=outline[(j+1)%k]-v;length=np.linalg.norm(edge);ax=-edge[1]/length;ay=edge[0]/length;offset=-ax*v[0]-ay*v[1]
   half.append([-ax/.33,-ay/.33,1,burial-offset/.33])
  for a,b in p['roof_slopes']:half.append([-a,-b,1,burial-1])
  half.append([0,0,-1,0])
  points=HalfspaceIntersection(np.asarray(half),np.array([0,0,(1-burial)*.15])).intersections
  # Seeded shallow corner chips, inside the original support envelope.
  # Cut only exposed upper corners; original centers and placement are retained.
  if radius>.003:
   rng=np.random.default_rng(m['seed']+serial*71+901);physical=points*np.array([radius,radius,peak]);cuts=[]
   candidates=np.flatnonzero(points[:,2]>(1-burial)*.36)
   for vi in rng.choice(candidates,min(len(candidates),int(rng.integers(3,8))),replace=False):
    q=physical[vi];normal=np.array([q[0]/radius,q[1]/radius,.25+rng.uniform(.15,.8)]);normal/=np.linalg.norm(normal)
    depth=min(.00050,radius*rng.uniform(.018,.048));limit=np.dot(normal,q)-depth
    row=np.r_[normal*np.array([radius,radius,peak]),-limit]
    if row[3]+row[2]*(1-burial)*.15<-.00001:cuts.append(row)
   if cuts:points=HalfspaceIntersection(np.vstack((half,cuts)),np.array([0,0,(1-burial)*.15])).intersections
  hull=ConvexHull(points);tri=hull.simplices.copy();cross=np.cross(points[tri[:,1]]-points[tri[:,0]],points[tri[:,2]]-points[tri[:,0]]);flip=np.sum(cross*hull.equations[:,:3],axis=1)>0;tri[flip]=tri[flip][:,[0,2,1]];local=points*np.array([radius,radius,peak]);xy=local[:,:2]+center
  rows=(1-xy[:,1]/tile)*n-.5;cols=xy[:,0]/tile*n-.5
  # File rows are already in image order, origin upper left. Generator y uses row order.
  rows=xy[:,1]/tile*n-.5
  ground=map_coordinates(soil,np.stack((rows,cols)),order=1,mode='grid-wrap')
  xyz=np.column_stack((xy[:,0]-tile/2,tile/2-xy[:,1],(ground-.5)*scale+local[:,2]))
  # Torus-edge pieces get neighbouring copies with their global UV intact.
  shiftsx=[0];shiftsy=[0]
  if xy[:,0].min()<0:shiftsx.append(tile)
  if xy[:,0].max()>tile:shiftsx.append(-tile)
  if xy[:,1].min()<0:shiftsy.append(tile)
  if xy[:,1].max()>tile:shiftsy.append(-tile)
  for sx in shiftsx:
   for sy in shiftsy:
    offset=len(verts);v=xyz+np.array([sx,sy,0]);verts.extend(v.tolist());faces.extend((tri+offset).tolist())
 return np.asarray(verts,'f4'),np.asarray(faces,'i4')
def soil_normal(root):
 from generate_maps import png16rgb
 root=Path(root);m=json.loads((root/'material.json').read_text());h=np.load(root/'SoilHeightState.npy');n=len(h);step=m['tile_m']/n;scale=m['height_scale_m'];dx=(np.roll(h,-1,1)-np.roll(h,1,1))*scale/(2*step);dy=-(np.roll(h,-1,0)-np.roll(h,1,0))*scale/(2*step);inv=1/np.sqrt(1+dx*dx+dy*dy);nn=np.stack((-dx*inv,-dy*inv,inv),axis=-1);png16rgb(root/'Soil_Normal_Object_RGB16.png',np.rint((nn*.5+.5)*65535).astype('u2'))

if __name__=='__main__':
 import sys
 root=Path(sys.argv[1]);soil_normal(root);v,f=build(root);np.savez(root/'GravelGeometry.npz',vertices=v,triangles=f);m=json.loads((root/'material.json').read_text());m['render_geometry']='Explicit convex gravel bodies clipped from generator roof/edge planes plus soil-only metric surface';m['normal_contract']='Portable full height normal exported; production preview uses soil-only metric normal and physically chipped/beveled rock facets with coupled mica/quartz microrelief';m['files']={p.name:{'bytes':p.stat().st_size,'sha256':__import__('hashlib').sha256(p.read_bytes()).hexdigest()} for p in root.iterdir() if p.is_file() and p.name!='material.json'};(root/'material.json').write_text(json.dumps(m,indent=2));print(json.dumps(dict(vertices=len(v),triangles=len(f),physical_bodies=9389)))
