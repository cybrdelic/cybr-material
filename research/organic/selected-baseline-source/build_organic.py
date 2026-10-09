"""Editable structural replacements; authored procedural geometry, no photo inputs."""
import bpy,math,random,json,hashlib,sys
from pathlib import Path
import numpy as np
from mathutils import Vector
ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT/'inputs/CYBR_structures_portable_review/source'
sys.path.insert(0,str(OLD))
from expansion.geometry import mesh,plain,box,Tubes
from expansion.studies import aim,area

def configure(s,objects,width,view,lighting,size,samples,device):
 scene=bpy.context.scene;scene.render.engine="CYCLES";scene.cycles.samples=samples;scene.cycles.use_denoising=False
 scene.render.resolution_x=size;scene.render.resolution_y=size;scene.render.resolution_percentage=100
 scene.render.image_settings.file_format="PNG";scene.render.image_settings.color_mode="RGB";scene.render.image_settings.color_depth="16"
 scene.view_settings.view_transform="AgX";scene.view_settings.look="AgX - Medium High Contrast"
 scene.world=bpy.data.worlds.new("Capture / studio world");scene.world.use_nodes=True;scene.world.node_tree.nodes["Background"].inputs[0].default_value=(.18,.21,.25,1);scene.world.node_tree.nodes["Background"].inputs[1].default_value=.22
 scene.unit_settings.system="METRIC";scene.unit_settings.scale_length=1
 box("Studio / physical ground",(0,0,-.003),(width*30,width*30,.004),plain("Studio / warm neutral",(.15,.145,.13),.82))
 target=(0,0,.012 if s["family"]=="Grass" else .005)
 bpy.ops.object.camera_add();camera=bpy.context.object;camera.name="Capture / hero";scene.camera=camera;camera.data.type="ORTHO";camera.data.clip_start=.0001;camera.data.clip_end=100
 camera.location=(width*.86,-width*1.13,width);camera.data.ortho_scale=width*1.42;aim(camera,target)
 gain=width*width
 area("Capture / soft key",(-width*.7,-width*.65,width*1.4),target,width*1.05,36*gain,(1,.92,.82))
 area("Capture / grazing rim",(width*.55,width*.25,width*.60),target,width*.75,22*gain,(.80,.90,1))
 area("Capture / fill",(width*.8,-width*.55,width*.55),target,width*1.2,5*gain)
 scene["recipe_id"]=s["id"];scene["physical_scale_m"]=width;scene["capture_view"]="hero";scene["capture_lighting"]="soft"
 return scene

def material(name,rgb,rough=.85,sheen=0):
 m=plain(name,rgb,rough);p=m.node_tree.nodes.get('Principled BSDF');p.inputs['Sheen Weight'].default_value=sheen
 return m

def bark():
 a=np.load(ROOT/'data/bark_r4.npz');h=a['height'];c=a['color'];n=h.shape[0];w=.18
 y,x=np.mgrid[0:n,0:n];v=np.stack(((x/(n-1)-.5)*w,(y/(n-1)-.5)*w,h),-1).reshape(-1,3).tolist()
 inds=np.arange(n*n).reshape(n,n);q=inds[:-1,:-1].reshape(-1)
 f=np.stack((q,q+1,q+n+1,q+n),-1).tolist()
 m=material('Bark / weathered fractured cork',(.2,.15,.1),.94)
 nodes=m.node_tree.nodes;links=m.node_tree.links;p=nodes.get('Principled BSDF');attr=nodes.new('ShaderNodeVertexColor');attr.layer_name='Cork color';links.new(attr.outputs['Color'],p.inputs['Base Color'])
 tex=nodes.new('ShaderNodeTexNoise');tex.inputs['Scale'].default_value=3600;tex.inputs['Roughness'].default_value=.8
 bump=nodes.new('ShaderNodeBump');bump.inputs['Strength'].default_value=.28;bump.inputs['Distance'].default_value=.00017;links.new(tex.outputs['Fac'],bump.inputs['Height']);links.new(bump.outputs['Normal'],p.inputs['Normal'])
 o=mesh('Bark / branched craggy ridge islands',v,f,m);ca=o.data.color_attributes.new(name='Cork color',type='FLOAT_COLOR',domain='POINT');rgba=np.concatenate((c.reshape(-1,3),np.ones((n*n,1))),1).astype('f4');ca.data.foreach_set('color',rgba.ravel())
 for p in o.data.polygons:p.use_smooth=True
 # Separate dark cortex under the displaced surface creates a real cut section.
 under=box('Bark / dark inner cortex',(0,0,.0005),(w,w,.001),material('Bark / inner cortex',(.06,.036,.020),.98))
 sol=o.modifiers.new('Bark cut edge thickness','SOLIDIFY');sol.thickness=.0008
 o['primary_ridge_islands']=125;o['secondary_fracture_cells']=560;o['structure']='nonuniform junction network, interrupted checks, brittle multiscale relief'
 return [o,under],.216,{'id':'11_bark','name':'Mature fir bark r4','family':'Bark','seed':604},'r4'

def bark5():
 a=np.load(ROOT/'data/bark_r5.npz');h=a['height'];c=a['color'];n=h.shape[0];w=.18
 y,x=np.mgrid[0:n,0:n];v=np.stack(((x/(n-1)-.5)*w,(y/(n-1)-.5)*w,h),-1).reshape(-1,3).tolist()
 inds=np.arange(n*n).reshape(n,n);q=inds[:-1,:-1].reshape(-1)
 f=np.stack((q,q+1,q+n+1,q+n),-1).tolist()
 m=material('Bark / weathered fractured cork',(.2,.15,.1),.94)
 nodes=m.node_tree.nodes;links=m.node_tree.links;p=nodes.get('Principled BSDF');attr=nodes.new('ShaderNodeVertexColor');attr.layer_name='Cork color';links.new(attr.outputs['Color'],p.inputs['Base Color'])
 tex=nodes.new('ShaderNodeTexNoise');tex.inputs['Scale'].default_value=3600;tex.inputs['Roughness'].default_value=.8
 bump=nodes.new('ShaderNodeBump');bump.inputs['Strength'].default_value=.18;bump.inputs['Distance'].default_value=.00007;links.new(tex.outputs['Fac'],bump.inputs['Height']);links.new(bump.outputs['Normal'],p.inputs['Normal'])
 o=mesh('Bark / branched craggy ridge islands',v,f,m);ca=o.data.color_attributes.new(name='Cork color',type='FLOAT_COLOR',domain='POINT');rgba=np.concatenate((c.reshape(-1,3),np.ones((n*n,1))),1).astype('f4');ca.data.foreach_set('color',rgba.ravel())
 for p in o.data.polygons:p.use_smooth=True
 # Separate dark cortex under the displaced surface creates a real cut section.
 under=box('Bark / dark inner cortex',(0,0,.0005),(w,w,.001),material('Bark / inner cortex',(.06,.036,.020),.98))
 sol=o.modifiers.new('Bark cut edge thickness','SOLIDIFY');sol.thickness=.0008
 o['primary_ridge_islands']=82;o['secondary_fracture_cells']=560;o['structure']='nonuniform junction network, interrupted checks, brittle multiscale relief'
 return [o,under],.216,{'id':'11_bark','name':'Mature fir bark r5','family':'Bark','seed':605},'r5'

def bark6():
 a=np.load(ROOT/'data/bark_r5.npz');h=a['height'];c=a['color'];n=h.shape[0];w=.18
 y,x=np.mgrid[0:n,0:n];v=np.stack(((x/(n-1)-.5)*w,(y/(n-1)-.5)*w,h),-1).reshape(-1,3).tolist()
 inds=np.arange(n*n).reshape(n,n);q=inds[:-1,:-1].reshape(-1)
 f=np.stack((q,q+1,q+n+1,q+n),-1).tolist()
 # Close the cut section down to a flat inner cortex: no suspended relief shell.
 boundary=list(range(n))+[j*n+n-1 for j in range(1,n)]+list(range(n*n-2,n*(n-1)-1,-1))+[j*n for j in range(n-2,0,-1)]
 lower=[]
 for aidx in boundary:
  lower.append(len(v));v.append((v[aidx][0],v[aidx][1],-.001))
 for k in range(len(boundary)):
  k2=(k+1)%len(boundary);f.append((boundary[k],lower[k],lower[k2],boundary[k2]))
 f.append(tuple(reversed(lower)))
 m=material('Bark / weathered fractured cork',(.2,.15,.1),.94)
 nodes=m.node_tree.nodes;links=m.node_tree.links;p=nodes.get('Principled BSDF');attr=nodes.new('ShaderNodeVertexColor');attr.layer_name='Cork color';links.new(attr.outputs['Color'],p.inputs['Base Color'])
 tex=nodes.new('ShaderNodeTexNoise');tex.inputs['Scale'].default_value=3600;tex.inputs['Roughness'].default_value=.8
 bump=nodes.new('ShaderNodeBump');bump.inputs['Strength'].default_value=.18;bump.inputs['Distance'].default_value=.00007;links.new(tex.outputs['Fac'],bump.inputs['Height']);links.new(bump.outputs['Normal'],p.inputs['Normal'])
 o=mesh('Bark / branched craggy ridge islands',v,f,m);ca=o.data.color_attributes.new(name='Cork color',type='FLOAT_COLOR',domain='POINT');rgba=np.concatenate((c.reshape(-1,3),np.ones((n*n,1))),1).astype('f4');rgba=np.concatenate((rgba,np.tile([.010,.006,.003,1],(len(boundary),1))),axis=0);ca.data.foreach_set('color',rgba.ravel())
 for p in o.data.polygons:p.use_smooth=p.index<(n-1)**2
 o['primary_ridge_islands']=82;o['secondary_fracture_cells']=560;o['structure']='nonuniform junction network, interrupted checks, brittle multiscale relief'
 return [o],.216,{'id':'11_bark','name':'Mature fir bark r6','family':'Bark','seed':605},'r6'

def bark7():
 objects,width,settings,revision=bark6();rng=random.Random(607)
 a=np.load(ROOT/'data/bark_r5.npz');h=a['height'];c=a['color'];n=h.shape[0];w=.18
 dy,dx=np.gradient(h,w/(n-1),w/(n-1));verts=[];faces=[];colors=[];flakes=0
 def sample(arr,x,y):
  i=max(0,min(n-1,round((x/w+.5)*(n-1))));j=max(0,min(n-1,round((y/w+.5)*(n-1))))
  return arr[j,i]
 for trial in range(15000):
  if flakes>=2200:break
  x=rng.uniform(-w*.49,w*.49);y=rng.uniform(-w*.49,w*.49);sx=float(sample(dx,x,y));sy=float(sample(dy,x,y));hc=float(sample(h,x,y))
  if abs(sx)+abs(sy)>.65 or hc<.0040:continue
  # Fragments occupy crowns only, and remain much smaller than main ridge islands.
  length=rng.uniform(.0017,.0065);breadth=length*rng.uniform(.22,.58);yaw=rng.gauss(0,.28)
  sides=rng.choice((6,7,8));q=len(verts);outline=[];lift=rng.uniform(.00004,.00024);peel=rng.uniform(.00003,.00035)
  col=np.array(sample(c,x,y))*rng.uniform(.70,1.28)
  for k in range(sides):
   angle=math.tau*k/sides+rng.uniform(-.16,.16);radius=rng.uniform(.72,1.08)
   px=math.cos(angle)*breadth*.5*radius;py=math.sin(angle)*length*.5*radius
   xx=x+px*math.cos(yaw)-py*math.sin(yaw);yy=y+px*math.sin(yaw)+py*math.cos(yaw)
   zz=max(float(sample(h,xx,yy))+.000022,hc+sx*(xx-x)+sy*(yy-y)+lift+peel*(py/length+.5))
   outline.append((xx,yy,zz));verts.append((xx,yy,zz));colors.append((*col,1))
  for xx,yy,zz in outline:
   verts.append((xx,yy,float(sample(h,xx,yy))-.000015));colors.append((col[0]*.46,col[1]*.40,col[2]*.35,1))
  faces.append(tuple(q+k for k in range(sides)))
  for k in range(sides):k2=(k+1)%sides;faces.append((q+k,q+sides+k,q+sides+k2,q+k2))
  flakes+=1
 m=material('Bark / brittle exposed lamellae',(.05,.03,.014),.95)
 node=m.node_tree.nodes.new('ShaderNodeVertexColor');node.layer_name='Fracture color';m.node_tree.links.new(node.outputs['Color'],m.node_tree.nodes.get('Principled BSDF').inputs['Base Color'])
 ob=mesh('Bark / irregular small fractured cork lamellae',verts,faces,m);ca=ob.data.color_attributes.new(name='Fracture color',type='FLOAT_COLOR',domain='POINT');ca.data.foreach_set('color',np.array(colors,dtype='f4').ravel())
 ob['fragment_count']=flakes;ob['fragment_length_range_m']=[.0017,.0065];ob['construction']='random crown-only cork fragments, no rows and no repeating shingle lattice';objects.append(ob)
 settings.update(name='Mature fir bark r7',seed=607)
 return objects,width,settings,'r7'

def bark8():
 objects,width,settings,revision=bark7()
 for m in [bpy.data.materials.get('Bark / weathered fractured cork'),bpy.data.materials.get('Bark / brittle exposed lamellae')]:
  nodes=m.node_tree.nodes;links=m.node_tree.links;p=nodes.get('Principled BSDF')
  source=p.inputs['Base Color'].links[0].from_socket
  texco=nodes.new('ShaderNodeTexCoord');stretch=nodes.new('ShaderNodeVectorMath');stretch.operation='MULTIPLY';stretch.inputs[1].default_value=(1,.42,1)
  links.new(texco.outputs['Object'],stretch.inputs[0])
  grain=nodes.new('ShaderNodeTexNoise');grain.inputs['Scale'].default_value=2900;grain.inputs['Detail'].default_value=3.2;grain.inputs['Roughness'].default_value=.78
  links.new(stretch.outputs[0],grain.inputs['Vector'])
  bump=nodes.new('ShaderNodeBump');bump.inputs['Strength'].default_value=.60;bump.inputs['Distance'].default_value=.00020;links.new(grain.outputs['Fac'],bump.inputs['Height']);links.new(bump.outputs['Normal'],p.inputs['Normal'])
  tone=nodes.new('ShaderNodeTexNoise');tone.inputs['Scale'].default_value=330;tone.inputs['Detail'].default_value=3;links.new(stretch.outputs[0],tone.inputs['Vector'])
  ramp=nodes.new('ShaderNodeValToRGB');ramp.color_ramp.elements[0].position=.18;ramp.color_ramp.elements[0].color=(.54,.49,.44,1);ramp.color_ramp.elements[1].position=.82;ramp.color_ramp.elements[1].color=(1.12,1.10,1.06,1);links.new(tone.outputs['Fac'],ramp.inputs[0])
  mul=nodes.new('ShaderNodeMixRGB');mul.blend_type='MULTIPLY';mul.inputs[0].default_value=1;links.new(source,mul.inputs[1]);links.new(ramp.outputs[0],mul.inputs[2]);links.new(mul.outputs[0],p.inputs['Base Color'])
  rough=nodes.new('ShaderNodeMapRange');rough.inputs['From Min'].default_value=0;rough.inputs['From Max'].default_value=1;rough.inputs['To Min'].default_value=.82;rough.inputs['To Max'].default_value=.97;links.new(grain.outputs['Fac'],rough.inputs[0]);links.new(rough.outputs[0],p.inputs['Roughness'])
 fragment=objects[-1];bev=fragment.modifiers.new('Minute fractured edge rounding','BEVEL');bev.width=.000035;bev.segments=2
 settings.update(name='Mature fir bark r8',seed=607)
 return objects,width,settings,'r8'

sys.path.insert(0,str(Path(__file__).resolve().parent))
from organic_shapes import grass,carpet,fleece

def main():
 key=sys.argv[sys.argv.index('--')+1] if '--' in sys.argv else 'bark';bpy.ops.wm.read_factory_settings(use_empty=True)
 objects,width,s,rev=globals()[key]();scene=configure(s,objects,width,'hero','soft',960,64,'CPU')
 scene.render.engine='CYCLES';scene['organic_revision']=rev
 scene['procedural_only']=True;scene['specimen_construction']='Explicit geometry; source in organic rebuild bundle'
 # Keep original camera/light recipe for direct specimen comparison.
 blend=ROOT/'scenes'/f'{s["id"]}_{rev}_hero.blend';bpy.ops.wm.save_as_mainfile(filepath=str(blend))
 hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'source').glob('*.py')}
 receipt={'recipe':s['id'],'revision':rev,'scene':str(blend),'scene_sha256':hashlib.sha256(blend.read_bytes()).hexdigest(),'source_sha256':hashes,'rendered':False,'visual_acceptance':'pending','physical_capture_width_m':width,'objects':[{'name':o.name,'vertices':len(o.data.vertices) if o.type=='MESH' else 0,'faces':len(o.data.polygons) if o.type=='MESH' else 0} for o in objects]}
 (ROOT/'receipts'/f'{s["id"]}_{rev}_build.json').write_text(json.dumps(receipt,indent=2));print('READY',json.dumps(receipt),flush=True)
if __name__=='__main__':main()
