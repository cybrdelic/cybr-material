import bpy,math,sys,json,hashlib,bmesh
from pathlib import Path
R=Path('/workspace/shared/material-slabs');NATIVE=Path('/workspace/shared/material-native4k');sys.path.insert(0,str(R/'source'));from studio import configure
cid=sys.argv[sys.argv.index('--')+1];SRC=Path('/workspace/shared/material-technical-rebuild/scenes')/(cid+'_r5_hero.blend');P=NATIVE/'sets'/cid;meta=json.load(open(P/'material.json'));assert meta['resolution']==[4096,4096]
bpy.ops.wm.read_factory_settings(use_empty=True)
with bpy.data.libraries.load(str(SRC),link=False) as (a,b):b.materials=[n for n in a.materials if n.startswith('PETG / ')][:1]
m=b.materials[0];assert m
for node in m.node_tree.nodes:
 if node.type=='TEX_IMAGE' and node.image:
  name=Path(node.image.filepath).name;im=bpy.data.images.load(str(P/name));assert list(im.size)==[4096,4096];im.colorspace_settings.name=node.image.colorspace_settings.name;node.image=im;im.pack()
# Unrolled continuation of selected cylindrical fused-bead wall; same bead formula and seam cross-section.
w=.028;h=.026;wall=.0012;pitch=.0002;beadwidth=.00045;nx=128;ny=round(h/pitch)*8;v=[];uv=[];f=[]
for inner in (False,True):
 for j in range(ny+1):
  y=j/ny*h;phase=(y/pitch)%1;bead=math.sqrt(max(0,1-((phase-.5)/.52)**2));relief=(bead-.65)*beadwidth*.155
  for i in range(nx+1):
   x=(i/nx-.5)*w;seam=0 if inner else .000055*math.exp(-.5*((x+.006)/(.014*.035))**2)*(.65+.35*bead)
   v.append((x,y-h/2,wall+relief+seam-(wall if inner else 0)));uv.append(((x+w/2)/meta['tile_m'],y/meta['tile_m']))
stride=(nx+1)*(ny+1)
for side in (0,1):
 for j in range(ny):
  for i in range(nx):
   q=side*stride+j*(nx+1)+i;face=(q,q+1,q+nx+2,q+nx+1);f.append(face if side==0 else tuple(reversed(face)))
bound=list(range(nx+1))+[j*(nx+1)+nx for j in range(1,ny+1)]+[ny*(nx+1)+i for i in range(nx-1,-1,-1)]+[j*(nx+1) for j in range(ny-1,0,-1)]
for a,b in zip(bound,bound[1:]+bound[:1]):f.append((a,a+stride,b+stride,b))
mesh=bpy.data.meshes.new(cid+' / closed unrolled bead-wall mesh');mesh.from_pydata(v,[],f);mesh.update();ob=bpy.data.objects.new(cid+' / true1.2mm printed panel',mesh);bpy.context.scene.collection.objects.link(ob);mesh.materials.append(m);layer=mesh.uv_layers.new(name='UVMap');layer.active_render=True
for poly in mesh.polygons:
 poly.use_smooth=poly.index<2*nx*ny
 for k in poly.loop_indices:layer.data[k].uv=uv[mesh.loops[k].vertex_index]
bm=bmesh.new();bm.from_mesh(mesh);assert all(len(e.link_faces)==2 for e in bm.edges);bm.free()
ob['layer_height_m']=pitch;ob['wall_thickness_m']=wall;ob['no_double_relief']=True;ob['source_bead_formula']='sqrt(max(0,1-((phase-.5)/.52)^2)); relief=(bead-.65)*.00045*.155; original55um start/stop seam';report=configure(cid,[ob]);report.update(source_scene=str(SRC),source_sha256=hashlib.sha256(SRC.read_bytes()).hexdigest(),maps_native_resolution=4096,layer_height_m=pitch,wall_thickness_m=wall,layer_count=round(h/pitch),vertices=len(v),faces=len(f),manifold=True,selected_appearance='r5 PETG shader and exact source bead/seam formula; flattened wall topology only',optical_limit='Continuous fused shell; internal voids omitted. Clear variant remains translucent/frosted and is not certified optically clear.',map_manifest=str(P/'material.json'))
s=bpy.context.scene;s.name='Slab / '+cid+' / native4K';out=R/'scenes'/(cid+'_native4k_slab.blend');bpy.ops.wm.save_as_mainfile(filepath=str(out),compress=True);report['scene']=str(out);report['scene_sha256']=hashlib.sha256(out.read_bytes()).hexdigest();(R/'receipts'/(cid+'_native4k_slab.json')).write_text(json.dumps(report,indent=2));print('PETG_SLAB_READY',report['scene_sha256'],flush=True)
