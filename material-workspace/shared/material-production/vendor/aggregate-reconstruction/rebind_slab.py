"""Rebind selected immutable SI geometry to native4K PBR without rerendering."""
import bpy,sys,json,hashlib
from pathlib import Path
R=Path('/workspace/shared/material-aggregate-rebuild');N=R/'native4k';recipe=sys.argv[sys.argv.index('--')+1];row=next(x for x in json.loads((R/'provenance/CURRENT_CANDIDATES.json').read_text()) if x['recipe']==recipe);src=Path(row['scene']);meta_path=N/'maps'/recipe/'material_4k.json';assert meta_path.exists(),'Native set incomplete: no final manifest';native=json.loads(meta_path.read_text());assert native['resolution']==[4096,4096]
for channel,record in native['maps'].items():
 path=meta_path.parent/record['file'];assert path.exists(),path;hh=hashlib.sha256()
 with path.open('rb') as stream:
  for chunk in iter(lambda:stream.read(1024*1024),b''):hh.update(chunk)
 assert hh.hexdigest()==record['sha256'],('Incomplete/mixed native set',channel)
bpy.ops.wm.open_mainfile(filepath=str(src));scene=bpy.context.scene;print('SOURCE_OBJECTS',[(o.name,o.type,list(o.dimensions)) for o in scene.objects],flush=True)
# Remove only old studio support. Never include its30x floor in specimen bounds.
for o in list(scene.objects):
 if o.type in {'LIGHT','CAMERA'} or o.name.startswith(('Studio /','Capture /','Study / seamless ground')):bpy.data.objects.remove(o,do_unlink=True)
objects=[o for o in scene.objects if o.type in {'MESH','CURVE','CURVES','SURFACE'} and not o.hide_render];before={o.name:(len(o.data.vertices) if o.type=='MESH' else None,[list(v) for v in o.matrix_world]) for o in objects}
changed=[]
for mat in {m for o in objects if hasattr(o.data,'materials') for m in o.data.materials if m and m.use_nodes}:
 for node in mat.node_tree.nodes:
  if node.type!='TEX_IMAGE' or not node.image:continue
  stem=Path(node.image.filepath).stem;candidate=N/'maps'/recipe/(stem+'.png')
  # Selected rendering keeps16-bit normals;8-bit same-resolution copies are
  # separate compatibility/delivery images, never a hidden replacement.
  assert candidate.exists(),('Referenced source channel has no native counterpart',mat.name,stem)
  if candidate.exists():
   im=bpy.data.images.load(str(candidate),check_existing=True);im.colorspace_settings.name='sRGB' if stem=='BaseColor' else 'Non-Color';node.image=im;changed.append({'material':mat.name,'channel':stem,'map':str(candidate),'resolution':list(im.size)})
assert changed
for x in changed:assert x['resolution']==[4096,4096]
sys.path.insert(0,'/workspace/shared/material-slabs/source');from studio import configure
studio=configure(recipe,objects);after={o.name:(len(o.data.vertices) if o.type=='MESH' else None,[list(v) for v in o.matrix_world]) for o in objects};assert before==after
# Native maps stay external in this working slab to avoid9 redundant packed copies.
# The per-material archive carries their exact bytes and portable conventions.
out=N/'scenes'/(recipe+'_native4k_slab.blend');scene['native4k_specimen']=True;scene['native4k_maps_manifest']=str(N/'maps'/recipe/'material_4k.json');scene['source_selected_sha256']=row['sha256'];bpy.ops.wm.save_as_mainfile(filepath=str(out));report={'id':recipe,'selected_source':str(src),'selected_source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'scene':str(out),'scene_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'physical_tile_m':row['tile_m'],'geometry_and_transforms_unchanged':True,'maps':changed,'studio':studio,'images_external_for_bounded_disk':True,'rendered':False};(N/'receipts'/(recipe+'_slab.json')).write_text(json.dumps(report,indent=2)+'\n');print('NATIVE_SLAB_READY',json.dumps(report),flush=True)
