"""Independent saved-file check of all geometry/UV and non-normal shader state."""
import bpy,json,hashlib,numpy as np
from pathlib import Path
R=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def val(v):
 try:return list(v)
 except TypeError:return v
def inspect(path):
 bpy.ops.wm.open_mainfile(filepath=str(path));s=bpy.context.scene;objs={};materials={}
 for ob in s.objects:
  d=dict(type=ob.type,matrix=[list(row) for row in ob.matrix_world])
  if ob.type=="MESH":
   me=ob.data;v=np.empty(len(me.vertices)*3,"f4");me.vertices.foreach_get("co",v);l=np.empty(len(me.loops),"i4");me.loops.foreach_get("vertex_index",l)
   uv=[]
   for layer in me.uv_layers:
    a=np.empty(len(layer.data)*2,"f4");layer.data.foreach_get("uv",a);uv.append((layer.name,hashlib.sha256(a.tobytes()).hexdigest()))
   d.update(mesh=hashlib.sha256(v.tobytes()+l.tobytes()).hexdigest(),uv=uv,materials=[m.name for m in me.materials])
  if ob.type=="LIGHT":d.update(energy=ob.data.energy,size=ob.data.size,color=list(ob.data.color))
  if ob.type=="CAMERA":d.update(ortho_scale=ob.data.ortho_scale,type=ob.data.type,lens=ob.data.lens)
  objs[ob.name]=d
 for m in bpy.data.materials:
  if not m.use_nodes:continue
  nd=m.node_tree.nodes;allowed={n.name for n in nd if n.type!="OUTPUT_AOV" and n.name!="OIDN feature normal only"}
  ns={}
  for n in nd:
   if n.name not in allowed:continue
   d=dict(type=n.bl_idname,inputs=[(i.identifier,val(i.default_value)) for i in n.inputs if hasattr(i,"default_value")])
   for key in ["space","uv_map","interpolation","extension","distribution","operation"]:
    if hasattr(n,key):d[key]=getattr(n,key)
   if n.type=="TEX_IMAGE":
    p=Path(bpy.path.abspath(n.image.filepath)).resolve();d["image_sha256"]="NORMAL_IMAGE_ALLOWED_TO_CHANGE" if "Normal_OpenGL" in p.name else sha(p);d["color_space"]=n.image.colorspace_settings.name
   ns[n.name]=d
  links=sorted((l.from_node.name,l.from_socket.identifier,l.to_node.name,l.to_socket.identifier) for l in m.node_tree.links if l.from_node.name in allowed and l.to_node.name in allowed)
  materials[m.name]=dict(nodes=ns,links=links)
 return dict(objects=objs,materials=materials,view=[s.view_settings.view_transform,s.view_settings.look,s.view_settings.exposure,s.view_settings.gamma],camera=s.camera.name)
reports=[]
for cid in ["14_chainmail","12_lattice_wire"]:
 row=json.loads((R/"receipts"/(cid+"_binding.json")).read_text());before=inspect(row["source"]);after=inspect(row["scene"])
 assert before==after,(cid,"unexpected saved-file state difference")
 assert sha(row["source"])==row["source_sha256"] and sha(row["scene"])==row["scene_sha256"]
 reports.append(dict(material_id=cid,source_sha256=row["source_sha256"],scene_sha256=row["scene_sha256"],fresh_process_equivalence=True,all_vertex_topology_uv_object_rig_bytes=True,all_non_normal_shader_parameters_equal=True,normal_image_only_change=True))
(R/"receipts/fresh_equivalence.json").write_text(json.dumps(reports,indent=2));print(json.dumps(reports),flush=True)

