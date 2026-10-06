"""Continuity correction for the visible rounded glaze surface only.
No solved height, packet layout, geometry, shader coefficient or light changes.
The radial developed chart is exact on ideal cylindrical edges; rounded corners
have unavoidable metric distortion and are not an isometric coating-flow solve.
"""
import bpy,numpy as np,json,hashlib,math
from pathlib import Path
R=Path(__file__).resolve().parents[1];EXP=R.parent;WORK=EXP.parents[2];ref=json.loads((EXP/'porcelain_arrival_poisson/experiment.json').read_text());src=WORK/ref['source_scene']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(src)==ref['source_sha256'];bpy.ops.wm.open_mainfile(filepath=str(src));bpy.context.preferences.filepaths.save_version=0
s=bpy.context.scene;ob=s.objects['07 new / finite60micron mean glaze shell'];me=ob.data;uv=me.uv_layers.active;period=.0192;edge=.0006;ix=iy=.096/2-edge;zc=.0036/2-edge
verts=np.array([tuple(v.co) for v in me.vertices],dtype='f4');topology=np.array([l.vertex_index for l in me.loops],dtype='i4');before_uv=np.array([tuple(u.uv) for u in uv.data],dtype='f4');keep={p.index for p in me.polygons if p.normal.z>.02}
def jumps():
 edges={};rows=[]
 for p in me.polygons:
  if p.index not in keep:continue
  ids=list(p.loop_indices)
  for a,b in zip(ids,ids[1:]+ids[:1]):
   va,vb=me.loops[a].vertex_index,me.loops[b].vertex_index;k=tuple(sorted((va,vb)));edges.setdefault(k,[]).append((p.index,{va:np.array(uv.data[a].uv),vb:np.array(uv.data[b].uv)}))
 for k,values in edges.items():
  if len(values)!=2:continue
  dif=np.array([values[0][1][v]-values[1][1][v] for v in k]);wrap=dif-np.rint(dif);error=float(np.max(abs(wrap)))
  if error>1e-6:rows.append({'vertices':k,'faces':[v[0] for v in values],'wrapped_phase_jump':error,'raw_delta':dif.tolist()})
 return rows
before=jumps();assert before,'No source chart jump found; stop rather than invent correction'
new={}
for p in me.polygons:
 if p.index not in keep:continue
 for k in p.loop_indices:
  vid=me.loops[k].vertex_index
  if vid not in new:
   x,y,z=map(float,verts[vid]);cx=min(max(x,-ix),ix);cy=min(max(y,-iy),iy);dx,dy=x-cx,y-cy;rho=math.hypot(dx,dy)
   if rho<1e-9:U,V=x,y
   else:
    rise=z-zc;theta=math.atan2(rho,rise);radius=math.hypot(rho,rise);scale=radius*theta/rho;U,V=cx+dx*scale,cy+dy*scale
   new[vid]=(U/period+.5,V/period+.5)
  uv.data[k].uv=new[vid]
after=jumps();assert not after
new_uv=np.array([tuple(u.uv) for u in uv.data],dtype='f4');flat=[k for p in me.polygons if p.normal.z>.999 for k in p.loop_indices];assert np.array_equal(before_uv[flat],new_uv[flat]),'Flat top chart changed'
areas=[];cyl_ratios=[]
for p in me.polygons:
 if p.index not in keep:continue
 a=np.array([tuple(uv.data[k].uv) for k in p.loop_indices]);area=sum(np.cross(a[i]-a[0],a[i+1]-a[0]) for i in range(1,len(a)-1))/2;areas.append(float(area));assert area>0,'Degenerate or reversed visible chart'
 if p.normal.z<.999:
  # Real bevel-strip midlines avoid offset corner vertices. Compare the UV
  # interpolation actually used by the face to its geometry interpolation.
  for axis in [0,1]:
   mids=[];indices=list(p.loop_indices)
   for ka,kb in zip(indices,indices[1:]+indices[:1]):
    va,vb=me.loops[ka].vertex_index,me.loops[kb].vertex_index;pa,pb=verts[va].astype(float),verts[vb].astype(float)
    if abs(pa[axis]-pb[axis])>.08:
     mids.append(((pa+pb)/2,(np.array(uv.data[ka].uv)+np.array(uv.data[kb].uv))/2))
   if len(mids)==2:
    dist=np.linalg.norm(mids[0][0]-mids[1][0]);chart=np.linalg.norm((mids[0][1]-mids[1][1])*period);cyl_ratios.append(float(chart/dist))
assert cyl_ratios and min(cyl_ratios)>.98 and max(cyl_ratios)<1.02,cyl_ratios
assert np.array_equal(verts,np.array([tuple(v.co) for v in me.vertices],dtype='f4'));assert np.array_equal(topology,np.array([l.vertex_index for l in me.loops],dtype='i4'))
out=R/'scenes/07_porcelain_poisson_continuous_visible_atlas.blend';bpy.ops.wm.save_as_mainfile(filepath=str(out),compress=True);assert sha(src)==ref['source_sha256']
report={'source':str(src),'source_sha256':ref['source_sha256'],'scene':str(out),'scene_sha256':sha(out),'audited_outer_positive_Z_faces':len(keep),'discontinuous_shared_edges_before':len(before),'maximum_wrapped_phase_jump_before':max(x['wrapped_phase_jump'] for x in before),'phase_jumps_before':before,'discontinuous_shared_edges_after':0,'flat_top_UVs_bit_exact':True,'geometry_and_topology_bit_exact':True,'changed_loop_count':int(np.count_nonzero(np.any(before_uv!=new_uv,axis=1))),'minimum_positive_chart_face_area':min(areas),'straight_edge_chart_to_chord_ratio_range':[min(cyl_ratios),max(cyl_ratios)],'scope':'Outer visible top/rounded edge only. Inner and coating-termination charts retained. Corner radial development is continuous but not isometric; no curved-edge coating-flow solve claimed.','shader_coefficients_images_camera_lights_untouched':True,'visual_acceptance':False}
(R/'receipts/chart_audit.json').write_text(json.dumps(report,indent=2));ref.update(experiment_id='porcelain_continuous_atlas',source_scene=str(out.relative_to(WORK)),source_sha256=sha(out),purpose='Isolate a documented texture-chart phase discontinuity on the same Poisson deposition state and exact mesh',print_on_image='Porcelain · 28mm edge · continuous visible glaze chart');ref['pass_gates'].update(visible_chart_continuity=True,flat_top_chart_unchanged=True,geometry_unchanged=True,visual_acceptance=False);(R/'experiment.json').write_text(json.dumps(ref,indent=2));print('CONTINUOUS_ATLAS_READY',json.dumps({k:report[k] for k in ['scene_sha256','discontinuous_shared_edges_before','maximum_wrapped_phase_jump_before','straight_edge_chart_to_chord_ratio_range']}),flush=True)
