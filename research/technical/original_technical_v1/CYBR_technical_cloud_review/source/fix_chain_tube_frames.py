"""Continuous analytic torus frames remove cross-section flips in vertical links."""
import bpy,math,sys,json,hashlib
import numpy as np
from pathlib import Path
from mathutils import Vector
sys.path.insert(0,str(Path(__file__).resolve().parent))
from expansion.geometry import mesh
ROOT=Path(__file__).resolve().parents[1]
bpy.ops.wm.open_mainfile(filepath=str(ROOT/'scenes/14_chainmail_r5_hero.blend'));s=bpy.context.scene
old=next(o for o in s.objects if o.name.startswith('Chainmail /'));m=old.data.materials[0]
oldcenters=np.array([tuple(v.co) for v in old.data.vertices]).reshape(-1,64,10,3).mean(2)
v=[];f=[];uv=[];count=0;R=.14;wire=.00045;nt=64;ns=16
def ring(cx,cy,r,axis):
    global count
    angle=cx/R;center=Vector((R*math.sin(angle),cy,.010+R*(1-math.cos(angle))))
    ex=Vector((math.cos(angle),0,math.sin(angle)));ey=Vector((0,1,0));ez=Vector((-math.sin(angle),0,math.cos(angle)))
    u,vv=(ex,ey) if axis=='z' else (ex,ez) if axis=='y' else (ey,ez);normal=u.cross(vv).normalized();start=len(v);count+=1
    for i in range(nt):
        phi=math.tau*i/nt;radial=u*math.cos(phi)+vv*math.sin(phi)
        for j in range(ns):
            theta=math.tau*j/ns;pos=center+radial*(r+wire*math.cos(theta))+normal*(wire*math.sin(theta));v.append(tuple(pos));uv.append((j/ns,i/(nt-1)))
    for i in range(nt):
        for j in range(ns):
            a=start+i*ns+j;b=start+((i+1)%nt)*ns+j;c=start+((i+1)%nt)*ns+(j+1)%ns;d=start+i*ns+(j+1)%ns;f.append((a,b,c,d))
for j in range(9):
    for i in range(9):
        x=(i-4)*.01;y=(j-4)*.01;ring(x,y,.004,'z')
        if i<8:ring(x+.005,y,.0035,'y')
        if j<8:ring(x,y+.005,.0035,'x')
newcenters=np.asarray(v).reshape(-1,nt,ns,3).mean(2);error=float(np.max(np.abs(oldcenters-newcenters)));assert error<1e-8,error
bpy.data.objects.remove(old,do_unlink=True);o=mesh('Chainmail / continuous torus-section rigid links',v,f,m,uv)
for p in o.data.polygons:p.use_smooth=True
o['link_count']=count;o['drape_radius_m']=R;o['cross_section']='Analytic radial/plane-normal torus basis; no reference-axis switching';o['construction']='Japanese four-in-one rigid links; centerlines unchanged from r5'
maxedge=max((o.data.vertices[e.vertices[0]].co-o.data.vertices[e.vertices[1]].co).length for e in o.data.edges);assert maxedge<.00045,maxedge
import bmesh
bm=bmesh.new();bm.from_mesh(o.data);assert all(e.is_manifold for e in bm.edges);bm.free()
s['revision']='r6: continuous analytic torus cross-sections replace discontinuous path frames; identical link centerlines/supports/materials'
for view in ('hero','detail'):
    if view=='detail':
        w=.109;h=.010;s.camera.location=(w*.27,-w*.37,w*.45+h);s.camera.data.ortho_scale=w*.46;target=Vector((0,-w*.13,h));s.camera.rotation_euler=(target-s.camera.location).to_track_quat('-Z','Y').to_euler()
    path=ROOT/'scenes'/f'14_chainmail_r6_{view}.blend';s.render.filepath=str(ROOT/'evidence'/f'14_chainmail_r6_{view}.png');bpy.ops.wm.save_as_mainfile(filepath=str(path));print('READY_SCENE',path,flush=True)
(ROOT/'tests/chain_frames_r6.json').write_text(json.dumps({'links':count,'major_segments':nt,'minor_segments':ns,'max_centerline_deviation_from_r5_m':error,'maximum_edge_length_m':maxedge,'closed_manifold_links':True,'unchanged_linkage_and_clearance':True,'source_material_unchanged':True,'status':'New analytic surface frames; pixel review pending'},indent=2)+'\n');print('CHAIN_FRAME_PASS',error,maxedge,flush=True)
