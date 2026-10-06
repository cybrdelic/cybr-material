"""Replace protruding r4 flange with a shared, continuous rounded seam profile."""
import bpy,math,sys,json,hashlib
from pathlib import Path
from mathutils import Vector
sys.path.insert(0,str(Path(__file__).resolve().parent))
from expansion.geometry import mesh
ROOT=Path(__file__).resolve().parents[1]

def outline(inset,z):
    hx=.06-inset;hy=.04-inset;r=.013-inset;v=[]
    for cx,cy,a in [(hx-r,hy-r,0),(-hx+r,hy-r,math.pi/2),(-hx+r,-hy+r,math.pi),(hx-r,-hy+r,3*math.pi/2)]:
        for i in range(17):
            q=a+i/16*math.pi/2;v.append((cx+r*math.cos(q),cy+r*math.sin(q),z))
    return v

def solid(name,rings,mat):
    verts=[];faces=[]
    for inset,z in rings:verts.extend(outline(inset,z))
    n=68
    for j in range(len(rings)-1):
        for i in range(n):a=j*n+i;b=j*n+(i+1)%n;faces.append((a,b,b+n,a+n))
    faces.append(tuple(reversed(range(n))));faces.append(tuple(range((len(rings)-1)*n,len(rings)*n)))
    ob=mesh(name,verts,faces,mat)
    for f in ob.data.polygons[:-2]:f.use_smooth=True
    return ob

for id in ('20_plastic_smooth','21_plastic_texture'):
    path=ROOT/'scenes'/f'{id}_r4_hero.blend';bpy.ops.wm.open_mainfile(filepath=str(path));s=bpy.context.scene
    m=bpy.data.objects['Plastic / curved molded shell'].data.materials[0];seam=bpy.data.objects['Plastic / 0.35mm mold parting line'].data.materials[0]
    for ob in list(s.objects):
        if ob.name.startswith('Plastic /'):bpy.data.objects.remove(ob,do_unlink=True)
    bottom=[(.003*(1-math.sin(t)),.003*(1-math.cos(t))) for t in [i/8*math.pi/2 for i in range(9)]]+[(0,.01165)]
    top=[(0,.012),(0,.024)]+[(.012*(1-math.cos(t)),.024+.012*math.sin(t)) for t in [i/16*math.pi/2 for i in range(1,17)]]
    solid('Plastic / lower half rounded profile',bottom,m);solid('Plastic / upper half rounded profile',top,m)
    connector=solid('Plastic / recessed seam connector',[(.0003,.01165),(.0003,.012)],seam);connector['seam_gap_m']=.00035;connector['inset_m']=.0003
    ground=bpy.data.objects['Studio / floor'];ground.location.z=-.0025
    fill=bpy.data.objects['Capture / fill'];w=.145;target=(0,0,.018);fill.location=(1.1*w,.5*w,.9*w+.018);fill.rotation_euler=(Vector(target)-fill.location).to_track_quat('-Z','Y').to_euler();fill.data.energy*=.5
    s['revision']='r5: closed molded halves share rounded outline; recessed 0.35mm seam, flush bottom, moved fill; materials unchanged'
    for view in ('hero','detail'):
        if view=='detail':
            s.camera.location=(w*.27,-w*.37,w*.45+.018);s.camera.data.ortho_scale=w*.38;target=(0,-w*.13,.018);s.camera.rotation_euler=(Vector(target)-s.camera.location).to_track_quat('-Z','Y').to_euler()
        out=ROOT/'scenes'/f'{id}_r5_{view}.blend';s.render.filepath=str(ROOT/'evidence'/f'{id}_r5_{view}.png');bpy.ops.wm.save_as_mainfile(filepath=str(out));print('READY_SCENE',out,flush=True)
