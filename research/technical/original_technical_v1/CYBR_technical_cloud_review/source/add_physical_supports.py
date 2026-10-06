"""Support the authored chain drape at its outer rings and ground PETG specimens."""
import bpy,math,sys
from pathlib import Path
from mathutils import Vector
sys.path.insert(0,str(Path(__file__).resolve().parent))
from expansion.geometry import box,plain
ROOT=Path(__file__).resolve().parents[1]
for id in ('14_chainmail','22_petg','23_petg_transparent'):
    bpy.ops.wm.open_mainfile(filepath=str(ROOT/'scenes'/f'{id}_r4_hero.blend'));s=bpy.context.scene;s.objects['Studio / floor'].location.z=-.0025
    if id=='14_chainmail':
        m=plain('Fixture / matte dark steel',(.023,.029,.035),.5,metal=.7);R=.14;wire=.00045;rod=.0018
        for sign in (-1,1):
            angle=sign*.04/R;ex=Vector((math.cos(angle),0,math.sin(angle)));ez=Vector((-math.sin(angle),0,math.cos(angle)))
            ring_center=Vector((R*math.sin(angle),0,.010+R*(1-math.cos(angle))))
            edge=ring_center+ex*(sign*.004);center=edge-ez*(wire+rod)
            bpy.ops.mesh.primitive_cylinder_add(vertices=48,radius=rod,depth=.095,location=center,rotation=(math.pi/2,0,0));o=bpy.context.object;o.name='Fixture / support under outer ring extremities';o.data.materials.append(m)
            for p in o.data.polygons:p.use_smooth=True
            for y in (-.043,.043):box('Fixture / grounded rail support',(center.x,y,center.z/2),(.0034,.0034,center.z),m,.0003)
        s['support_note']='Rigid authored links rest on two rails at outermost ring extremities; links are unchanged. Static drape, not cloth simulation.'
    else:
        s.objects['PETG / unified printed wall rim and solid floor'].location.z=-.0012
        for o in s.objects:
            if o.name.startswith('Optics / external background'):o.location.z-=.002
        s['support_note']='Closed specimen outer floor sits on studio surface; materials and wall geometry unchanged.'
    s['revision']='r5: physical supports/contact, unchanged specimen materials'
    out=ROOT/'scenes'/f'{id}_r5_hero.blend';s.render.filepath=str(ROOT/'evidence'/f'{id}_r5_hero.png');bpy.ops.wm.save_as_mainfile(filepath=str(out));print('READY_SCENE',out,flush=True)
