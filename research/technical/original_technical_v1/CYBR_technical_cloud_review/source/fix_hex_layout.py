"""Match pointy-top hex orientation to the existing staggered triangular center lattice."""
import bpy,math,json
from pathlib import Path
from mathutils import Vector
ROOT=Path(__file__).resolve().parents[1];src=ROOT/'scenes/15_future_ceramic_r4_hero.blend';bpy.ops.wm.open_mainfile(filepath=str(src));s=bpy.context.scene
plates=[o for o in s.objects if o.name.startswith('Engineered / discrete ceramic')];m=plates[0].data.materials[0]
for o in plates:bpy.data.objects.remove(o,do_unlink=True)
p=.016;r=p/1.8;apothem=r*math.sqrt(3)/2;plates=[]
for j in range(-3,4):
    for i in range(-3,4):
        x=i*p+(j%2)*p/2;y=j*p*math.sqrt(3)/2
        if abs(x)+apothem>.056:continue
        # Blender's unrotated six-sided cylinder has vertices at 30,90,... degrees.
        bpy.ops.mesh.primitive_cylinder_add(vertices=6,radius=r,depth=.003,location=(x,y,.006));o=bpy.context.object;o.name='Engineered / separated pointy-top ceramic hex';o.data.materials.append(m)
        mod=o.modifiers.new('Rounded sintered ceramic edges','BEVEL');mod.width=.00035;mod.segments=4;o.modifiers.new('Weighted corner normals','WEIGHTED_NORMAL');plates.append(o)
bpy.context.view_layer.update()
def polygon(o):
    vv=sorted(set((round((o.matrix_world@v.co).x,9),round((o.matrix_world@v.co).y,9)) for v in o.data.vertices))
    return sorted(vv,key=lambda q:math.atan2(q[1]-o.location.y,q[0]-o.location.x))
def gap(a,b):
    largest=-1
    for poly in (a,b):
        for p,q in zip(poly,poly[1:]+poly[:1]):
            nx=-(q[1]-p[1]);ny=q[0]-p[0];norm=math.hypot(nx,ny);nx/=norm;ny/=norm
            aa=[x*nx+y*ny for x,y in a];bb=[x*nx+y*ny for x,y in b]
            largest=max(largest,min(bb)-max(aa),min(aa)-max(bb))
    return largest
mingap=1
for i,a in enumerate(plates):
    for b in plates[:i]:
        g=gap(polygon(a),polygon(b));assert g>0,(a.name,b.name,g)
        mingap=min(mingap,g)
assert len(plates)==45,len(plates)
s['revision']='r5: pointy-top orientation consistent with triangular center lattice; positive plate gaps, recovered edge plates, material unchanged'
out=ROOT/'scenes/15_future_ceramic_r5_hero.blend';s.render.filepath=str(ROOT/'evidence/15_future_ceramic_r5_hero.png');bpy.ops.wm.save_as_mainfile(filepath=str(out))
(ROOT/'tests/hex_layout_r5.json').write_text(json.dumps({'plate_count':len(plates),'minimum_unbeveled_separation_m':mingap,'all_pairs_strictly_separated':True,'method':'Separating-axis test of world-space convex hex footprints; bevels recede into these footprints','source_material_unchanged':True},indent=2)+'\n')
print('HEX_LAYOUT_PASS',len(plates),mingap,flush=True);print('READY_SCENE',out,flush=True)
