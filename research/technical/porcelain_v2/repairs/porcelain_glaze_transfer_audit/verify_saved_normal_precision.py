"""Double-precision read-only normal check; avoids mathutils.angle float rounding."""
import bpy, hashlib, json, math, sys
from pathlib import Path
from mathutils import Vector
out=Path(sys.argv[sys.argv.index('--')+1])
r=json.loads((out/'build_receipt.json').read_text())
rows=[]
for row in r['scenes']:
    if row['treatment']!='candidate':continue
    path=Path(row['scene']);digest=hashlib.sha256(path.read_bytes()).hexdigest();assert digest==row['sha256']
    bpy.ops.wm.open_mainfile(filepath=str(path),load_ui=False)
    me=bpy.context.scene.objects['07 new / finite60micron mean glaze shell'].data
    half=len(me.vertices)//2;errors=[]
    for p in me.polygons:
        ids=[me.loops[k].vertex_index for k in p.loop_indices]
        for k in p.loop_indices:
            vid=me.loops[k].vertex_index
            if p.normal.z>.02 or all(v>=half for v in ids):
                co=me.vertices[vid+half if p.normal.z>.02 else vid].co
                v=Vector((co.x-min(max(co.x,-.0474),.0474),co.y-min(max(co.y,-.0474),.0474),co.z-.0012))
                for i in range(3):
                    if abs(v[i])<1e-8:v[i]=0
                v.normalize()
                if p.normal.z<=.02:v=-v
            else:v=p.normal.copy()
            a=tuple(me.corner_normals[k].vector);b=tuple(v)
            cross=(a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
            errors.append(math.degrees(math.atan2(math.sqrt(sum(x*x for x in cross)),sum(x*y for x,y in zip(a,b)))))
    assert max(errors)<.02,max(errors)
    assert hashlib.sha256(path.read_bytes()).hexdigest()==digest
    row['normal_verification']['max_saved_custom_normal_error_degrees']=max(errors)
    row['normal_verification']['angular_measurement']='double-precision atan2 of cross norm and dot'
    rows.append({'scene':str(path),'sha256':digest,'max_saved_normal_error_degrees':max(errors)})
(out/'normal_precision.json').write_text(json.dumps(rows,indent=2))
(out/'build_receipt.json').write_text(json.dumps(r,indent=2))
print('PRECISE_NORMAL_CHECK',json.dumps(rows),flush=True)
