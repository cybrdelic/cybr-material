"""Scene evidence checks, no renderer invoked. Geometry/shaders do not certify image quality."""
import bpy,sys,json,hashlib,collections,math
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def edge_counts(o):
    counts=collections.Counter()
    for p in o.data.polygons:
        ids=list(p.vertices)
        for a,b in zip(ids,ids[1:]+ids[:1]):counts[tuple(sorted((a,b)))]+=1
    return sum(v==1 for v in counts.values()),sum(v>2 for v in counts.values())
report={}
for path in sorted((ROOT/'scenes').glob('*r4_hero.blend')):
    bpy.ops.wm.open_mainfile(filepath=str(path));s=bpy.context.scene;id=s['recipe_id'];record={'scene':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'status':'Geometry/material checks only; visual review pending'}
    record['meshes']=sum(o.type=='MESH' for o in s.objects)
    missing=[im.filepath for im in bpy.data.images if im.source=='FILE' and not im.packed_file and not Path(bpy.path.abspath(im.filepath)).is_file()];assert not missing,missing
    record['all_dependencies_available']=True
    if id in ('18_lcd','19_crt'):
        cover=next(o for o in s.objects if o.name.startswith('Display / planar coverglass') or o.name.startswith('Display / curved faceplate'))
        b,n=edge_counts(cover);assert (b,n)==(0,0)
        p=cover.data.materials[0].node_tree.nodes['Principled BSDF'];assert abs(p.inputs['IOR'].default_value-1.52)<1e-5
        assert p.inputs['Transmission Weight'].default_value==1 and p.inputs['Alpha'].default_value==1
        pixels=bpy.data.objects['Display / independently colored real RGB islands'];nx,ny=pixels['pixel_grid'];assert len(pixels.data.polygons)==nx*ny*3
        assert not any(n.type=='TEX_IMAGE' for m in pixels.data.materials for n in m.node_tree.nodes)
        record.update(glass_thickness_m=cover['glass_thickness_m'],glass_boundary_edges=b,glass_nonmanifold_edges=n,subpixel_count=nx*ny*3,pixel_pitch_m=pixels['pixel_pitch_m'],bow_m=cover['center_bow_m'],emission_mask_resampling=False)
    if id in ('22_petg','23_petg_transparent'):
        o=bpy.data.objects['PETG / unified printed wall rim and solid floor'];b,n=edge_counts(o);assert (b,n)==(0,0)
        p=o.data.materials[0].node_tree.nodes['Principled BSDF'];assert abs(p.inputs['IOR'].default_value-1.57)<1e-5
        if id=='23_petg_transparent':assert p.inputs['Transmission Weight'].default_value>.95
        record.update(boundary_edges=b,nonmanifold_edges=n,wall_thickness_m=o['wall_thickness_m'],base_thickness_m=o['base_thickness_m'],layer_pitch_m=o['layer_height_m'],seam_peak_m=o['seam_peak_m'],ior=p.inputs['IOR'].default_value,alpha=p.inputs['Alpha'].default_value)
    if id=='14_chainmail':
        o=next(o for o in s.objects if o.name.startswith('Chainmail /'));b,n=edge_counts(o);assert (b,n)==(0,0)
        arr=np.array([tuple(v.co) for v in o.data.vertices]).reshape(-1,64,10,3).mean(2)
        # Exact all-point separation between all distinct rigid-ring centerlines.
        mind=1;pair=None
        for i in range(len(arr)):
            for j in range(i):
                if np.linalg.norm(arr[i].mean(0)-arr[j].mean(0))>.010:continue
                d=float(np.sqrt(((arr[i,:,None,:]-arr[j,None,:,:])**2).sum(-1)).min())
                if d<mind:mind=d;pair=(j,i)
        clearance=mind-.0009;assert clearance>0,(clearance,pair)
        def link(a,b):
            da=np.roll(a,-1,0)-a;db=np.roll(b,-1,0)-b;ma=a+da/2;mb=b+db/2;r=ma[:,None]-mb[None,:]
            return float((np.einsum('ijk,ijk->ij',r,np.cross(da[:,None,:],db[None,:,:]))/(np.linalg.norm(r,axis=-1)**3)).sum()/(4*math.pi))
        edges=[]
        for i in range(len(arr)):
            for j in range(i):
                if np.linalg.norm(arr[i].mean(0)-arr[j].mean(0))>.010:continue
                lk=link(arr[i],arr[j])
                if abs(lk)>.85:edges.append((i,j))
        connected={0}
        while True:
            nxt=connected|{b for a,b in edges if a in connected}|{a for a,b in edges if b in connected}
            if len(nxt)==len(connected):break
            connected=nxt
        assert len(connected)==len(arr),(len(connected),len(arr))
        record.update(links=len(arr),linked_pairs=len(edges),one_connected_network=True,centerline_sampled_wire_clearance_m=clearance,boundary_edges=b,nonmanifold_edges=n)
    report[id]=record;print('TECHNICAL_CHECK_PASS',id,flush=True)
(ROOT/'tests'/'technical_r4_checks.json').write_text(json.dumps(report,indent=2)+'\n')
