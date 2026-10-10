# Explicit external inputs replace machine-specific historical file locations.
from pathlib import Path as _ResearchInputPath
import sys as _research_input_sys
for _research_parent in _ResearchInputPath(__file__).resolve().parents:
    if (_research_parent / "research_inputs.py").is_file():
        _research_input_sys.path.insert(0, str(_research_parent))
        break
from research_inputs import required_input as _required_research_input
import bpy, json, hashlib, importlib.util, resource
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STUDIO=Path(str(_required_research_input("studio-source")))
BASELINE=Path(str(_required_research_input("baseline-scene")))
assert hashlib.sha256(BASELINE.read_bytes()).hexdigest()=='55df658fd44a0fb17bab14b46692c1c252b2d2658a967910e8c281f82e743879'
bpy.ops.wm.open_mainfile(filepath=str(BASELINE))
spec=importlib.util.spec_from_file_location('frozen_studio',STUDIO);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
receipt=mod.configure('24_carpet')
bpy.context.view_layer.update()
scene=bpy.context.scene
def attrs(ob,names):
    return {n:(list(getattr(ob,n)) if hasattr(getattr(ob,n),'__len__') and not isinstance(getattr(ob,n),str) else getattr(ob,n)) for n in names}
report={'label':'Reconstructed baseline rig from recovered studio helper, not historical capture-byte identity','receipt':receipt,'source_sha256':hashlib.sha256(BASELINE.read_bytes()).hexdigest(),'studio_sha256':hashlib.sha256(STUDIO.read_bytes()).hexdigest(),'objects':{}}
for ob in scene.objects:
    if not ob.name.startswith('Slab /'):continue
    d={'type':ob.type,'matrix_world':list(map(list,ob.matrix_world))}
    if ob.type=='CAMERA':d['camera']=attrs(ob.data,['type','ortho_scale','clip_start','clip_end','lens','shift_x','shift_y'])
    if ob.type=='LIGHT':d['light']=attrs(ob.data,['type','energy','shape','size','color','use_shadow'])
    if ob.type=='LIGHT':
        assert (ob.matrix_world.translation-ob.location).length<1e-10
        assert ob.location.length>.001
    if ob.type=='MESH':
        d['vertices']=[list(v.co) for v in ob.data.vertices];d['faces']=[list(p.vertices) for p in ob.data.polygons]
        p=ob.data.materials[0].node_tree.nodes.get('Principled BSDF');d['material']={'name':ob.data.materials[0].name,'Base Color':list(p.inputs['Base Color'].default_value),'Roughness':p.inputs['Roughness'].default_value}
    report['objects'][ob.name]=d
bg=scene.world.node_tree.nodes['Background'];report['world']={'name':scene.world.name,'color':list(bg.inputs[0].default_value),'strength':bg.inputs[1].default_value}
report['view']=attrs(scene.view_settings,['view_transform','look','exposure','gamma'])
report['render']=attrs(scene.render,['engine','resolution_x','resolution_y','resolution_percentage','film_transparent'])
report['cycles']=attrs(scene.cycles,['device','use_denoising','use_light_tree','max_bounces','transmission_bounces'])
report['unit']=attrs(scene.unit_settings,['system','scale_length'])
report['image']=attrs(scene.render.image_settings,['file_format','color_mode','color_depth'])
report['rss_mib']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
out=ROOT/'receipts/reconstructed_baseline_rig.json';out.write_text(json.dumps(report,indent=2));print(str(out),report['rss_mib'],flush=True)
