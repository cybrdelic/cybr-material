"""Portable clean Cycles views. Build defaults to CPU-only; render is explicit."""
import argparse,hashlib,json,math,sys,time
from pathlib import Path
import bpy
from mathutils import Vector
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from expansion.geometry import build,box,plain
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'expansion'

def aim(o,target):o.rotation_euler=(Vector(target)-o.location).to_track_quat('-Z','Y').to_euler()
def area(name,location,target,size,power,color=(1,1,1)):
    data=bpy.data.lights.new(name,'AREA');data.energy=power;data.shape='DISK';data.size=size;data.color=color
    o=bpy.data.objects.new(name,data);bpy.context.collection.objects.link(o);o.location=location;aim(o,target);return o

def configure(s,objects,width,view,lighting,size,samples,device):
    scene=bpy.context.scene;scene.render.engine='CYCLES';scene.cycles.samples=samples
    scene.cycles.adaptive_min_samples=128;scene.cycles.adaptive_threshold=.008
    scene.cycles.use_denoising=True;scene.cycles.denoiser='OPENIMAGEDENOISE';scene.cycles.denoising_input_passes='RGB_ALBEDO_NORMAL'
    scene.cycles.max_bounces=12;scene.cycles.transmission_bounces=12;scene.cycles.transparent_max_bounces=12
    scene.cycles.seed=s['seed']+{'hero':0,'grazing':1,'detail':2}[view]
    scene.render.resolution_x=size;scene.render.resolution_y=size;scene.render.resolution_percentage=100
    scene.render.image_settings.file_format='PNG';scene.render.image_settings.color_mode='RGB';scene.render.image_settings.color_depth='16'
    scene.view_settings.view_transform='AgX';scene.view_settings.look='AgX - Medium High Contrast';scene.view_settings.exposure=0
    if not scene.world:scene.world=bpy.data.worlds.new('Capture / studio world')
    scene.world.use_nodes=True;scene.world.node_tree.nodes['Background'].inputs[0].default_value=(.18,.21,.25,1);scene.world.node_tree.nodes['Background'].inputs[1].default_value=.22
    if device!='CPU':
        prefs=bpy.context.preferences.addons['cycles'].preferences;prefs.compute_device_type=device;prefs.get_devices()
        for d in prefs.devices:d.use=d.type==device
        scene.cycles.device='GPU'
    else:scene.cycles.device='CPU';scene.render.threads_mode='FIXED';scene.render.threads=4
    scene.unit_settings.system='METRIC';scene.unit_settings.scale_length=1
    ground=plain('Studio / warm neutral',(.15,.145,.13),.82)
    box('Studio / physical ground',(0,0,-.003),(width*30,width*30,.004),ground)
    # Physical reference stripe behind PETG makes transmission legible.
    if s.get('transmission',0):
        marker=plain('Studio / transmission reference',( .03,.045,.065),.5)
        box('Studio / solid transmission marker',(-.002,width*.36,.011),(.004,.001,.022),marker,.0001)
    target=(0,0,.012 if s['family'] in ('PETG','Grass') else .005)
    bpy.ops.object.camera_add();camera=bpy.context.object;camera.name='Capture / '+view;scene.camera=camera;camera.data.type='ORTHO';camera.data.clip_start=.0001;camera.data.clip_end=100
    if view=='hero':camera.location=(width*.86,-width*1.13,width*1.0);camera.data.ortho_scale=width*1.42
    elif view=='grazing':camera.location=(width*.8,-width*1.1,width*.43);camera.data.ortho_scale=width*1.1
    else:camera.location=(width*.26,-width*.40,width*.40);camera.data.ortho_scale=width*.45;target=(0,-width*.12,target[2])
    aim(camera,target)
    # Power scales with specimen area: identical illuminance across study scales.
    gain=width*width
    area('Capture / soft key',(-width*.7,-width*.65,width*1.4),target,width*1.05,36*gain,(1,.92,.82))
    area('Capture / grazing rim',(width*.55,width*.25,width*.60),target,width*.75,22*gain,(.80,.90,1))
    area('Capture / fill',(width*.8,-width*.55,width*.55),target,width*1.2,5*gain)
    if lighting=='hard':
        lamp=bpy.data.objects['Capture / soft key'];lamp.data.size=width*.09
    elif lighting=='grazing':
        lamp=bpy.data.objects['Capture / soft key'];lamp.location=(-width,-width*.2,width*.2);aim(lamp,target)
    scene['recipe_id']=s['id'];scene['physical_scale_m']=width;scene['capture_view']=view;scene['capture_lighting']=lighting
    return scene

def receipt(scene,s,objects,path,status,elapsed=0):
    meshes=[o for o in objects if o.type=='MESH']
    data=dict(recipe=s['id'],name=s['name'],status=status,seed=s['seed'],cycles_seed=scene.cycles.seed,
              view=scene['capture_view'],lighting=scene['capture_lighting'],size=scene.render.resolution_x,
              samples=scene.cycles.samples,min_samples=scene.cycles.adaptive_min_samples,threshold=scene.cycles.adaptive_threshold,
              denoiser='OpenImageDenoise with albedo and normal guides',device=scene.cycles.device,
              blender_version=bpy.app.version_string,elapsed_seconds=elapsed,
              source_baseline='de4146ac759164e29323b152d37f64417c00c3fb',
              camera=dict(location=list(scene.camera.location),rotation=list(scene.camera.rotation_euler),ortho_scale=scene.camera.data.ortho_scale),
              lights=[dict(name=o.name,location=list(o.location),power_W=o.data.energy,size_m=o.data.size) for o in scene.objects if o.type=='LIGHT'],
              topology=dict(objects=len(meshes),vertices=sum(len(o.data.vertices) for o in meshes),faces=sum(len(o.data.polygons) for o in meshes)),
              settings=s,source_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(__file__).parent.glob('*.py')})
    if path.exists():data.update(file=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),bytes=path.stat().st_size)
    return data

def main():
    p=argparse.ArgumentParser();p.add_argument('--only',nargs='+',required=True);p.add_argument('--view',choices=['hero','grazing','detail'],default='hero')
    p.add_argument('--lighting',choices=['soft','grazing','hard'],default='soft');p.add_argument('--size',type=int,default=1024);p.add_argument('--samples',type=int,default=512)
    p.add_argument('--device',choices=['CPU','OPTIX','CUDA'],default='CPU');p.add_argument('--render',action='store_true');p.add_argument('--parameters',type=Path)
    a=p.parse_args(sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else [])
    if a.render and (a.size<768 or a.samples<384):raise ValueError('Clean review captures need >=768 pixels, >=384 samples and 128 minimum')
    for id in a.only:
        bpy.ops.wm.read_factory_settings(use_empty=True)
        s,objects,width=build(id,a.parameters);scene=configure(s,objects,width,a.view,a.lighting,a.size,a.samples,a.device)
        stem=f'{id}_{a.view}_{a.lighting}';scenes=OUT/'scenes';renders=OUT/'renders';receipts=OUT/'receipts'
        for folder in (scenes,renders,receipts):folder.mkdir(parents=True,exist_ok=True)
        blend=scenes/(stem+'.blend');scene.render.filepath=str(renders/(stem+'.png'))
        for image in bpy.data.images:
            if image.source=='FILE':image.pack()
        bpy.ops.wm.save_as_mainfile(filepath=str(blend));data=receipt(scene,s,objects,blend,'editable scene, not rendered')
        if a.render:
            start=time.time();bpy.ops.render.render(write_still=True);elapsed=time.time()-start
            data=receipt(scene,s,objects,Path(scene.render.filepath),'rendered, awaiting visual review',elapsed)
        (receipts/(stem+'.json')).write_text(json.dumps(data,indent=2)+'\n')
        print('EXPANSION_READY',stem,data['status'],data['topology'],flush=True)
if __name__=='__main__':main()
