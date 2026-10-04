"""Build displaced inspection patches reused by the Cycles scene renderer.

The standalone legacy preview CLI uses Eevee and labels its outputs accordingly.
Geometry uses Height_Macro; the shader uses residual Normal_Micro_OpenGL.
These two bands reconstruct the authored relief without doubling it.
"""
import argparse,importlib.util,json,sys
from pathlib import Path
import bpy
argv=sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else []
p=argparse.ArgumentParser();p.add_argument('--draft',action='store_true');p.add_argument('--only',nargs='*',default=[])
p.add_argument('--size',type=int,default=1536);p.add_argument('--samples',type=int,default=96);p.add_argument('--save',action='store_true')
p.add_argument('--reuse-previews',action='store_true',help='Build scenes but preserve already reviewed previews')
args=p.parse_args(argv)
ROOT=Path(__file__).resolve().parents[1]
sys.argv=['blender','--']+(['--draft'] if args.draft else [])
module=importlib.util.spec_from_file_location('atelier',ROOT/'source/build_blender.py');atelier=importlib.util.module_from_spec(module);module.loader.exec_module(atelier)
atelier.args.samples=args.samples
CENTERS={
 '01_calacatta_oro':(.29,.32),'02_roman_travertine':(.51,.49),
 '03_american_walnut':(.50,.50),'04_fumed_oak':(.50,.50),
 '05_champagne_brass':(.46,.43),'06_blackened_steel':(.46,.52),
 '07_bone_porcelain':(.49,.44),'08_saddle_leather':(.46,.46),
 '09_lime_plaster':(.56,.50),'10_natural_linen':(.52,.49)}

def macro(spec):
 span=spec['preview_diameter_m'];side=span*1.8
 scene=atelier.studio_scene('CYBR MATERIAL 3 / Macro / '+spec['name'],args.size,args.size,world_strength=.13)
 scene.world.node_tree.nodes['Background'].inputs['Color'].default_value=(.86,.89,.94,1)
 scene.view_settings.look='AgX - Medium High Contrast';scene.view_settings.exposure=-.1
 mat=atelier.material(spec)
 mat.asset_clear()
 shader=next(node for node in mat.node_tree.nodes if node.type=='GROUP')
 shader.inputs['Displacement Mode'].default_value=True
 bpy.ops.mesh.primitive_grid_add(x_subdivisions=385,y_subdivisions=385,size=side)
 ob=bpy.context.object;ob.name=spec['name']+' / physically displaced inspection patch'
 ob.data.materials.append(mat)
 for poly in ob.data.polygons:poly.use_smooth=True
 uv=ob.data.uv_layers.active;center=CENTERS[spec['id']]
 for loop in uv.data:loop.uv=(center[0]+(loop.uv.x-.5)*side/spec['tile_m'],center[1]+(loop.uv.y-.5)*side/spec.get('tile_y_m',spec['tile_m']))
 texture=bpy.data.textures.new('Macro relief / '+spec['name'],type='IMAGE')
 texture.image=bpy.data.images.load(str(atelier.MAP_ROOT/spec['id']/'Height_Macro.png'),check_existing=True)
 texture.image.colorspace_settings.name='Non-Color';texture.use_interpolation=True;texture.extension='EXTEND' if spec.get('tiling')=='finite-cut' else 'REPEAT'
 mod=ob.modifiers.new('Measured macro relief; shader adds residual micro relief','DISPLACE')
 mod.texture=texture;mod.texture_coords='UV';mod.uv_layer=uv.name;mod.direction='Z'
 mod.mid_level=.5;mod.strength=spec['height_scale_m']
 atelier.camera(scene,(span*.12,-span*1.20,span*1.85),(0,0,0),span*1.06)
 scene.camera.data.clip_start=span*.005;scene.camera.data.clip_end=span*100
 atelier.area(scene,'Grazing window / long silk',(-1.5*span,-.7*span,.56*span),(0,0,0),40*span*span,.52*span,1.15*span,color=(1,.96,.89))
 atelier.area(scene,'Soft daylight fill',(span*.8,span*.2,span*1.7),(0,0,0),12*span*span,1.4*span,color=(.86,.92,1))
 atelier.area(scene,'Reflection strip',(span*.25,span*1.3,span*1.05),(0,0,0),9*span*span,.22*span,1.5*span,color=(1,.985,.95))
 scene['inspection_width_m']=span;scene['physical_tile_m']=spec['tile_m'];scene['relief_scale_m']=spec['height_scale_m']
 scene['geometry_height']='Height_Macro.png';scene['normal']='Normal_Micro_OpenGL.png'
 return scene

def main():
 bpy.ops.wm.read_factory_settings(use_empty=True)
 scenes=[];metadata={}
 for spec in atelier.SPECS:
  if args.only and spec['id'] not in args.only:continue
  if not (atelier.MAP_ROOT/spec['id']/'Height_Macro.png').exists():continue
  print('MACRO_BEGIN '+spec['id'],flush=True)
  scene=macro(spec);scenes.append(scene)
  suffix='-draft' if args.draft else ''
  dest=ROOT/'previews'/(spec['id']+suffix+'.png')
  scene.render.filepath=str(dest)
  if not (args.reuse_previews and dest.exists()):bpy.ops.render.render(write_still=True)
  metadata[spec['id']]=dict(renderer='Blender 4.3 Eevee',size=args.size,samples=args.samples,inspection_width_m=spec['preview_diameter_m'],postprocessing='Native AgX display transform and antialiasing; no denoising or blur',geometry_height='Height_Macro.png',residual_normal='Normal_Micro_OpenGL.png',uv_center=CENTERS[spec['id']])
  print('MACRO_READY '+str(dest),flush=True)
  # Release cached image pixels between scenes; datablocks retain their
  # relative file references and reload when an inspection scene is opened.
  for im in bpy.data.images:
   if im.source=='FILE':im.buffers_free()
 if args.save and scenes:
  for im in bpy.data.images:
   if im.source=='FILE':im.filepath='//../materials/'+Path(im.filepath).parent.name+'/'+Path(im.filepath).name
  for tex in bpy.data.textures:
   if tex.image and tex.image.source=='FILE':tex.image.filepath='//../materials/'+Path(tex.image.filepath).parent.name+'/'+Path(tex.image.filepath).name
  bpy.context.window.scene=scenes[2] if len(scenes)>2 else scenes[0]
  text=bpy.data.texts.new('READ ME / Macro Inspection')
  text.write('CYBR MATERIAL 3 / Natural Detail + Wear\nSwitch scenes to inspect all ten surfaces.\n'
             'Geometry uses Height_Macro via a UV displacement modifier at the measured meter scale.\n'
             'Displacement Mode selects the residual micro normal; it is enabled on these samples.\n'
             'Do not combine the full normal with this displacement: that would double the relief.\n'
             'Each camera shows a measured close-up width; scene custom properties record it.\n'
             'Preview rendering: Eevee, native antialiasing, no smoothing or denoising.\n'
             'These are custom procedural surfaces, not photogrammetry or scans.\n')
  for screen in bpy.data.screens:
   for area in screen.areas:
    if area.type=='VIEW_3D':
     space=area.spaces.active;space.shading.type='MATERIAL'
     space.shading.use_scene_lights=True;space.shading.use_scene_world=True
     space.region_3d.view_perspective='CAMERA';space.clip_start=.00001
  bpy.ops.wm.save_as_mainfile(filepath=str(ROOT/'blender/CYBR_Macro_Inspection.blend'),compress=True)
 name='macro_rendering-draft.json' if args.draft else 'macro_rendering.json'
 (ROOT/'previews'/name).write_text(json.dumps(metadata,indent=2)+'\n')
 print('MACROS_COMPLETE',flush=True)

if __name__=='__main__':main()
