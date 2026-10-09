"""One supported selected-material interface: python -m materials --help."""
import argparse,datetime,importlib.util,json,os,shutil,subprocess,sys,time
from pathlib import Path
from .core import *

def run(command,log,accepted=(0,)):
 with Path(log).open('x') as f:
  p=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'})
 if p.returncode not in accepted:raise ContractError('Command failed; inspect '+str(log))
 return p.returncode

def destination(path):
 path=Path(path).resolve()
 if path.exists():raise ContractError('Output directory exists; use a new capture/build ID: '+str(path))
 if path.is_relative_to(workspace()) and any(path.is_relative_to(resolve(m['selected']['scene']).parent) for m in registry()['materials'] if m['selected']):
  raise ContractError('Cannot write inside selected source directories')
 path.mkdir(parents=True);return path

def build(identifier,out):
 m=material(identifier);source=validate_selection(m);out=destination(out)
 run([os.environ.get('BLENDER_BIN','blender'),'-b','--threads','1','--python-exit-code','1','--python',str(ROOT/'materials/inspect_scene.py'),'--',source['scene'],source['scene_sha256'],str(out/'inspection.json')],out/'inspect.log')
 inspection=read(out/'inspection.json')
 validate_inspection(m,inspection)
 b={'schema_version':1,'id':m['id'],'registry_sha256':sha(ROOT/'registry.json'),'source':source,'inspection':inspection,'contract':m.get('map_contract'),'quality':m['quality'],'build_kind':'validated immutable selected scene and all actual shader dependencies; no geometry or image duplication'}
 write_new(out/'build.json',b);return out/'build.json'

def generate(identifier,out,n):
 from .generation import evaluate
 import numpy as np
 m=material(identifier);validate_selection(m);out=destination(out);fields,report=evaluate(m['id'],n)
 np.savez_compressed(out/'fields.npz',**fields);write_new(out/'generation.json',report)
 return out/'generation.json'

def build_experiment(name,out):
 from .experiments import experiment
 r,path=experiment(name);out=destination(out);source=resolve(r['source_scene'])
 run([os.environ.get('BLENDER_BIN','blender'),'-b','--threads','1','--python-exit-code','1','--python',str(ROOT/'materials/inspect_scene.py'),'--',str(source),r['source_sha256'],str(out/'inspection.json')],out/'inspect.log')
 d=read(out/'inspection.json')
 contract=r.get('contract',{})
 expected=contract.get('map_resolution',[4096,4096])
 if d['unit_scale']!=1:raise ContractError('Experiment SI unit contract failed')
 overrides={str(resolve(k)):v for k,v in contract.get('image_resolutions',{}).items()}
 for im in d['images']:
  if im['resolution']!=overrides.get(str(Path(im['path']).resolve()),expected):raise ContractError('Experiment native map dimensions failed: '+im['path'])
 b={'schema_version':1,'id':r['material_id'],'experiment_id':name,'experiment_manifest_sha256':sha(path),'selected':False,'quality_acceptance':False,'source':{'scene':str(source),'scene_sha256':r['source_sha256']},'inspection':d,'purpose':r['purpose'],'pass_gates':r['pass_gates']}
 write_new(out/'build.json',b);return out/'build.json'

def capture(build_path,out,resolution=None,samples=None,mode=None,fixed_samples=False,component_filter=None):
 if mode is not None and mode not in ('guided','unfiltered','component-recombined'):raise ContractError('Unsupported capture mode')
 if component_filter is not None and component_filter not in ('retain_transmission','filter_both'):raise ContractError('Unsupported component filter policy')
 if samples is not None and samples<1:raise ContractError('Samples must be positive')
 if fixed_samples and samples is None:raise ContractError('--fixed-samples requires --samples')
 b=read(build_path)
 if 'experiment_id' in b:
  from .experiments import verify_experiment_build
  r=verify_experiment_build(b);m=material(r['material_id']);s={'scene':r['source_scene'],'scene_sha256':r['source_sha256'],'capture':r['capture']}
 else:m=verify_build(b);s=m['selected']
 j=dict(s['capture']);mode=mode if mode is not None else j.get('capture_mode','guided')
 if mode not in ('guided','unfiltered','component-recombined'):raise ContractError('Unsupported selected capture mode')
 policy=component_filter if component_filter is not None else j.get('component_filter_policy')
 if policy is not None:
  if mode!='component-recombined':raise ContractError('Component filter policy requires component-recombined mode')
  if policy not in ('retain_transmission','filter_both'):raise ContractError('Unsupported selected component filter policy')
  j['component_filter_policy']=policy
 engine=j.get('engine','CYCLES')
 if mode=='unfiltered' and engine not in ('CYCLES','BLENDER_WORKBENCH'):raise ContractError('Unfiltered mode supports Cycles or Workbench geometry only')
 if mode=='component-recombined' and engine!='CYCLES':raise ContractError('Component capture supports Cycles only')
 if fixed_samples and engine!='CYCLES':raise ContractError('Fixed sampling applies only to Cycles')
 out=destination(out)
 j.update(source=str(resolve(s['scene'])),source_sha256=s['scene_sha256'],output=str(out/'INTERNAL_PANEL.png'),color_depth=16,dither_intensity=0,save_guides=not(mode=='unfiltered' and engine=='BLENDER_WORKBENCH'),capture_mode=mode,rss_cap_mib=j.get('rss_cap_mib',3072))
 if resolution:j.update(width=resolution,height=resolution)
 if samples is not None:j['samples']=samples
 if fixed_samples:j['use_adaptive_sampling']=False
 write_new(out/'job.json',j)
 pipeline=resolve('shared/cloud-eevee')
 run([sys.executable,str(pipeline/'run_panel_queue.py'),str(out/'job.json')],out/'capture.log')
 previous_path=sys.path[:]
 try:
  sys.path.insert(0,str(pipeline))
  spec=importlib.util.spec_from_file_location('material_capture_guard',pipeline/'cache_validation.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
 finally:sys.path[:]=previous_path
 final_name={'guided':'CLEAN_OIDN_PANEL.png','unfiltered':'UNFILTERED_PANEL.png','component-recombined':'COMPONENT_RECOMBINED_PANEL.png'}[mode]
 result=mod.validate_cache(j,out/(final_name+'.json'),pipeline)
 if not result['valid']:raise ContractError('Post-capture verification failed: '+str(result['reasons']))
 if 'experiment_id' in b:verify_experiment_build(b)
 else:verify_build(b)
 receipt=read(out/(final_name+'.json'))
 report={'id':m['id'],'fresh_render':True,'build_manifest':str(Path(build_path).resolve()),'image_path':str(out/final_name),({'guided':'clean_image','unfiltered':'unfiltered_image','component-recombined':'component_image'}[mode]):str(out/final_name),'capture_mode':mode,'noise_acceptance':False,'render_receipt':str(out/(final_name+'.json')),'capture_verified':result,'render_seconds':receipt['render_seconds'],'resolution':receipt['resolution'],'underlying_source_maps':sorted({tuple(im['resolution']) for im in b['inspection']['images']}),'quality_acceptance':False,'status':'integration capture; noise/visual acceptance pending, never automatic promotion'}
 if mode=='component-recombined':report['component_filter_policy']=j.get('component_filter_policy','retain_transmission')
 report.update(experiment_id=b.get('experiment_id'),selected=False if 'experiment_id' in b else True)
 write_new(out/'capture_result.json',report);return report

def doctor():
 rows=[]
 for m in registry()['materials']:
  try:r=validate_selection(m,deep=True);rows.append({'id':m['id'],'selected_assets_valid':True,'procedural_rebuild':m['procedural_rebuild']})
  except (ContractError,OSError) as e:rows.append({'id':m['id'],'selected_assets_valid':False,'blocker':str(e)})
 return {'materials':rows,'selected_count':sum(r['selected_assets_valid'] for r in rows),'catalogue_count':len(rows),'free_disk_bytes':shutil.disk_usage(ROOT).free}

def smoke(out):
 out=destination(out);reports=[];start=time.time()
 for identifier in ['05_champagne_brass','32_concrete_polished']:
  generate(identifier,out/(identifier+'_fields'),64)
  b=build(identifier,out/(identifier+'_build'))
  reports.append(capture(b,out/(identifier+'_capture'),128,8))
 result={'status':'passed','fresh_selected_captures':reports,'elapsed_seconds':time.time()-start,'not_claimed':['visual realism acceptance','native4K render pixels','full4096 regenerated-byte parity','other33 material integration qualification']}
 write_new(out/'smoke_result.json',result);return result

def main():
 p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='cmd',required=True)
 s.add_parser('list');s.add_parser('doctor')
 q=s.add_parser('parity');q.add_argument('material');q.add_argument('--out',required=True)
 q=s.add_parser('experiment');q.add_argument('name');q.add_argument('--out',required=True)
 q=s.add_parser('recover-native');q.add_argument('material',choices=['05','32']);q.add_argument('--out',required=True);q.add_argument('--geometry')
 q=s.add_parser('recover-geometry');q.add_argument('material',choices=['32']);q.add_argument('--out',required=True)
 q=s.add_parser('reconstruct-brass');q.add_argument('--native',required=True);q.add_argument('--out',required=True)
 q=s.add_parser('finish-brass');q.add_argument('--build',required=True);q.add_argument('--admission',required=True)
 q=s.add_parser('qualify-brass');q.add_argument('--build',required=True)
 q=s.add_parser('reconstruct-concrete');q.add_argument('--native',required=True);q.add_argument('--out',required=True)
 for command in ['build','generate']:
  q=s.add_parser(command);q.add_argument('material');q.add_argument('--out',required=True)
  if command=='generate':q.add_argument('--resolution',type=int,default=64)
 q=s.add_parser('capture');q.add_argument('build');q.add_argument('--out',required=True);q.add_argument('--resolution',type=int);q.add_argument('--samples',type=int);q.add_argument('--mode',choices=['guided','unfiltered','component-recombined'],default=None);q.add_argument('--fixed-samples',action='store_true');q.add_argument('--component-filter',choices=['retain_transmission','filter_both'],default=None)
 q=s.add_parser('smoke');q.add_argument('--out',required=True)
 args=p.parse_args();exit_code=0
 try:
  if args.cmd=='list':result=[{'id':m['id'],'selected_version':m['selected']['version'] if m['selected'] else None,'quality':m['quality'],'protected':m['protected_appearance']} for m in registry()['materials']]
  elif args.cmd=='doctor':result=doctor()
  elif args.cmd=='parity':
   m=material(args.material);out=destination(args.out);exit_code=run([sys.executable,'-m','materials.native_parity',m['id'],str(out)],out/'parity.log',accepted=(0,2));result=read(out/'parity.json')
  elif args.cmd=='reconstruct-brass':
   from .rebuild_brass import prepare
   result=prepare(args.native,args.out)
  elif args.cmd=='reconstruct-concrete':
   from .recover_aggregate import rebuild_concrete_scene
   result=rebuild_concrete_scene(args.native,args.out)
  elif args.cmd=='qualify-brass':
   from .rebuild_brass import qualify_reconstruction
   result=qualify_reconstruction(args.build)
  elif args.cmd=='finish-brass':
   from .rebuild_brass import finish
   result=finish(args.build,args.admission)
  elif args.cmd=='recover-geometry':
   from .recover_aggregate import recover_concrete_geometry
   result=recover_concrete_geometry(args.out)
  elif args.cmd=='recover-native':
   if args.material=='32':
    if not args.geometry:raise ContractError('Concrete native recovery requires --geometry with verified GeometryHeight.npy')
    from .recover_aggregate import recover_concrete_native
    result=recover_concrete_native(args.out,args.geometry)
   else:
    out=destination(args.out);run([sys.executable,'-m','materials.recover_native',str(out)],out/'generation.log');result=read(out/'native_recovery.json')
  elif args.cmd=='experiment':result={'build':str(build_experiment(args.name,args.out))}
  elif args.cmd=='build':result={'build':str(build(args.material,args.out))}
  elif args.cmd=='generate':result={'generation':str(generate(args.material,args.out,args.resolution))}
  elif args.cmd=='capture':result=capture(args.build,args.out,args.resolution,args.samples,args.mode,args.fixed_samples,args.component_filter)
  else:result=smoke(args.out)
  print(json.dumps(result,indent=2))
 except (ContractError,OSError,ValueError) as e:
  print(json.dumps({'status':'failed','error':str(e)}),file=sys.stderr);return 1
 return exit_code
if __name__=='__main__':raise SystemExit(main())
