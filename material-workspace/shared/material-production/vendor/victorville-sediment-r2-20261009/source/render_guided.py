"""Build the requested grazing-light views, capture Cycles guides, run CPU OIDN."""
import argparse,subprocess,os,json,hashlib,shutil
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--view',choices=['macro','detail','both'],default='both');p.add_argument('--resolution',type=int,default=960);p.add_argument('--samples',type=int,default=128);p.add_argument('--oidn-bin',default=os.environ.get('OIDN_BIN'));a=p.parse_args()
if not a.oidn_bin:raise SystemExit('Provide --oidn-bin pointing to official OpenImageDenoise oidnDenoise executable.')
R=Path(__file__).resolve().parents[1];maps=R/'native4096';out=R/'deliverables';out.mkdir(exist_ok=True);env=dict(os.environ,OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2',OIDN_BIN=str(Path(a.oidn_bin).resolve()))
(R/'receipts').mkdir(parents=True,exist_ok=True)
for view in ['macro','detail'] if a.view=='both' else [a.view]:
 log=R/(view+'_build.log')
 with log.open('w') as f:subprocess.run(['blender','-b','-t','2','--python-exit-code','1','--python',str(Path(__file__).resolve().parent/'build_scene.py'),'--',str(maps),str(out),view,str(a.samples),str(a.resolution),'build_only'],stdout=f,stderr=subprocess.STDOUT,env=env,check=True)
 source=out/f'CYBR_Victorville_{view}.blend';capture=R/'captures'/view;capture.mkdir(parents=True,exist_ok=True)
 job=dict(source=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),engine='CYCLES',width=a.resolution,height=a.resolution,samples=a.samples,capture_mode='guided',save_guides=True,light_tree=False,color_depth=16,dither_intensity=0,use_adaptive_sampling=False,rss_cap_mib=4608,time_limit_seconds=900,sync_limit_seconds=180,output=str(capture/'INTERNAL_PANEL.png'))
 jobpath=capture/'job.json';jobpath.write_text(json.dumps(job,indent=2))
 with (R/(view+'_guided.log')).open('w') as f:subprocess.run(['bash',str(R/'pipeline/clean_render.sh'),str(jobpath)],stdout=f,stderr=subprocess.STDOUT,env=env,check=True)
 final=out/f'CYBR_Victorville_Dust_Sand_R2_{view.title()}.png';shutil.copy2(capture/'CLEAN_OIDN_PANEL.png',final);shutil.copy2(capture/'CLEAN_OIDN_PANEL.png.json',R/'receipts'/f'{view}_guided_receipt.json');print(str(final),flush=True)
