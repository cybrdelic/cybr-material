"""Rebuild, validate and package the complete CYBR suite."""
import argparse,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--skip-maps',action='store_true');p.add_argument('--skip-renders',action='store_true');p.add_argument('--blender',default='blender');args=p.parse_args()
def run(*cmd):
    print('BUILD_STAGE '+' '.join(str(v) for v in cmd),flush=True)
    subprocess.run([str(v) for v in cmd],cwd=ROOT,check=True)
for folder in ['materials','exports','blender','docs/images/after','path_traced/renders','path_traced/verification','build']:(ROOT/folder).mkdir(parents=True,exist_ok=True)
if not args.skip_maps:run(sys.executable,'source/generate_materials.py','--resolution','4096')
run(sys.executable,'source/pack_unity.py','--pipeline','hdrp')
run(sys.executable,'source/pack_unity.py','--pipeline','urp')
run(sys.executable,'source/verify_maps.py')
if not args.skip_renders:
    run(args.blender,'-b','--python','source/render_path_traced.py','--','--kind','macros','--size','1024','--samples','768','--min-samples','128','--threshold','.005')
    run(args.blender,'-b','--python','source/render_path_traced.py','--','--kind','scenes','--size','1280','--samples','768','--min-samples','128','--threshold','.008')
    run(args.blender,'-b','--python','source/render_path_traced.py','--','--kind','architecture','--size','2048','--samples','512','--min-samples','128','--threshold','.01','--denoise')
run(args.blender,'-b','--python','source/build_blender.py','--','--samples','768')
for file in ['CYBR_Cycles_Details.blend','CYBR_Cycles_Still_Lifes.blend','CYBR_Cycles_Architecture.blend']:run(args.blender,'-b','blender/'+file,'--python','source/verify_path_traced.py')
run(sys.executable,'source/package_suite.py')
run(sys.executable,'source/verify_deliverable.py')
run(sys.executable,'source/create_archive.py')
print('SUITE_BUILD_COMPLETE',flush=True)
