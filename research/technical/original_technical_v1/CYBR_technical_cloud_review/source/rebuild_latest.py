"""Rebuild all current editable studies into a NEW directory without rendering.

Usage: python3 source/rebuild_latest.py --blender /path/to/blender --output /new/rebuild
Requires a Blender 4.3+ executable. The supplied frozen scene/evidence files are never overwritten.
"""
import argparse,json,shutil,subprocess,sys
from pathlib import Path

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--blender',required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    src=Path(__file__).resolve().parents[1];out=a.output.resolve()
    if out.exists():raise SystemExit('Output must be a new directory; frozen evidence must not be overwritten.')
    out.mkdir(parents=True);shutil.copytree(src/'source',out/'source',ignore=shutil.ignore_patterns('__pycache__'));shutil.copytree(src/'expansion/maps',out/'expansion/maps')
    for folder in ('scenes','tests','evidence'):(out/folder).mkdir()
    def run(script,args=()):
        command=[a.blender,'-b','--factory-startup','--python-exit-code','1','--python',str(out/'source'/script)]
        if args:command+=['--',*args]
        print('BUILD',script,flush=True)
        with (out/'tests'/('replay_'+Path(script).stem+'.log')).open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
    ids=['12_lattice_wire','14_chainmail','15_future_ceramic','16_future_ribbed','18_lcd','19_crt','20_plastic_smooth','21_plastic_texture','22_petg','23_petg_transparent','25_tile_ceramic']
    run('technical_studies.py',['--only',*ids])
    for script in ['refine_display_lighting.py','refine_plastic_joint.py','ceramic_glaze_proof.py','add_physical_supports.py','make_support_details.py','fix_hex_layout.py','clear_petg_transmission_proof.py','clear_petg_backlit_control.py','fix_chain_tube_frames.py','petg_seam_proof.py','check_technical.py','check_revision5.py','check_petg_optics.py','check_control_consistency.py']:
        run(script)
    index={
        '12_lattice_wire':['12_lattice_wire_r4_hero.blend'],
        '14_chainmail':['14_chainmail_r6_hero.blend','14_chainmail_r6_detail.blend'],
        '15_future_ceramic':['15_future_ceramic_r5_hero.blend'],
        '16_future_ribbed':['16_future_ribbed_r4_hero.blend'],
        '18_lcd':['18_lcd_r5_hero.blend'],
        '19_crt':['19_crt_r5_hero.blend'],
        '20_plastic_smooth':['20_plastic_smooth_r5_hero.blend','20_plastic_smooth_r5_detail.blend'],
        '21_plastic_texture':['21_plastic_texture_r5_hero.blend','21_plastic_texture_r5_detail.blend'],
        '22_petg':['22_petg_r5_hero.blend','22_petg_r5_detail.blend','22_petg_r6_seamproof.blend'],
        '23_petg_transparent':['23_petg_transparent_r5_hero.blend','23_petg_transparent_r6_sideproof.blend','23_petg_transparent_r7_backlit_control.blend'],
        '25_tile_ceramic':['25_tile_ceramic_r4_hero.blend','25_tile_ceramic_r5_grazing.blend']}
    for names in index.values():
        for name in names:assert (out/'scenes'/name).is_file(),name
    (out/'LATEST_SCENES.json').write_text(json.dumps({'scenes':index,'render_status':'Built and checked only; not rendered','clear_petg_limit':'As-printed rough transmission; backlit control is explicitly diagnostic, not a clear-window claim'},indent=2)+'\n')
    print('REBUILD_COMPLETE',out,flush=True)

if __name__=='__main__':main()
