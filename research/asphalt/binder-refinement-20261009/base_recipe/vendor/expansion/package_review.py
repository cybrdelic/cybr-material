"""Small portable source/scene/map bundles, with explicit unreviewed state."""
import hashlib,json,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'expansion'
GROUPS={
 'structures':['12_lattice_wire','13_lattice_expanded','15_future_ceramic','16_future_ribbed','17_blanket','24_carpet','27_grass'],
 'surfaces_displays':['18_lcd','19_crt','20_plastic_smooth','21_plastic_texture','25_tile_ceramic','26_tile_terrazzo','28_asphalt'],
 'concrete_cement':['29_concrete_cast','30_concrete_board_formed','31_concrete_brushed','32_concrete_polished','33_concrete_exposed_aggregate','34_concrete_weathered','35_cement_paste'],
}
def package(name,ids,scenes=True,revision=None):
    files=list((ROOT/'source/expansion').glob('*.py'))+list((ROOT/'source/expansion').glob('*.json'))
    files += [ROOT/'source/surface_math.py']
    files += [OUT/p for p in ('README.md','INTEGRATION.md','BASE_HASHES.json','STATUS.json','geometry_checks.json','BARK_REFERENCE_AND_REVIEW.md') if (OUT/p).exists()]
    scene_paths={}
    for id in ids:
        files+=list((OUT/'maps'/id).glob('*'))
        if scenes:
            stem=f'{id}_hero_soft'+('_'+revision if revision else '')
            if not (OUT/'scenes'/(stem+'.blend')).exists():stem=f'{id}_hero_soft'
            files +=[OUT/'scenes'/(stem+'.blend'),OUT/'receipts'/(stem+'.json')]
            scene_paths[id]=f'expansion/scenes/{stem}.blend'
    assert all(p.is_file() for p in files)
    path=OUT/f'CYBR_{name}_portable_review.zip'
    with zipfile.ZipFile(path,'w',zipfile.ZIP_DEFLATED,6) as z:
        for p in files:z.write(p,p.relative_to(ROOT).as_posix())
        z.writestr('SNAPSHOT.json',json.dumps(dict(baseline='de4146ac759164e29323b152d37f64417c00c3fb',
            recipes=ids,source_scenes=scene_paths,status='Major overhaul CPU candidate; final acceptance pending matched visual review',
            cloud_use='Separate EEVEE adapter may preview geometry; raster transmission is not bulk optical evidence',
            source_scope='Additive module snapshot; original ten assets not modified',
            files={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files}),indent=2))
    with zipfile.ZipFile(path) as z:assert z.testzip() is None
    print(path.name,path.stat().st_size,len(files),hashlib.sha256(path.read_bytes()).hexdigest(),flush=True)
if __name__=='__main__':
    for name,ids in GROUPS.items():package(name,ids)
    package('opaque_PETG_native4K',['22_petg'],False)
