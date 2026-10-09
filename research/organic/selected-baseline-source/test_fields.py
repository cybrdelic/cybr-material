"""Numerical field and immutable provenance tests; not an appearance score."""
from pathlib import Path
import hashlib,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
rows=[]
scope_path=ROOT/'PORTABLE_SCOPE.json';scope=json.loads(scope_path.read_text()) if scope_path.exists() else {};excluded=set(scope.get('excluded_workspace_paths',[]));skipped=[]
for rev in ('r4','r5'):
 p=ROOT/'data'/f'bark_{rev}.npz';a=np.load(p);h=a['height'];c=a['color']
 assert h.shape==(641,641) and c.shape==(641,641,3)
 assert np.isfinite(h).all() and np.isfinite(c).all()
 assert 0<h.min()<h.max()<.02 and 0<c.min()<c.max()<1
 assert h.std()>.0004
 rows.append({'field':p.name,'pass':True,'height_min_m':float(h.min()),'height_max_m':float(h.max()),'height_std_m':float(h.std()),'content_sha256':hashlib.sha256(h.tobytes()+c.tobytes()).hexdigest()})
baseline=json.loads((ROOT/'receipts/INPUT_SHA256.json').read_text())
for path,digest in baseline.items():
 p=ROOT/path
 if not p.exists() and path in excluded:skipped.append(path);continue
 assert hashlib.sha256(p.read_bytes()).hexdigest()==digest,path
for p in (ROOT/'receipts').glob('*_build.json'):
 d=json.loads(p.read_text());short={'11_bark':'bark','17_blanket':'fleece','27_grass':'grass','24_carpet':'carpet'}[d['recipe']]
 snap=ROOT/'snapshots'/f'{short}_{d["revision"]}'/'source'
 for name,hsh in d['source_sha256'].items():assert hashlib.sha256((snap/name).read_bytes()).hexdigest()==hsh,(short,name)
 scene=ROOT/'scenes'/Path(d['scene']).name
 if not scene.exists() and str(scene.relative_to(ROOT)) in excluded:skipped.append(str(scene.relative_to(ROOT)))
 else:assert hashlib.sha256(scene.read_bytes()).hexdigest()==d['scene_sha256']
report={'fields':rows,'frozen_baseline_files':len(baseline),'available_baseline_unchanged':True,'available_scene_and_all_snapshot_hashes_match':True,'excluded_from_portable_subset':skipped}
(ROOT/'receipts/source_tests.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
