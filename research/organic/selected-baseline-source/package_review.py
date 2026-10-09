"""Portable private review archive. Includes candidates and labeled failures."""
from pathlib import Path
import json,hashlib,zipfile
ROOT=Path(__file__).resolve().parents[1]
scenes=['11_bark_r8_hero.blend','27_grass_r3_hero.blend','17_blanket_r5_hero.blend','24_carpet_r5_hero.blend']
files=[ROOT/'README.md',ROOT/'organic_rebuild.patch']
for folder in ('source','data','receipts','snapshots','renders','failed_comparisons'):
 files.extend(p for p in (ROOT/folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix not in ('.log',))
files.extend(ROOT/'scenes'/n for n in scenes)
# Only baseline code needed by the new generators, never unchanged old scenes/maps.
files.extend(p for p in (ROOT/'inputs/CYBR_structures_portable_review/source').rglob('*') if p.is_file() and '__pycache__' not in p.parts)
files=sorted(set(files));manifest={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files if p.exists()}
assert len(manifest)==len(files)
name=ROOT/'CYBR_organic_structural_rebuild_review.zip'
with zipfile.ZipFile(name,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in files:z.write(p,str(p.relative_to(ROOT)))
 z.writestr('MANIFEST_SHA256.json',json.dumps(manifest,indent=2))
 baseline=json.loads((ROOT/'receipts/INPUT_SHA256.json').read_text())
 excluded=[p for p in baseline if p not in manifest]+['scenes/'+p.name for p in (ROOT/'scenes').glob('*.blend') if p.name not in scenes]
 z.writestr('PORTABLE_SCOPE.json',json.dumps({'scope':'Current four scenes plus historical source and failed-pixel evidence','excluded_workspace_paths':excluded},indent=2))
with zipfile.ZipFile(name) as z:
 assert z.testzip() is None
 assert len(z.namelist())==len(files)+2
report={'path':str(name),'bytes':name.stat().st_size,'sha256':hashlib.sha256(name.read_bytes()).hexdigest(),'entries':len(files)+2,'zip_integrity_pass':True,'visual_acceptance':'See receipts/visual_QA.json; this archive includes failed comparisons.'}
(ROOT/'package_receipt.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
