"""Freeze only complete ramp batches; leave an active later batch untouched."""
from pathlib import Path
import argparse,json,hashlib,zipfile,re,numpy as np
P=Path(__file__).resolve().parents[1];ap=argparse.ArgumentParser();ap.add_argument('--through',type=int,required=True);ap.add_argument('--revision',type=int,required=True);a=ap.parse_args();records={}
with zipfile.ZipFile(P/'Textile_suspended_spinning_r01.zip') as z:
 for name in z.namelist():
  if name!='PREPARATION_FREEZE.json':records[name]=z.read(name)
for f in sorted((P/'source').glob('*')):
 if f.suffix in {'.py','.cpp'}:records['suspended-spinning/source/'+f.name]=f.read_bytes()
for f in sorted(P.glob('*.md')):records['suspended-spinning/'+f.name]=f.read_bytes()
for f in sorted((P/'receipts').glob('*.json')):
 if not any(t in f.name for t in ['library','archive']):records['suspended-spinning/receipts/'+f.name]=f.read_bytes()
accepted=[];summaries=[]
for i in range(1,a.through+1):
 folder=f'ramp{i:02}';summary=json.loads((P/folder/'receipts/progress.json').read_text());assert summary['status']!='running';summaries.append(summary)
 for f in sorted((P/folder).rglob('*')):
  if f.is_file() and f.suffix in {'.npz','.json','.log'}:records['suspended-spinning/'+f.relative_to(P).as_posix()]=f.read_bytes()
 records['suspended-spinning/receipts/'+folder+'.log']=(P/'receipts'/f'{folder}.log').read_bytes()
 for f in sorted((P/folder/'receipts').glob('accepted*.json')):accepted.append(json.loads(f.read_text()))
latest=Path(summaries[-1]['latest_accepted_state']);last=accepted[-1];initial=json.loads((P/'receipts/qualified_suspended_restart.json').read_text());work=sum(h['step_work_J'] for h in accepted);de=last['stored_energy_J']-initial['stored_energy_J'];ledger=dict(accepted_steps=len(accepted),last_state=str(latest.relative_to(P)),last_state_sha256=hashlib.sha256(latest.read_bytes()).hexdigest(),angle_rad=last['twist_rad'],turns_per_m=last['twist_turns_per_m'],omega_R=last['dimensionless_root_twist'],force_residual=last['force_residual'],stability=last['stability'],cumulative_boundary_work_J=work,stored_change_J=de,relative_work_gap=abs(work-de)/max(abs(work),abs(de),1e-18),ramp_seconds=sum(h['elapsed_s'] for h in summaries),maximum_state_stock_error=max(h['relative_stock_error'] for h in accepted),maximum_state_penetration_m=max(h['capsule_penetration_m'] for h in accepted),all_committed_gates_pass=all(h['history_accepted'] for h in accepted),scope='Fixed-stock suspended normal-contact control history at 128 segments per fibre. Finite friction and spatial convergence remain open.')
ledger_bytes=json.dumps(ledger,indent=2).encode();(P/'receipts'/f'history_through_{a.through:02}.json').write_bytes(ledger_bytes);records[f'suspended-spinning/receipts/history_through_{a.through:02}.json']=ledger_bytes
# Exact old comparison inputs, without copying its Blender scene.
old=P.parent/'spinning-carriage/ramp/batch3'
for part in ['data/accepted_004.npz','receipts/accepted_004.json']:records['spinning-carriage/ramp/batch3/'+part]=(old/part).read_bytes()
for f in sorted((P/'visible-candidate').rglob('*')):
 if f.is_file() and f.suffix in {'.py','.json','.blend','.log'} and '__pycache__' not in f.parts:records['suspended-spinning/'+f.relative_to(P).as_posix()]=f.read_bytes()
for folder,base in [(P/'resolution64',P.parent),(P.parent/'finite-friction-prototype',P.parent),(P.parent/'contact-migration-gate',P.parent),(P.parent/'generalized-friction-gate',P.parent),(P.parent/'flexible-friction-bridge',P.parent),(P.parent/'flexible-friction-process',P.parent)]:
 for f in sorted(folder.rglob('*')):
  if f.is_file() and f.suffix in {'.py','.json','.npz','.md','.log'} and '__pycache__' not in f.parts:records[f.relative_to(base).as_posix()]=f.read_bytes()
f=P.parent/'compact-friction/source/friction_compact.py';records['compact-friction/source/friction_compact.py']=f.read_bytes()
manifest={k:dict(sha256=hashlib.sha256(v).hexdigest(),bytes=len(v)) for k,v in records.items()};m=json.dumps(manifest,indent=2,sort_keys=True).encode();name=f'FREEZE_R{a.revision:02}.json';(P/name).write_bytes(m);out=P/f'Textile_suspended_spinning_r{a.revision:02}.zip'
with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for k,v in records.items():z.writestr(k,v)
 z.writestr(name,m)
r=dict(file=str(out),bytes=out.stat().st_size,sha256=hashlib.sha256(out.read_bytes()).hexdigest(),manifest_sha256=hashlib.sha256(m).hexdigest(),files=len(records),ledger=ledger);(P/'receipts'/f'archive_r{a.revision:02}_receipt.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
