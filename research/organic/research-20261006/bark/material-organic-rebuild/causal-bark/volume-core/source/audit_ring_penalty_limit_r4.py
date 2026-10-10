from pathlib import Path
import json,time,resource
from audit_ring_penalty_mesh_r4 import solve
from source_receipt_r3 import hashes
R=Path(__file__).resolve().parents[1];start=time.time();dep=hashes();rows=[];failure=None
try:
 for n in (16,32,64):
  for eta in (.005,.0025,.00125):
   row=solve(n,eta);rows.append(row);print('CASE',json.dumps({k:v for k,v in row.items() if k!='trace'}),flush=True)
except Exception as e:failure=type(e).__name__+': '+str(e)
out={'completed':failure is None,'failure':failure,'cases':rows,'source_sha256':dep,'dependencies_unchanged':dep==hashes(),'elapsed_s':time.time()-start,'peak_rss_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'scope':'Intact ring numerical-continuity limit. No fracture path or natural-cork calibration.'};(R/'receipts/ring_penalty_limit_r4.json').write_text(json.dumps(out,indent=2));print('RESULT',json.dumps({k:v for k,v in out.items() if k not in ('source_sha256','cases')}),flush=True)
