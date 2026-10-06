"""Copy completed rendered attempts and verify frozen source receipts without altering pixels."""
from pathlib import Path
import json,hashlib,shutil
ROOT=Path(__file__).resolve().parents[1]
PREV=Path('/workspace/shared/cloud-eevee/previews')
ids={'12','14','15','16','18','19','20','21','22','23','25'}
dst=ROOT/'evidence'/'attempts';dst.mkdir(parents=True,exist_ok=True);report=[]
for p in sorted(PREV.glob('*.png')):
    if p.name.split('_')[0] not in ids:continue
    receipt=p.with_suffix(p.suffix+'.json')
    if not receipt.exists():continue
    data=json.loads(receipt.read_text());source=Path(data['job']['source']);actual=hashlib.sha256(source.read_bytes()).hexdigest();assert actual==data['source_sha256'],source
    shutil.copy2(p,dst/p.name);shutil.copy2(receipt,dst/receipt.name)
    report.append({'image':'evidence/attempts/'+p.name,'image_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'source_scene':source.name,'source_sha256':actual,'source_frozen_hash_matches_receipt':True,'engine':data['engine'],'resolution':data['resolution'],'render_seconds':data['render_seconds']})
(ROOT/'tests/rendered_attempts.json').write_text(json.dumps(report,indent=2)+'\n');print('EVIDENCE_COLLECTED',len(report),flush=True)
