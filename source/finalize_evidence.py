"""Retain byte hashes and actual render settings for baseline and proof PNGs."""
import json,hashlib,struct
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def image_record(p):
    header=p.read_bytes()[:26]
    w,h,bits,channels=struct.unpack('>IIBB',header[16:26])
    assert bits==16 and channels==2,p
    return dict(file=p.relative_to(ROOT).as_posix(),resolution=[w,h],bits=bits,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size)
baseline={p.stem:image_record(p) for p in (ROOT/'evidence').glob('*.png')}
assert len(baseline)==17
(ROOT/'evidence/baseline_manifest.json').write_text(json.dumps(dict(commit='d8b35a52c80667f88c336c9150f64f82cefcf54f',renders=baseline),indent=2)+'\n')
draft=json.loads((ROOT/'path_traced/draft_manifest.json').read_text())
controlled={}
for p in (ROOT/'path_traced/controlled_surfaces').glob('*.png'):
    r=draft[p.stem]
    r['file']=r['file'].replace('\\','/')
    assert r['file'].startswith('controlled_surfaces/') and not r['denoising']
    assert r['maximum_samples']==768 and r['minimum_samples']==128
    controlled[p.stem]=dict(r,**image_record(p))
assert len(controlled)==10
(ROOT/'path_traced/controlled_manifest.json').write_text(json.dumps(controlled,indent=2)+'\n')
print('EVIDENCE_VERIFIED / 17 original native PNGs and 10 matched surface proofs')
