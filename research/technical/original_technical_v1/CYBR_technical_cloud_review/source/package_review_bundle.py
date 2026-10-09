"""Package this user-facing review and its frozen evidence; no Library writes."""
import argparse,hashlib,json,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.parent.mkdir(parents=True,exist_ok=True)
files=[ROOT/'REVIEW.md',ROOT/'LATEST_SCENES.json']
for name in ('source','scenes','tests','evidence'):
    files.extend(q for q in (ROOT/name).rglob('*') if q.is_file() and '__pycache__' not in q.parts and q.suffix!='.blend1')
ids=['12_lattice_wire','14_chainmail','15_future_ceramic','16_future_ribbed','18_lcd','19_crt','20_plastic_smooth','21_plastic_texture','22_petg','23_petg_transparent','25_tile_ceramic']
for id in ids:files.extend(q for q in (ROOT/'expansion/maps'/id).rglob('*') if q.is_file())
files=sorted(set(files));manifest={str(q.relative_to(ROOT)):hashlib.sha256(q.read_bytes()).hexdigest() for q in files}
mp=ROOT/'SHA256_MANIFEST.json';mp.write_text(json.dumps(manifest,indent=2)+'\n');files.append(mp)
with zipfile.ZipFile(a.output,'w',zipfile.ZIP_DEFLATED,6) as z:
    for q in files:z.write(q,'CYBR_technical_cloud_review/'+str(q.relative_to(ROOT)))
with zipfile.ZipFile(a.output) as z:assert z.testzip() is None
print(json.dumps({'path':str(a.output),'bytes':a.output.stat().st_size,'files':len(files),'sha256':hashlib.sha256(a.output.read_bytes()).hexdigest()}),flush=True)
