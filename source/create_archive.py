"""Package the complete suite, excluding caches, drafts and Git internals."""
import argparse,hashlib,json,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--output',default='CYBR_Material_V3.zip');args=p.parse_args()
OUT=ROOT/'build';OUT.mkdir(exist_ok=True)
files=[]
for file in sorted(ROOT.rglob('*')):
    if not file.is_file():continue
    rel=file.relative_to(ROOT)
    if any(part in rel.parts for part in ('.git','build','draft_materials','drafts','metadata','__pycache__')):continue
    if file.suffix in ('.log','.blend1','.pyc') or 'draft' in file.name.lower():continue
    if file.name=='SHA256SUMS.txt':continue
    files.append(file)
for name in ['README.md','OPEN_ME.html','docs/IMPORT_GUIDE.md','docs/CRITICAL_REVIEW.md','blender/CYBR_Material_Atelier.blend','blender/CYBR_Cycles_Details.blend','blender/CYBR_Cycles_Still_Lifes.blend','blender/CYBR_Cycles_Architecture.blend']:assert (ROOT/name).is_file(),name
def sha(file):
    result=hashlib.sha256()
    with file.open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):result.update(chunk)
    return result.hexdigest()
checksum=ROOT/'SHA256SUMS.txt';checksum.write_text('\n'.join(sha(file)+'  '+file.relative_to(ROOT).as_posix() for file in files)+'\n')
files.append(checksum);dest=OUT/Path(args.output).name
with zipfile.ZipFile(dest,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=3,allowZip64=True) as archive:
    for file in files:archive.write(file,arcname='cybr-material/'+file.relative_to(ROOT).as_posix())
with zipfile.ZipFile(dest) as archive:assert archive.testzip() is None
release=dict(file=dest.name,files=len(files),bytes=dest.stat().st_size,sha256=sha(dest))
(OUT/'release.json').write_text(json.dumps(release,indent=2)+'\n')
(OUT/'CYBR_Material_V3.zip.sha256').write_text(release['sha256']+'  '+dest.name+'\n')
print(json.dumps(release),flush=True)
