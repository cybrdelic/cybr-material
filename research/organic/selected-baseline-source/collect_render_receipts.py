"""Verify render provenance, never infer appearance approval from file presence."""
from pathlib import Path
import hashlib,json,shutil
from PIL import Image
ROOT=Path(__file__).resolve().parents[1];PRE=Path('/workspace/shared/cloud-eevee/previews')
expected={'11_bark_r8':'11_bark_r8_hero.blend','27_grass_r3':'27_grass_r3_hero.blend','17_blanket_r5':'17_blanket_r5_hero.blend','24_carpet_r5':'24_carpet_r5_hero.blend'}
out=ROOT/'renders';out.mkdir(exist_ok=True);rows=[]
for stem,scene in expected.items():
 digest=hashlib.sha256((ROOT/'scenes'/scene).read_bytes()).hexdigest()
 for p in sorted(PRE.glob(stem+'*.png')):
  rec=p.with_suffix('.png.json')
  if not rec.exists():continue
  d=json.loads(rec.read_text());assert d['source_sha256']==digest,(p,'scene mismatch')
  with Image.open(p) as im:
   im.verify()
  with Image.open(p) as im:size=list(im.size)
  assert size==d['resolution'];assert d['missing_images']==[]
  row={'render':p.name,'image_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'scene_sha256':digest,'resolution':size,'engine':d['engine'],'render_seconds':d['render_seconds'],'provenance_pass':True,'visual_acceptance':'see visual_QA.json'}
  shutil.copy2(p,out/p.name);shutil.copy2(rec,out/rec.name);rows.append(row)
(ROOT/'receipts/checked_render_receipts.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows,indent=2))
