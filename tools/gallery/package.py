"""Assemble a static-host directory without a suite rebuild or release mutation."""
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'build/material-gallery'
if OUT.exists():
    shutil.rmtree(OUT)
OUT.mkdir(parents=True)
for directory in ('gallery','materials','path_traced','blender','docs','audit','source'):
    shutil.copytree(ROOT/directory,OUT/directory,ignore=shutil.ignore_patterns('__pycache__','.*'))
registry='material-workspace/shared/material-production/registry.json'
(OUT/registry).parent.mkdir(parents=True,exist_ok=True)
shutil.copyfile(ROOT/registry,OUT/registry)
for filename in ('SHA256SUMS.txt','README.md'):
    shutil.copyfile(ROOT/filename,OUT/filename)
(OUT/'evidence').mkdir(exist_ok=True)
shutil.copyfile(ROOT/'evidence/baseline_manifest.json',OUT/'evidence/baseline_manifest.json')
(OUT/'index.html').write_text('''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>CYBR MATERIAL</title><meta http-equiv="refresh" content="0;url=gallery/index.html"></head><body><a href="gallery/index.html">Open CYBR MATERIAL gallery</a></body></html>''')
(OUT/'.nojekyll').touch()
total=sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file())
assert total < 900*1024*1024, f'Gallery exceeds conservative Pages budget: {total}'
print(f'Static site assembled at {OUT}: {total/1048576:.1f} MiB; native PNGs copied unchanged.')
