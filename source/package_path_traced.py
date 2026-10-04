"""Collect verified Cycles renders into an offline gallery and render pack."""
import argparse
import hashlib
import html
import json
import re
import struct
import zipfile
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'path_traced'
SPECS = json.loads((ROOT / 'materials/manifest.json').read_text())['materials']
STUDIES = [('11_stone_and_timber', 'Stone and Timber'),
           ('12_leather_and_linen', 'Leather and Linen'),
           ('13_metal_and_ceramic', 'Metal and Ceramic')]


def inspect():
    manifest = json.loads((OUT / 'render_manifest.json').read_text())
    for file in (OUT / 'metadata').glob('*.json'):
        if '_draft' not in file.name:
            manifest[file.stem] = json.loads(file.read_text())
    expected = [s['id'] + '_detail' for s in SPECS] + [key for key, _ in STUDIES]
    for key in expected:
        record = manifest[key]
        file = OUT / record['file']
        assert file.parent.name == 'renders', key
        assert file.name == key + '.png', key
        assert record['engine'] == 'CYCLES', key
        assert record['denoising'] is False, key
        with file.open('rb') as stream:
            header = stream.read(26)
        assert header[:8] == b'\x89PNG\r\n\x1a\n', key
        width, height, bit_depth, color_type = struct.unpack('>IIBB', header[16:26])
        assert bit_depth == 16 and color_type == 2, (key, bit_depth, color_type)
        assert [width, height] == record['resolution'], key
        with Image.open(file) as image:
            image.verify()
        record['color_depth'] = '16-bit RGB PNG'
        record['sha256'] = hashlib.sha256(file.read_bytes()).hexdigest()
        record['bytes'] = file.stat().st_size
    clean = {key: manifest[key] for key in expected}
    (OUT / 'render_manifest.json').write_text(json.dumps(clean, indent=2) + '\n')
    return clean


STYLE = '''
:root{--paper:#efece4;--ink:#262824;--muted:#73766d;--line:#d3d1c7}*{box-sizing:border-box}
body{margin:0;background:var(--paper);color:var(--ink);font:14px/1.6 Arial,Helvetica,sans-serif}
a{color:inherit}button{font:inherit;cursor:pointer}header,main,footer{max-width:1440px;margin:auto;padding-left:4vw;padding-right:4vw}
header{padding-top:26px;padding-bottom:26px;display:flex;align-items:center;justify-content:space-between;border-bottom:1px solid var(--line)}
.logo{font:28px Georgia,serif;letter-spacing:.16em;text-decoration:none}.meta{font-size:10px;letter-spacing:.13em;text-transform:uppercase;color:var(--muted)}
header a:last-child{font-size:12px;text-decoration:none}h1,h2,h3{font-family:Georgia,serif;font-weight:400;line-height:1.15}
h1{font-size:clamp(38px,5vw,65px);margin:12px 0 16px}h2{font-size:34px;margin:0}h3{font-size:21px;margin:8px 0 2px}
.intro{display:flex;justify-content:space-between;align-items:end;padding:45px 0 28px;gap:30px}.intro p{color:var(--muted);max-width:410px;margin:0}
.hero{width:100%;display:block;aspect-ratio:4/3;object-fit:cover}.picture{border:0;padding:0;background:none;color:inherit;text-align:left;display:block;width:100%}
.picture img{width:100%;display:block}.caption{display:flex;justify-content:space-between;margin-top:12px;color:var(--muted);font-size:12px}
.studies{display:grid;grid-template-columns:1fr 1fr;gap:22px;margin-top:34px}.studies img{aspect-ratio:4/3;object-fit:cover}
.section{border-top:1px solid var(--line);margin-top:70px;padding-top:32px;display:flex;justify-content:space-between;align-items:baseline;gap:20px}
.grid{display:grid;grid-template-columns:repeat(5,1fr);gap:28px 16px;margin:30px 0 55px}.grid img{aspect-ratio:1}.grid .meta{margin-top:12px}
footer{border-top:1px solid var(--line);padding-top:22px;padding-bottom:28px;color:var(--muted);font-size:11px;display:flex;justify-content:space-between;gap:20px}
button:focus-visible,a:focus-visible{outline:2px solid #8d724b;outline-offset:5px}dialog{border:0;padding:0;background:var(--paper);color:var(--ink);width:min(96vw,1450px);max-width:96vw;max-height:96vh}dialog::backdrop{background:#111a}
.dialog-head{display:flex;align-items:center;justify-content:space-between;gap:15px;padding:12px 20px;border-bottom:1px solid var(--line)}.actions{display:flex;align-items:center;gap:14px}
.actions button{border:1px solid var(--line);padding:7px 12px;background:none;color:inherit}.actions a{font-size:12px}.image-frame{max-height:81vh;overflow:auto;background:#252824;display:flex;justify-content:center;align-items:flex-start}
.image-frame img{display:block;max-width:100%;max-height:81vh;object-fit:contain}.image-frame.native{display:block}.image-frame.native img{max-width:none;max-height:none;margin:auto}.close{font-size:21px!important;border:0!important;padding:3px 10px!important}
@media(max-width:900px){.grid{grid-template-columns:repeat(3,1fr)}.intro{display:block}.intro p{margin-top:18px}.intro{padding-top:30px}.section{margin-top:45px}.caption{display:block}.caption span{display:block}.dialog-head{flex-wrap:wrap}}
@media(max-width:560px){.grid{grid-template-columns:repeat(2,1fr);gap:22px 12px}.studies{grid-template-columns:1fr}h3{font-size:19px}.section{display:block}.section .meta{display:block;margin-top:9px}footer{display:block}.actions{gap:10px}.actions button{padding:7px}}
'''


def gallery():
    def card(key, name, meta, cls=''):
        return f'<button class="picture {cls}" data-image="renders/{key}.png" data-name="{html.escape(name, quote=True)}"><img src="renders/{key}.png" alt="{html.escape(name, quote=True)} — Cycles path-traced render" loading="lazy"><div class="meta">{meta}</div><h3>{html.escape(name)}</h3></button>'
    cards = ''.join(card(s['id'] + '_detail', s['name'].split(' / ')[0],
                         s['id'][:2] + ' / ' + s['family'] + ' / ' + str(round(s['preview_diameter_m'] * 1000)) + ' mm view') for s in SPECS)
    content = f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>CYBR / Cycles Studies</title><style>{STYLE}</style></head><body>
<header><a class="logo" href="OPEN_RENDERS.html">CYBR</a><a href="README.txt">Render notes ↗</a></header><main>
<div class="intro"><div><span class="meta">Cycles / Path-traced studies</span><h1>Natural detail &amp; wear.</h1></div><p>Ten surfaces in close-up and three still-life scenes. Open any render to inspect it at its native resolution or save the 16-bit PNG.</p></div>
<button class="picture" data-image="renders/11_stone_and_timber.png" data-name="Stone and Timber"><img class="hero" src="renders/11_stone_and_timber.png" alt="Path-traced travertine sculpture, marble tray, oak cup, walnut boards and plaster" fetchpriority="high"></button>
<div class="caption"><span>01 / Stone and Timber</span><span>Calacatta · Travertine · Walnut · Fumed oak · Lime plaster</span></div>
<div class="studies">{card('12_leather_and_linen','Leather and Linen','02 / Stitched hide, woven flax and aged brass')}{card('13_metal_and_ceramic','Metal and Ceramic','03 / Crazed glaze, brushed brass and worn oxide')}</div>
<div class="section"><h2>Surface studies</h2><span class="meta">10 close-ups / physical displacement</span></div><div class="grid">{cards}</div>
</main><footer><span>Original procedural PBR surfaces · Blender Cycles · AgX color</span><span>Native renders · No denoising, blur, or added grain · <a href="render_manifest.json">Render settings</a></span></footer>
<dialog id="viewer"><div class="dialog-head"><span id="title"></span><div class="actions"><button id="zoom" aria-pressed="false">View 100%</button><a id="save" download>Save PNG ↓</a><button class="close" id="close" aria-label="Close render">×</button></div></div><div class="image-frame" id="frame"><img id="image" alt=""></div></dialog>
<script>
const dialog=document.getElementById('viewer'),frame=document.getElementById('frame'),img=document.getElementById('image'),zoom=document.getElementById('zoom');
document.querySelectorAll('[data-image]').forEach(button=>button.addEventListener('click',()=>{{img.src=button.dataset.image;img.alt=button.dataset.name;document.getElementById('title').textContent=button.dataset.name;document.getElementById('save').href=button.dataset.image;frame.classList.remove('native');zoom.setAttribute('aria-pressed','false');zoom.textContent='View 100%';dialog.showModal()}}));
zoom.addEventListener('click',()=>{{const native=frame.classList.toggle('native');zoom.setAttribute('aria-pressed',String(native));zoom.textContent=native?'Fit image':'View 100%'}});
document.getElementById('close').addEventListener('click',()=>dialog.close());dialog.addEventListener('click',event=>{{if(event.target===dialog){{const box=dialog.getBoundingClientRect();if(event.clientX<box.left||event.clientX>box.right||event.clientY<box.top||event.clientY>box.bottom)dialog.close()}}}});
</script></body></html>'''
    (OUT / 'OPEN_RENDERS.html').write_text(content, encoding='utf-8')


def contacts():
    serif=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSerif.ttf',54)
    sans=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',18)
    small=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',14)
    canvas=Image.new('RGB',(2000,1180),'#efece4');draw=ImageDraw.Draw(canvas)
    draw.text((40,35),'C Y B R  M A T E R I A L',font=serif,fill='#262824')
    draw.text((40,120),'CYCLES / PATH-TRACED SURFACE STUDIES',font=small,fill='#73766d')
    for i,spec in enumerate(SPECS):
        x=40+(i%5)*392;y=185+(i//5)*450
        with Image.open(OUT/'renders'/(spec['id']+'_detail.png')) as source:
            canvas.paste(source.convert('RGB').resize((360,360),Image.Resampling.LANCZOS),(x,y))
        draw.text((x,y+374),spec['id'][:2]+' / '+spec['name'].split(' / ')[0],font=sans,fill='#262824')
        draw.text((x,y+407),str(round(spec['preview_diameter_m']*1000))+' MM VIEW  /  '+spec['family'].upper(),font=small,fill='#73766d')
    draw.text((40,1140),'Actual 3D path tracing · Native 4K procedural maps · Physical relief · No denoising or added grain',font=small,fill='#73766d')
    canvas.save(OUT/'CYBR_Cycles_Detail_Collection.jpg',quality=97,subsampling=0)
    canvas=Image.new('RGB',(1280,1715),'#efece4');draw=ImageDraw.Draw(canvas)
    draw.text((32,27),'C Y B R',font=serif,fill='#262824')
    draw.text((32,110),'CYCLES / THREE STILL-LIFE STUDIES',font=small,fill='#73766d')
    for i,(key,title) in enumerate(STUDIES):
        y=165+i*510
        with Image.open(OUT/'renders'/(key+'.png')) as source:
            canvas.paste(source.convert('RGB').resize((608,456),Image.Resampling.LANCZOS),(32,y))
            canvas.paste(source.convert('RGB').crop((320,240,960,720)).resize((552,414),Image.Resampling.LANCZOS),(680,y+21))
        draw.text((32,y+465),f'{i+1:02d} / '+title,font=sans,fill='#262824')
        draw.text((680,y+465),'DETAIL / CENTER CROP',font=small,fill='#73766d')
    canvas.save(OUT/'CYBR_Cycles_Scenes.jpg',quality=97,subsampling=0)


def archive(manifest):
    readme='''CYBR MATERIAL 3 / Cycles Path-Traced Renders

Open path_traced/OPEN_RENDERS.html in a browser.
The renders folder contains 10 close-ups (1024 x 1024) and 3 still lifes
(1280 x 960), all native 16-bit RGB PNGs rendered by Blender Cycles.
AgX Medium High Contrast color; no denoising, image blur, sharpening,
added film grain, or image upscaling. Contact sheets are display derivatives.
Actual settings and checksums are in render_manifest.json.

The Blender scenes are an add-on to CYBR MATERIAL 3. Put the blender and
path_traced directories beside the existing materials directory in the
extracted CYBR MATERIAL 3 suite. The scenes use relative 4K texture paths.
Open blender/CYBR_Cycles_Details.blend or CYBR_Cycles_Still_Lifes.blend.
Select a scene from Blender's scene selector and press F12 to render.

The complete suite archive includes these scenes, renders, and all textures.
Materials are original procedural surfaces, not scans.
'''
    (OUT/'README.txt').write_text(readme)
    files=[p for p in OUT.rglob('*') if p.is_file() and 'drafts' not in p.parts
           and 'metadata' not in p.parts and p.name!='SHA256SUMS.txt']
    files += [ROOT/'blender/CYBR_Cycles_Details.blend', ROOT/'blender/CYBR_Cycles_Still_Lifes.blend']
    files += [ROOT/'source/render_path_traced.py',ROOT/'source/verify_path_traced.py',ROOT/'source/package_path_traced.py']
    checksum=OUT/'SHA256SUMS.txt'
    checksum.write_text('\n'.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.relative_to(ROOT).as_posix() for p in sorted(files))+'\n')
    files.append(checksum)
    dest=ROOT.parent/'CYBR_Cycles_Path_Traced_Renders.zip'
    with zipfile.ZipFile(dest,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=3) as zipped:
        for p in sorted(files):zipped.write(p,arcname=ROOT.name+'/'+p.relative_to(ROOT).as_posix())
    with zipfile.ZipFile(dest) as zipped:assert zipped.testzip() is None
    print(json.dumps(dict(archive=str(dest),files=len(files),size_MiB=round(dest.stat().st_size/1024**2,2),renders=len(manifest)),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--gallery-only',action='store_true')
    args=parser.parse_args()
    gallery()
    if not args.gallery_only:
        manifest=inspect();contacts();archive(manifest)
