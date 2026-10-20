"""Build the static gallery from existing evidence; never regenerate source assets.

Run from any directory: python tools/gallery/build.py
Pillow is used only for explicitly labeled, downsampled display derivatives.
"""
from pathlib import Path
import hashlib
import html
import json

from PIL import Image

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'gallery'
REGISTRY = 'material-workspace/shared/material-production/registry.json'
PREVIEW_PIPELINE = 'webp-v2-gray16-unorm'
PREVIOUS_PREVIEWS = {}


def collect_preview_cache(value):
    if isinstance(value, dict):
        if all(key in value for key in ('path', 'sha256', 'previews')):
            PREVIOUS_PREVIEWS[value['path']] = value
        for child in value.values():
            collect_preview_cache(child)
    elif isinstance(value, list):
        for child in value:
            collect_preview_cache(child)


if (OUT / 'catalog.json').exists():
    previous_catalog = json.loads((OUT / 'catalog.json').read_text())
    if previous_catalog.get('preview_pipeline') == PREVIEW_PIPELINE:
        collect_preview_cache(previous_catalog)
CHANNELS = {
    'BaseColor': ('Base color', 'sRGB · RGB8'),
    'Roughness': ('Roughness', 'Linear · Gray8'),
    'Metallic': ('Metallic', 'Linear · Gray8'),
    'Normal_OpenGL': ('Normal +Y', 'Tangent space · RGB8'),
    'Normal_DirectX': ('Normal −Y', 'Tangent space · RGB8'),
    'Height': ('Height', 'Linear · Gray16'),
    'Height_Macro': ('Macro height', 'Linear · Gray16'),
    'Normal_Micro_OpenGL': ('Residual normal +Y', 'Tangent space · RGB8'),
    'AO': ('Ambient occlusion', 'Linear · Gray8'),
    'ORM': ('Packed ORM', 'Linear · R: AO / G: roughness / B: metallic'),
    'WearMask': ('Wear mask', 'Linear · Gray8'),
    'Opacity': ('Opacity', 'Linear · Gray8'),
}


def load(path):
    return json.loads((ROOT / path).read_text())


def esc(value):
    return html.escape(str(value), quote=True)


def asset(path, expected=None):
    p = ROOT / path
    raw = p.read_bytes()
    sha = hashlib.sha256(raw).hexdigest()
    if expected and expected != sha:
        raise ValueError(f'Accepted evidence hash changed: {path}')
    im = Image.open(p)
    bits = raw[24] if raw[:8] == b'\x89PNG\r\n\x1a\n' else 8
    return dict(path=path, sha256=sha, bytes=len(raw), width=im.width, height=im.height, bit_depth=bits)


def preview(a):
    """Make two bounded derivatives. Native links always point at original bytes."""
    previous = PREVIOUS_PREVIEWS.get(a['path'])
    if previous and previous['sha256'] == a['sha256'] and all((OUT / p).is_file() for p in previous['previews']):
        a['previews'] = previous['previews']
        return a
    key = hashlib.sha256(a['path'].encode()).hexdigest()[:16]
    outputs = []
    with Image.open(ROOT / a['path']) as source:
        if a['bit_depth'] == 16 and source.mode in ('I;16', 'I;16B', 'I'):
            display = source.convert('I').point(lambda v: v / 257).convert('RGB')
        else:
            display = source.convert('RGB')
        display.thumbnail((800, 800), Image.Resampling.LANCZOS)
    for width in (360, 800):
        dst = OUT / 'assets' / f'{key}-{width}.webp'
        im = display.copy()
        im.thumbnail((width, width), Image.Resampling.LANCZOS)
        im.save(dst, 'WEBP', quality=82, method=4)
        outputs.append('assets/' + dst.name)
    a['previews'] = outputs
    return a


def image(a, prefix='', alt='', eager=False):
    small, large = a['previews']
    return (f'<img src="{prefix}{small}" srcset="{prefix}{small} 360w, {prefix}{large} {min(800,a["width"])}w" '
            f'sizes="(max-width: 600px) 92vw, (max-width: 1000px) 46vw, 30vw" '
            f'width="{a["width"]}" height="{a["height"]}" alt="{esc(alt)}" '
            f'loading="{"eager" if eager else "lazy"}" decoding="async">')


def link(path, label, prefix='../'):
    return f'<a href="{prefix}{esc(path)}">{esc(label)}</a>'


def page(title, body, prefix=''):
    return f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="description" content="Inspect original CYBR MATERIAL procedural maps, physical scale and native Blender Cycles evidence. Individual images, no archive required.">
<title>{esc(title)} · CYBR MATERIAL</title><link rel="icon" href="{prefix}icon.svg" type="image/svg+xml"><link rel="stylesheet" href="{prefix}gallery.css">
<script src="{prefix}gallery.js" defer></script></head><body>
<a class="skip" href="#main">Skip to content</a>
<header><a class="wordmark" href="{prefix}index.html">CYBR <span>MATERIAL</span></a>
<nav aria-label="Main"><a href="{prefix}index.html#collection">Collection</a><a href="{prefix}index.html#selected">Selected studies</a><a href="{prefix}index.html#register">Registry</a><a href="{prefix}index.html#evidence">Evidence</a></nav></header>
<main id="main">{body}</main>
<footer><strong>Made from structure.</strong><p>Original procedural work. Actual Cycles rendering. Physical values are artistic estimates, not measured material calibration.</p>
<a href="https://github.com/cybrdelic/cybr-material">Source repository ↗</a> · <a href="{prefix}catalog.json">Image dimensions, settings &amp; SHA-256</a></footer></body></html>'''


def view(label, a, prefix='../', note=''):
    path = prefix + '../' + a['path']
    return f'''<figure class="view"><a href="{esc(path)}">{image(a,prefix,label)}</a>
<figcaption><strong>{esc(label)}</strong><span>{a['width']} × {a['height']} · {a['bit_depth']}-bit PNG{esc(note)}</span>
<a href="{esc(path)}">Open original image ↗</a></figcaption></figure>'''


def build():
    for sub in ('assets', 'materials', 'selected'):
        (OUT / sub).mkdir(parents=True, exist_ok=True)
    manifest = load('materials/manifest.json')
    renders = load('path_traced/render_manifest.json')
    proofs = load('path_traced/controlled_manifest.json')
    registry = load(REGISTRY)
    # Pin original evidence to the input checkpoint rather than a mutable release.
    baseline = '365bd17aed0e3bbe9d2e7896f2bfdce1b9b54a92'
    catalog = dict(schema='cybr.material.gallery.v1', preview_pipeline=PREVIEW_PIPELINE, baseline_commit=baseline,
                   baseline_source='https://github.com/cybrdelic/cybr-material/tree/'+baseline,
                   native_map_resolution=[4096,4096], materials=[], selected=[], scenes=[],
                   registry_path=REGISTRY, registry_sha256=hashlib.sha256((ROOT/REGISTRY).read_bytes()).hexdigest())
    cards = []
    direct_rows = []
    direct_details = []
    for spec in manifest['materials']:
        mid = spec['id']; key = mid + '_detail'
        beauty = preview(asset('path_traced/'+renders[key]['file'], renders[key]['sha256']))
        proof = preview(asset(proofs[key]['file'], proofs[key]['sha256']))
        maps = []
        for channel, (label, interpretation) in CHANNELS.items():
            a = asset(f'materials/{mid}/{channel}.png')
            assert [a['width'], a['height']] == [4096,4096]
            if channel.startswith('Height'):
                assert a['bit_depth'] == 16
            a.update(channel=channel, label=label, interpretation=interpretation)
            maps.append(preview(a))
        entry = dict(**spec, maps=maps, beauty=beauty, proof=proof,
                     render_settings=renders[key], proof_settings=proofs[key],
                     qualification='Rebuilt suite appearance; procedural approximation, not measured calibration')
        catalog['materials'].append(entry)
        name = spec['name'].split(' / ')[0]
        tile = spec['tile_m']*1000
        finite = spec.get('tiling') == 'finite-cut'
        cards.append(f'''<article class="card" data-family="{esc(spec['family'])}"><a href="materials/{mid}.html">{image(beauty,alt=name+' — actual Cycles specimen')}<div class="card-body"><span class="eyebrow">{esc(spec['family'])} / rebuilt suite</span><h3>{esc(name)} ↗</h3><p>4096² maps · {tile:g} mm {'finite cut' if finite else 'UV tile'}<br>Cycles specimen: {beauty['width']} × {beauty['height']}</p></div></a></article>''')
        channels_html = ''.join(f'''<article class="channel"><a href="../../{a['path']}">{image(a,'../',name+' — '+a['label'])}</a><h3>{a['label']}</h3><p>{a['interpretation']}<br>4096 × 4096 · {a['bytes']/1048576:.2f} MiB</p><a href="../../{a['path']}">Open native 4K PNG ↗</a><details><summary>SHA-256</summary><code>{a['sha256']}</code></details></article>''' for a in maps)
        settings=renders[key]
        body=f'''<div class="breadcrumb"><a href="../index.html#collection">Collection</a> / {esc(spec['family'])}</div>
<div class="material-heading"><p class="eyebrow">Rebuilt suite / {mid[:2]}</p><h1>{esc(name)}</h1><p class="intro">{esc(spec['description'])}</p></div>
<div class="facts"><div><b>4096 × 4096</b><span>Native texture maps</span></div><div><b>{tile:g} mm</b><span>{'Finite wood cut; not certified seamless' if finite else 'Physical width of one UV tile'}</span></div><div><b>{spec['height_scale_m']*1000:g} mm</b><span>Height displacement scale</span></div></div>
<section><div class="section-heading"><h2>Rendered appearance</h2><span>Actual 3D geometry · Blender Cycles</span></div><p>The specimen and matched surface proof are different presentations. Native texture resolution does not determine render resolution.</p><div class="views">{view('Fabricated specimen',beauty)}{view('Matched surface proof',proof,note=f" · {proofs[key]['inspection_width_m']*1000:g} mm inspection width")}</div>
<details><summary>Renderer settings &amp; source evidence</summary><p>Specimen: Blender {settings['blender']} · {settings['device']} · {settings['maximum_samples']} maximum / {settings['minimum_samples']} minimum samples · adaptive threshold {settings['adaptive_threshold']:.4g} · denoising: {str(settings['denoising']).lower()} · {esc(settings['view_transform'])}.</p><p>{link('path_traced/render_manifest.json','Specimen manifest','../../')} · {link('path_traced/controlled_manifest.json','Matched-proof manifest','../../')} · {link('blender/CYBR_Cycles_Details.blend','Editable specimen scenes','../../')} · {link('source/render_path_traced.py','Renderer source','../../')} · {link('docs/STUDY_REPAIR_REVIEW.md','Review and remaining limits','../../')}</p><p>Specimen SHA-256 <code>{beauty['sha256']}</code></p><p>Proof SHA-256 <code>{proof['sha256']}</code></p></details></section>
<section id="channels"><div class="section-heading"><h2>Native PBR channels</h2><span>12 individual PNGs · no ZIP required</span></div><p>Base color uses sRGB; data maps use linear / Non-Color. Displacement in meters = (height − 0.5) × {spec['height_scale_m']:g}. Pair macro height with residual normal for geometric displacement, or use the full normal without macro displacement.</p><p>Previews below are reduced display copies. Each link opens the unchanged native 4096 × 4096 map.</p><div class="channels">{channels_html}</div></section>
<section><h2>Scope of this appearance</h2><p>{'The exported wood maps are finite cuts. The Blender wood shader uses a separate 3D growth volume with per-object saw planes and coherent end grain; the maps do not encode that volume.' if finite else 'Procedural anatomy and finish models use artistic physical estimates. These are not measured scans or calibrated reflectance data.'}</p><p>{link(f'materials/{mid}/material.json','Physical material metadata','../../')} · {link(f'docs/images/comparisons/{key}.jpg','Original appearance comparison','../../')} · {link('evidence/baseline_manifest.json','Preserved baseline evidence','../../')}</p></section>'''
        (OUT/'materials'/f'{mid}.html').write_text(page(name,body,'../'))
        direct_rows.append(f'| {name} | [Individual panel](../gallery/materials/{mid}.html) | [4096² base color](../materials/{mid}/BaseColor.png) | [Native Cycles specimen](../{beauty["path"]}) · {beauty["width"]} × {beauty["height"]} |')
        channel_links = ' · '.join(f'[{a["label"]}](../{a["path"]})' for a in maps)
        direct_details.append(f'### {name}\n\nPhysical UV width: {tile:g} mm. Height scale: {spec["height_scale_m"]*1000:g} mm. {"Finite wood cut, not certified seamless." if finite else "Physical dimensions are artistic estimates."}\n\n{channel_links}\n\n[Matched surface proof](../{proof["path"]}) · {proof["width"]} × {proof["height"]}; {proofs[key]["inspection_width_m"]*1000:g} mm inspection width.\n')

    # Exact reviewed references; archive-only studies stay in the registry below.
    selected_specs = [
        ('21_plastic_texture','Molded stipple plastic','audit/previews/21_molded_stipple_plastic_full.png',
         'ec921941239f9954e6440cac0d7f48a5d305e1446f63ec054e131978ab32ab27',960,
         '80 × 80 × 8 mm complete part. Reviewed curved grain and edge highlight; physical calibration remains false.'),
    ]
    selected_cards=[]
    for mid,name,path,sha,native,description in selected_specs:
        a=preview(asset(path,sha));m=next(m for m in registry['materials'] if m['id']==mid)
        entry=dict(id=mid,name=name,image=a,native_render_resolution=[native,native],
                   labeled_presentation_resolution=[a['width'],a['height']],
                   review=m['preferred_visual_reference'],description=description)
        extra=''
        catalog['selected'].append(entry)
        selected_cards.append(f'''<article><h3>{esc(name)}</h3><p>{esc(description)}</p><p>Native render area {native} × {native}; labeled presentation {a['width']} × {a['height']}. Original annotations are preserved. Archived scenes/maps are required for replay; no standalone PBR set is implied.</p><div class="views">{view('Reviewed complete panel',a,'',note=f' · render area {native} × {native}')}{extra}</div><details><summary>Recorded source, capture &amp; checksum</summary><p>{link(REGISTRY,'Frozen registry')} · {link('audit/ARTIFACTS.json','Archive recovery references')}</p><pre>{esc(json.dumps(m['preferred_visual_reference'],indent=2))}</pre><code>{a['sha256']}</code></details></article>''')

    scene_cards=[]
    for key,settings in renders.items():
        if key.endswith('_detail') and key[:2] in [f'{i:02}' for i in range(1,11)]:continue
        a=preview(asset('path_traced/'+settings['file'],settings['sha256']))
        catalog['scenes'].append(dict(id=key,image=a,settings=settings))
        name=key[3:].replace('_',' ').title()
        scene_cards.append(f'<article class="card"><a href="../{a["path"]}">{image(a,alt=name)}<div class="card-body"><h3>{name} ↗</h3><p>Native Cycles render · {a["width"]} × {a["height"]}<br>{settings["maximum_samples"]} max samples · {"Guided OIDN" if settings["denoising"] else "No denoising"}</p></div></a></article>')

    rows=[]
    for m in registry['materials']:
        mid=m['id']; preferred=m.get('preferred_visual_reference',{})
        present=mid in [x[0] for x in selected_specs]
        experimental=mid.startswith(('24_','29_','30_','31_','32_','33_','34_','35_'))
        status='Experimental / appearance unqualified' if experimental else ('Reviewed appearance / image available' if present else 'Registry reference / archived inputs')
        rows.append(f'<tr><th scope="row">{esc(m["name"])}</th><td><span class="badge {"experimental" if experimental else ""}">{status}</span></td><td>{esc(m["quality"])}<br><small>{esc(m.get("dependencies_status",""))}</small></td></tr>')
    body=f'''<section class="hero"><p class="eyebrow">The surface collection / native evidence</p><h1>Structure.<br>Surface.<br><em>In full detail.</em></h1><p class="intro">Ten rebuilt procedural materials. Individual native 4K maps, physical scale and real Cycles appearances — ready to inspect, one image at a time.</p><div class="hero-links"><a class="button" href="#collection">Explore the materials ↓</a><a href="#selected">Reviewed engineering studies ↗</a></div><div class="hero-meta"><span>10 native 4K map sets</span><span>120 individual channels</span><span>27 original Cycles views</span></div></section>
<section id="collection"><div class="section-heading"><h2>The rebuilt collection</h2><span>01 — 10 / fabrication studies</span></div><p>Each material opens its own panel with all twelve channels, the fabricated specimen, the matched surface proof, physical scale and original checksums. Texture maps are 4096 × 4096. Specimen and proof renders are 1024 × 1024.</p><div class="filters" hidden><label for="family">Material family</label><select id="family"><option value="all">All materials</option>{''.join(f'<option>{esc(f)}</option>' for f in sorted(set(s['family'] for s in manifest['materials'])))}</select><span id="result-count" role="status">10 materials</span></div><div class="grid" id="material-grid">{''.join(cards)}</div></section>
<section id="selected"><div class="section-heading"><h2>Selected engineering studies</h2><span>Separate reviewed appearances</span></div><p>These exact reviewed images come from the engineering registry. They do not replace the rebuilt collection or establish measured material physics.</p>{''.join(selected_cards)}</section>
<section id="scenes"><div class="section-heading"><h2>Materials in context</h2><span>Still lifes / architecture</span></div><p>Original native renders; no upscaling. Open an image to inspect its full resolution. Camera, sample and denoising settings remain in the {link('path_traced/render_manifest.json','render manifest')}.</p><div class="grid">{''.join(scene_cards)}</div></section>
<section id="register"><div class="section-heading"><h2>Selection is not qualification.</h2><span>Engineering registry</span></div><p>The audit began with 35 identities and now includes a separate Victorville entry. Reviewed LCD/CRT and clear PETG appearances and protected organic anchors remain recorded, but their exact selected binary images/scenes are not all present here. Concrete formation, newer carpet constructions and organic research remain experimental. No newer experiment is promoted by this gallery.</p><details><summary>Inspect all {len(rows)} registry entries and availability</summary><div class="table-wrap"><table><caption>Frozen registry: selection, quality and availability are separate</caption><thead><tr><th>Material</th><th>Presentation status</th><th>Recorded limits</th></tr></thead><tbody>{''.join(rows)}</tbody></table></div></details><p>{link('audit/STATUS.md','Full audit status')} · {link(REGISTRY,'Exact selected registry')} · {link('audit/ARTIFACTS.json','Archived-input recovery inventory')}</p></section>
<section id="evidence"><div class="section-heading"><h2>Evidence stays intact.</h2><span>Source / scale / settings</span></div><div class="evidence-grid"><article><h3>Original images</h3><p>Native image links serve unchanged source bytes. Preview WebPs are downsampled display derivatives. A 4K map is never described as a 4K beauty render.</p><a href="catalog.json">Per-image SHA-256 &amp; dimensions ↗</a></article><article><h3>Editable construction</h3><p>Saved Blender scenes retain cameras, lighting, geometry and relative texture paths. Scene replay requires the materials directory alongside the Blender directory.</p>{link('blender/CYBR_Cycles_Details.blend','Specimen scenes')} · {link('blender/CYBR_Cycles_Still_Lifes.blend','Still lifes')} · {link('blender/CYBR_Cycles_Architecture.blend','Architecture')}</article><article><h3>Review &amp; reproduction</h3><p>Rebuilt studies, engineering selections and their accepted baselines remain distinct.</p>{link('docs/STUDY_REPAIR_REVIEW.md','Study review')} · {link('docs/IMPORT_GUIDE.md','PBR import guide')} · {link('docs/DOWNLOAD_RECEIPT_2026-10-12.json','Distribution receipt')} · {link('SHA256SUMS.txt','Original suite checksums')}</article></div></section>'''
    (OUT/'index.html').write_text(page('Surface collection',body))
    (OUT/'catalog.json').write_text(json.dumps(catalog,indent=2)+'\n')
    (ROOT/'docs/MATERIAL_GALLERY.md').write_text('# Individual material panels\n\n[Browse the static gallery](../gallery/index.html). Every channel opens independently; the complete suite ZIP is optional.\n\nNative maps are **4096 × 4096**; the existing specimen and matched-proof renders are **1024 × 1024**. Image previews are downsampled display copies.\n\n| Material | Panel | Native texture | Actual rendered appearance |\n|---|---|---|---|\n'+'\n'.join(direct_rows)+'\n\nThe gallery catalog records original image dimensions, bit depths, byte counts and SHA-256 checksums. Registry-based additions retain their separate review scope and archive dependencies.\n\n## Every native channel, directly\n\n'+'\n'.join(direct_details))
    print(f'Built {len(cards)} material pages, {len(maps)*len(cards)} map links, {len(selected_specs)} reviewed additions and {len(scene_cards)} context views.')


if __name__ == '__main__':
    build()
