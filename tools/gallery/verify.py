"""Verify gallery links, image dimensions, source hashes and loading contracts."""
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import urlsplit, unquote
import hashlib
import json

from PIL import Image

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'gallery'


class Document(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []
        self.images = []
        self.ids = set()

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if 'id' in a:
            assert a['id'] not in self.ids, f'Duplicate ID: {a["id"]}'
            self.ids.add(a['id'])
        if tag == 'a' and 'href' in a:
            self.links.append(a['href'])
        if tag == 'link' and 'href' in a:
            self.links.append(a['href'])
        if tag == 'script' and 'src' in a:
            self.links.append(a['src'])
        if tag == 'img':
            self.images.append(a)
            self.links.append(a['src'])
            for item in a.get('srcset', '').split(','):
                if item.strip():
                    self.links.append(item.strip().split()[0])


def main():
    documents = {}
    for p in OUT.rglob('*.html'):
        d = Document()
        d.feed(p.read_text())
        documents[p.resolve()] = d
        assert '<html lang="en">' in p.read_text()
        for im in d.images:
            assert im.get('alt') and im.get('width') and im.get('height'), (p, im)
            assert im.get('loading') in ('lazy', 'eager')
            assert im['src'].endswith('.webp'), f'Native image loaded as thumbnail: {im}'
            assert im.get('decoding') == 'async'
            image_path = (p.parent / im['src']).resolve()
            with Image.open(image_path) as preview:
                assert max(preview.size) <= 800
    for p, d in documents.items():
        for href in d.links:
            target = urlsplit(href)
            if target.scheme or target.netloc:
                continue
            dest = (p.parent / unquote(target.path)).resolve() if target.path else p
            assert dest.is_relative_to(ROOT), (p, href)
            assert dest.is_file(), f'Broken link: {p.relative_to(ROOT)} -> {href}'
            if target.fragment and dest in documents:
                assert target.fragment in documents[dest].ids, f'Missing anchor: {href}'

    catalog = json.loads((OUT / 'catalog.json').read_text())
    assert len(catalog['materials']) == 10
    assert len(catalog['selected']) >= 1
    assert len(catalog['scenes']) == 7
    assets = []
    for material in catalog['materials']:
        assert len(material['maps']) == 12
        assets += material['maps'] + [material['beauty'], material['proof']]
        assert [material['beauty']['width'], material['beauty']['height']] == material['render_settings']['resolution']
        assert [material['proof']['width'], material['proof']['height']] == material['proof_settings']['resolution']
        for a in material['maps']:
            assert [a['width'], a['height']] == [4096, 4096]
            if a['channel'].startswith('Height'):
                assert a['bit_depth'] == 16
    for material in catalog['selected']:
        assets.append(material['image'])
        if 'alternate_light' in material:
            assets.append(material['alternate_light'])
    assets += [s['image'] for s in catalog['scenes']]
    for a in assets:
        p = ROOT / a['path']
        raw = p.read_bytes()
        assert hashlib.sha256(raw).hexdigest() == a['sha256'], p
        assert len(raw) == a['bytes'], p
        with Image.open(p) as im:
            assert list(im.size) == [a['width'], a['height']], p
        for preview in a['previews']:
            assert (OUT / preview).is_file(), preview
    registry_path = ROOT / catalog['registry_path']
    assert hashlib.sha256(registry_path.read_bytes()).hexdigest() == catalog['registry_sha256']
    print(f'PASS: {len(documents)} pages, {len(assets)} original image hashes/dimensions, all local links and preview loading contracts.')


if __name__ == '__main__':
    main()
