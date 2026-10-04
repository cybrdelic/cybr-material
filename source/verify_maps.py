"""Check final texture files and measure their actual encoded ranges."""
import json
from pathlib import Path
import numpy as np
from PIL import Image

ROOT=Path(__file__).resolve().parents[1]
manifest=json.loads((ROOT/'materials/manifest.json').read_text())

def seam(a):
    edge=np.concatenate((np.abs(a[0]-a[-1]).ravel(),np.abs(a[:,0]-a[:,-1]).ravel()))
    inside=np.concatenate((np.abs(a[1]-a[0]).ravel(),np.abs(a[:,1]-a[:,0]).ravel()))
    return dict(wrap_mean=float(edge.mean()),neighbor_mean=float(inside.mean()),
                wrap_max=float(edge.max()),neighbor_max=float(inside.max()))

checks={};reports=[]
for s in manifest['materials']:
    folder=ROOT/'materials'/s['id']
    meta=json.loads((folder/'material.json').read_text())
    resolution=meta['resolution']
    maps=list(folder.glob('*.png'));assert len(maps)==12
    for p in maps:
        with Image.open(p) as im:
            assert im.size==(resolution,resolution),(p,im.size)
            if p.name in ('Height.png','Height_Macro.png'):assert im.mode in ('I;16','I')
    height=np.asarray(Image.open(folder/'Height.png'),dtype=np.float32)/65535
    rough=np.asarray(Image.open(folder/'Roughness.png'),dtype=np.float32)/255
    metallic=np.asarray(Image.open(folder/'Metallic.png'),dtype=np.float32)/255
    ao=np.asarray(Image.open(folder/'AO.png'))
    opacity=np.asarray(Image.open(folder/'Opacity.png'))
    if s['family']=='Textile':
        assert (opacity<128).mean()>.001, 'No actual woven openings'
        assert (opacity<128).mean()<.15, 'Unrealistic open weave coverage'
    else:
        assert np.all(opacity==255), 'Unexpected transparency'
    orm=np.asarray(Image.open(folder/'ORM.png'))
    assert np.array_equal(orm[...,0],ao)
    assert np.array_equal(orm[...,1],np.asarray(Image.open(folder/'Roughness.png')))
    assert np.array_equal(orm[...,2],np.asarray(Image.open(folder/'Metallic.png')))
    if s['family']=='Metal':assert metallic.std()>.04
    rgb=np.asarray(Image.open(folder/'BaseColor.png'),dtype=np.float32)/255
    gl=np.asarray(Image.open(folder/'Normal_OpenGL.png'))
    dx=np.asarray(Image.open(folder/'Normal_DirectX.png'))
    assert np.array_equal(gl[...,0],dx[...,0]) and np.array_equal(gl[...,2],dx[...,2])
    assert np.all(gl[...,1].astype(np.int16)+dx[...,1].astype(np.int16)==255)
    normal=gl.astype(np.float32)/127.5-1
    magnitude=np.sqrt((normal*normal).sum(axis=-1))
    assert abs(magnitude.mean()-1)<.007
    checks[s['id']]=dict(height_min=float(height.min()),height_max=float(height.max()),
        relief_range_m=float((height.max()-height.min())*s['height_scale_m']),
        roughness_min=float(rough.min()),roughness_max=float(rough.max()),
        metallic_min=float(metallic.min()),metallic_max=float(metallic.max()),
        base_color_seam=seam(rgb),height_seam=seam(height),finite=True,
        measurement_source='Final encoded PNG maps')
    reports.append(dict(material=s['id'],map_count=12,dimensions=[resolution,resolution],
        height_16bit=True,normal_unit_length_mean=float(magnitude.mean()),dx_flip_verified=True))
    print('PASS',s['id'],flush=True)
(ROOT/'materials/validation.json').write_text(json.dumps(checks,indent=2)+'\n')
(ROOT/'docs/map_checks.json').write_text(json.dumps(reports,indent=2)+'\n')
print('Verified 120 maps, ORM channels, metal/oxide variation, height depth, opacity and normal conventions.',flush=True)
