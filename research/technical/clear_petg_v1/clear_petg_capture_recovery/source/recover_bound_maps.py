"""Recover the selected bound bitmap fields at native resolution; no new material state."""
from pathlib import Path
import sys,json,hashlib,time,resource,gc
import numpy as np
from PIL import Image
R=Path(__file__).resolve().parents[1]
BASE=Path('/workspace/scratch/b4387906eb93')
SRC=BASE/'technical-recovery-20261006/restored/CYBR_technical_cloud_review/source'
PUB=BASE/'material-publication/audit/checkpoints/2026-10-05/sources/material-native4k'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(SRC/'expansion/maps.py')=='44d15a52e871088ea894f939c08e0bbe15b1d93a74fec19b5d11181b77270f2e'
assert sha(PUB/'source/generate_petg_bounded.py')=='70c9ff8c72e2fb78d703dce99d2354f1dde8976d7e40438537f6d4010e354bf1'
sys.path.insert(0,str(SRC));from surface_math import field,clamp,smooth,color
from expansion.catalog import get
s=get('23_petg_transparent');n=4096;old=json.loads((PUB/'sets/23_petg_transparent/material.json').read_text());t=time.monotonic();records=[]
def save(name,a):
 p=R/'maps'/(name+'.png');q=np.empty(a.shape,'u1')
 for j in range(0,n,64):q[j:j+64]=np.rint(a[j:j+64]*255).astype('u1')
 Image.fromarray(q).save(p,compress_level=4)
 expected=old['maps'][name]['sha256'];actual=sha(p);assert actual==expected,(name,actual,expected)
 records.append({'name':name,'path':str(p),'sha256':actual,'native_dimensions':[n,n],'byte_exact_archived':True,'color_space':old['maps'][name]['color_space']})
fine=field(n,64,s['seed']+1);broad=field(n,5,s['seed']+2,2)
rgb=clamp(color(s['color'],.005*broad)).astype('f4');save('BaseColor',rgb);del rgb,broad;gc.collect()
y=(np.arange(n,dtype='f4')[:,None]+.5)/n;phase=(y*s['tile_m']/s['layer_height_m'])%1;boundary=np.broadcast_to(1-smooth(np.minimum(phase,1-phase),.015,.12),(n,n)).copy();rough=clamp(s['roughness']+.14*boundary+.006*fine).astype('f4');save('Roughness',rough);del rough,boundary,fine;gc.collect()
save('Metallic',np.zeros((n,n),'f4'))
report={'source_generator':str(PUB/'source/generate_petg_bounded.py'),'source_generator_sha256':sha(PUB/'source/generate_petg_bounded.py'),'source_state':s,'native_resolution':[4096,4096],'maps':records,'scope':'Only existing shader-bound color, roughness and zero metallic; bead geometry remains explicit. No new height, bump or material coefficients.','wall_seconds':time.monotonic()-t,'peak_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024}
(R/'receipts/native_maps.json').write_text(json.dumps(report,indent=2)+'\n');print('CLEAR_NATIVE_FIELDS_RECOVERED',json.dumps(report),flush=True)
