from pathlib import Path
import json,hashlib,struct,zlib
import numpy as np
from PIL import Image
R=Path(__file__).resolve().parents[1];N=4096
(R/'tests').mkdir(exist_ok=True)
r=json.loads((R/'native4096/receipt.json').read_text())
assert r['native_resolution']==[N,N] and r['upscaled'] is False
assert r['fracture_model']['coarse_primitive_sha256']=='c4e263c446b1d626eede4f3b292868b764c300b23490481c50b7b7cb9b7b02c3'
assert r['fracture_model']['coarse_mask_sha256']=='9128c472f196cc5d2140abaad3916ec2a7ee0ee24b3627a2b2b1cf4c9b17d39c'
assert r['construction']['clipped_height_pixels']==0
assert r['fracture_model']['fine_state']['coverage_pixels_without_owner']==0
assert r['fracture_model']['fine_state']['exposure_pixels_without_relief']==0
for key,name in [('BaseColor','BaseColor.png'),('Roughness','Roughness.png'),('Metallic','Metallic.png'),('Normal_Object','Normal_Object_RGB16.png')]:
 assert hashlib.sha256((R/'native4096'/name).read_bytes()).hexdigest()==r['files'][key]
 assert Image.open(R/'native4096'/name).size==(N,N)
assert np.max(np.asarray(Image.open(R/'native4096/Metallic.png')))==0
h=np.load(R/'native4096/HeightState.npy',mmap_mode='r');g=np.load(R/'native4096/GeometryHeight.npy')
assert h.shape==(N,N) and g.shape==(513,513)
assert np.isfinite(h).all() and np.isfinite(g).all() and h.min()>=0 and h.max()<=1
# Decode the saved custom-filter-0 RGB16 PNG directly, bypassing image libraries'
# optional 8-bit conversion. Inspect an independent stratified sample of slopes.
p=(R/'native4096/Normal_Object_RGB16.png').read_bytes();k=8;stream=[]
while k<len(p):
 length=struct.unpack('>I',p[k:k+4])[0];typ=p[k+4:k+8];body=p[k+8:k+8+length];k+=length+12
 if typ==b'IHDR':assert struct.unpack('>IIBBBBB',body)==(N,N,16,2,0,0,0)
 if typ==b'IDAT':stream.append(body)
raw=np.frombuffer(zlib.decompress(b''.join(stream)),dtype='u1').reshape(N,1+N*6)
assert np.max(raw[:,0])==0
normal=raw[:,1:].copy().view('>u2').reshape(N,N,3)
y,x=np.meshgrid(np.arange(0,N,13),np.arange(0,N,17),indexing='ij')
x=x.ravel();y=y.ravel();decoded=normal[y,x].astype('f8')/65535*2-1
step=.25/N
sx=(h[y,(x+1)%N].astype('f8')-h[y,(x-1)%N].astype('f8'))*.004/(2*step)
sy=(h[(y+1)%N,x].astype('f8')-h[(y-1)%N,x].astype('f8'))*.004/(2*step)
expected=np.stack((-sx,sy,np.ones_like(sx)),axis=-1);expected/=np.linalg.norm(expected,axis=-1,keepdims=True)
error=float(np.max(np.abs(expected-decoded)));unit=float(np.max(np.abs(np.linalg.norm(decoded,axis=-1)-1)))
assert error<=1/65535+1e-10 and unit<3e-5
report={'pass':True,'original_native4096':True,'file_hashes_verified':True,'exact_reviewed_coarse_primitive_and_mask_hashes':True,'height_encoding_unclipped':True,'metallic_zero':True,'finite_height_state':True,'full_object_normal_independent_sample_count':len(x),'normal_channel_error_max':error,'normal_unit_length_error_max':unit,'shader_aov_proof':'Separate render required; see aov/native/receipt.json when available','selected':False,'physical_qualification':False}
(R/'tests/native_verification_receipt.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
