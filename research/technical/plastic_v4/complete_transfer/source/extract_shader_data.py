import bpy,numpy as np,sys,json,hashlib
from pathlib import Path
args=sys.argv[sys.argv.index('--')+1:];assert len(args)==4;out=Path(args[3]);assert not out.exists();fields={}
for name,path in zip(['native_cell_height','q','normal_field'],args[:3]):
 im=bpy.data.images.load(path,check_existing=False);w,h=im.size;a=np.empty(w*h*4,'f4');im.pixels.foreach_get(a);bpy.data.images.remove(im);fields[name]=a.reshape(h,w,4)[:,:,:3].copy();assert np.isfinite(fields[name]).all()
np.savez_compressed(out,**fields);out.with_suffix('.json').write_text(json.dumps({'input_hashes':{p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in args[:3]},'npz_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'no_filter_or_color_transform':True},indent=2)+'\n')
