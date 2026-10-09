"""Report physical layer-frequency contrast; never an automatic realism gate."""
import bpy,sys,json,numpy as np
from pathlib import Path
from mathutils import Vector
from bpy_extras.object_utils import world_to_camera_view
R=Path(__file__).resolve().parents[1];run=Path(sys.argv[sys.argv.index('--')+1]);rec=json.loads((R/'receipts/tangent_repair.json').read_text());bpy.ops.wm.open_mainfile(filepath=rec['scene']);s=bpy.context.scene;cam=s.objects['Slab / shared camera']
def read(path):
 im=bpy.data.images.load(str(path),check_existing=False);w,h=im.size;a=np.empty(w*h*4,'f4');im.pixels.foreach_get(a);bpy.data.images.remove(im);return a.reshape(h,w,4)[...,:3].astype('f8')
raw=read(run/'INTERNAL_PANEL_guides/beauty_0001.exr');clean=read(run/'INTERNAL_PANEL_guides/denoised_linear.exr');assert raw.shape==clean.shape;h,w,_=raw.shape
def bilinear(a,x,y):
 i=np.floor(x).astype(int);j=np.floor(y).astype(int);fx=x-i;fy=y-j
 assert i.min()>=0 and j.min()>=0 and i.max()+1<w and j.max()+1<h
 return (a[j,i]*(1-fx)*(1-fy)+a[j,i+1]*fx*(1-fy)+a[j+1,i]*(1-fx)*fy+a[j+1,i+1]*fx*fy)
ys=np.linspace(-.010,.010,4096,endpoint=False);basis=np.stack([np.ones_like(ys),ys,np.cos(2*np.pi*ys/.0002),np.sin(2*np.pi*ys/.0002)],1);rows=[];weights=np.array([.2126,.7152,.0722]);lr=raw@weights;lc=clean@weights
for x in [-.010,-.003,.004,.010]:
 uv=np.array([world_to_camera_view(s,cam,Vector((x,float(y),.0012)))[:2] for y in ys]);px=uv[:,0]*w-.5;py=uv[:,1]*h-.5
 r=bilinear(lr,px,py);c=bilinear(lc,px,py);ar=np.linalg.lstsq(basis,r,rcond=None)[0];ac=np.linalg.lstsq(basis,c,rcond=None)[0];amp_r=float(np.linalg.norm(ar[2:]));amp_c=float(np.linalg.norm(ac[2:]));rows.append({'x_m':x,'raw_linear_layer_harmonic':amp_r,'clean_linear_layer_harmonic':amp_c,'amplitude_retention':amp_c/amp_r if amp_r else None,'raw_mean':float(r.mean()),'clean_mean':float(c.mean())})
report={'physical_layer_pitch_m':.0002,'world_y_interval_m':[-.010,.010],'sample_count_per_line':4096,'lines':rows,'image_size':[w,h],'method':'Linear-light luminance demodulated at the known5,000cycles/m layer frequency along projected nominal front-plane lines. Constant/linear illumination trends are fitted jointly.','limitations':'Finite projection of curved bead fronts, multiple light paths and image sampling can shift phase and harmonic energy. Report is a supporting diagnostic; raw native visual review remains necessary. It cannot establish optical clarity or a calibrated PETG model.','automatic_acceptance':False}
(R/'receipts/layer_retention.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
