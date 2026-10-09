"""New porcelain glaze candidate: same qualified linear thin-film solver, source-calibrated finish."""
import numpy as np,json,hashlib,resource
from scipy.fft import rfft2,irfft2,rfftfreq,fftfreq
from PIL import Image
from pathlib import Path
R=Path(__file__).resolve().parents[1];N=4096;L=.0192
# Authored process parameters, not calibrated firing measurements.
P={'tile_m':L,'resolution':N,'pressed_powder_cell_m':L/160,'deposition_cell_m':L/80,'powder_skin_amplitude_m':4e-6,'wet_deposit_amplitude_m':6e-6,'nominal_glaze_thickness_m':60e-6,'integrated_capillary_mobility_m4':2e-19,'source_roughness':.17984575,'seed':15047}
rng=np.random.default_rng(P['seed']);particles=[]
def deposit(cells,amplitude):
 a=np.zeros((N,N),np.float32);step=N/cells
 for j in range(cells):
  for i in range(cells):
   cx=(i+.5+rng.uniform(-.16,.16))*step;cy=(j+.5+rng.uniform(-.16,.16))*step;rx=step*rng.uniform(.52,.73);ry=rx*rng.uniform(.85,1.15);height=amplitude*rng.uniform(.7,1.25)
   xx=np.arange(int(cx-rx*1.8),int(cx+rx*1.8+1));yy=np.arange(int(cy-ry*1.8),int(cy+ry*1.8+1));u=(xx[None,:]-cx)/rx;v=(yy[:,None]-cy)/ry;profile=np.exp(-2*(u*u+v*v)).astype('f4')*height;a[np.ix_(yy%N,xx%N)]+=profile;particles.append([cx/N*L,cy/N*L,rx/N*L,ry/N*L,height])
 return a
body=deposit(160,P['powder_skin_amplitude_m']);wet=deposit(80,P['wet_deposit_amplitude_m']);h0=body+wet+P['nominal_glaze_thickness_m'];del wet
kx=(2*np.pi*rfftfreq(N,L/N)).astype('f4');ky=(2*np.pi*fftfreq(N,L/N)).astype('f4');k4=(kx[None,:]**2+ky[:,None]**2)**2;F=rfft2(h0);H=irfft2(F*np.exp(-P['integrated_capillary_mobility_m4']*k4),s=(N,N)).astype('f4');Hhot=irfft2(F*np.exp(-2*P['integrated_capillary_mobility_m4']*k4),s=(N,N)).astype('f4');del F,k4
rms0=float(np.std(h0));rms1=float(np.std(H));rms2=float(np.std(Hhot));mean_error=float(abs(np.mean(H,dtype='f8')-np.mean(h0,dtype='f8')));assert rms2<rms1<rms0 and mean_error<1e-10;del Hhot,h0
body-=np.mean(body);H-=np.mean(H);body_scale=20e-6;glaze_scale=40e-6
assert np.max(abs(body))<body_scale*.5 and np.max(abs(H))<glaze_scale*.5
maps={}
def save(name,a,sixteen=False):
 assert np.isfinite(a).all() and a.min()>=0 and a.max()<=1
 out=np.empty(a.shape,'u2' if sixteen else 'u1')
 for j in range(0,N,64):out[j:j+64]=np.rint(a[j:j+64]*(65535 if sixteen else 255)).astype(out.dtype)
 p=R/'maps'/(name+'.png');Image.fromarray(out).save(p,compress_level=4);maps[name]={'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'resolution':[N,N],'bits':16 if sixteen else 8,'color_space':'Non-Color','range':[float(a.min()),float(a.max())]}
save('Substrate_Height',.5+body/body_scale,True);save('Glaze_Height',.5+H/glaze_scale,True)
dx=(np.roll(H,-1,1)-np.roll(H,1,1))*N/(2*L);dy=(np.roll(H,-1,0)-np.roll(H,1,0))*N/(2*L);slope=np.sqrt(dx*dx+dy*dy);rough=np.clip(.17984575+.008*(slope/(np.mean(slope)+1e-12)-1),.17184575,.19584575).astype('f4');save('Glaze_Roughness',rough)
enc=np.empty((N,N,3),'u1')
for j in range(0,N,64):
 normal=np.stack((-dx[j:j+64],dy[j:j+64],np.ones_like(dx[j:j+64])),axis=-1);normal/=np.linalg.norm(normal,axis=-1,keepdims=True);enc[j:j+64]=np.rint((normal*.5+.5)*255).astype('u1')
p=R/'maps/Glaze_Normal_OpenGL.png';Image.fromarray(enc).save(p,compress_level=4);maps['Glaze_Normal_OpenGL']={'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'resolution':[N,N],'bits':8,'color_space':'Non-Color','range':[float(enc.min()/255),float(enc.max()/255)]};del enc,normal,dx,dy,slope
# Particle positions/radii/heights are sufficient to reconstruct both native fields in order.
np.savez_compressed(R/'state/coating_state.npz',particles=np.array(particles,'f4'),pressed_count=np.int32(160*160),resolution=np.int32(N),tile_m=np.float32(L),integrated_mobility_m4=np.float64(P['integrated_capillary_mobility_m4']))
report={'parameters':P,'state':'Finite positive deposition kernels in SI coordinates; pressed skin is prescribed. Actual linear fourth-order capillary evolution is solved spectrally. No full firing, densification, shrinkage or viscosity-temperature solver.','shader_mapping':'Same intrinsic top color/coat values; derived surface normals and narrow empirical roughness response share the leveled glaze state. Matte substrate uses its pressed-skin field.','tests':{'mean_surface_height_conserved_m':mean_error,'initial_rms_m':rms0,'leveled_rms_m':rms1,'double_mobility_rms_m':rms2,'monotone_leveling':True,'finite_range_checks':True},'height_scale_m':{'Substrate_Height':body_scale,'Glaze_Height':glaze_scale},'maps':maps,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss};(R/'receipts/coating_state.json').write_text(json.dumps(report,indent=2));print('COATING_STATE_READY',report['tests'],report['peak_rss_kib'],flush=True)
