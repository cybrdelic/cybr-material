"""Reconstruct the fixed native tool and audited production state; no new material law."""
from pathlib import Path
import numpy as np,json,hashlib,time,resource,gc
from scipy.ndimage import gaussian_filter
from PIL import Image
from thermovisco_contact_audited import Parameters,simulate,asdict
from replica_geometry import sample_native
R=Path(__file__).resolve().parents[1];start=time.monotonic();g=json.loads((R/'inputs/mechanism_gates.json').read_text());b=json.loads((R/'inputs/production_native_bake.json').read_text());ledger=json.loads((R/'inputs/plastic_production_ledger.json').read_text());N=4096;WIDTH=.08;PIX=WIDTH/N
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
rng=np.random.default_rng(21052027);tool=np.zeros((N,N),np.float32)
for k in range(95000):
 cx,cy=rng.uniform(-.0003,WIDTH+.0003,2);radius=rng.uniform(.00012,.000245);depth=rng.uniform(.000020,.000044)
 ix0=max(0,int((cx-radius)/PIX));ix1=min(N,int((cx+radius)/PIX)+2);iy0=max(0,int((cy-radius)/PIX));iy1=min(N,int((cy+radius)/PIX)+2)
 if ix0>=ix1 or iy0>=iy1:continue
 x=((np.arange(ix0,ix1,dtype=np.float32)+.5)*PIX-cx)/radius;y=((np.arange(iy0,iy1,dtype=np.float32)+.5)*PIX-cy)/radius;rr=x[None,:]**2+y[:,None]**2;cap=depth*np.maximum(0,1-rr)**1.75;patch=tool[iy0:iy1,ix0:ix1];np.maximum(patch,cap,out=patch)
tool=gaussian_filter(tool,1.05,mode='reflect');normalized=tool/float(tool.max());tool_path=R/'state/native_tool_state4096.npy';np.save(tool_path,normalized);assert sha(tool_path)==g['input_tool_sha256'],'Native tool differs from preserved source';del tool;gc.collect();print('EXACT_TOOL_RECOVERED',sha(tool_path),flush=True)
p=Parameters();assert asdict(p)==g['parameters_all_assumed'];depths=np.linspace(0,g['source_tool_metric_max_m'],257);state=simulate(depths,p,False);legacy={k:state[k] for k in ledger['unchanged_saved_arrays_bit_exact']};state_path=R/'state/production_column_state.npz';np.savez_compressed(state_path,**legacy);assert sha(state_path)==b['column_state_sha256'],'Production arrays/archive differ from preserved state';np.savez_compressed(R/'state/audited_production_column_state.npz',**{k:v for k,v in state.items() if isinstance(v,np.ndarray)})
assert state['diagnostics']['max_full_cycle_mechanical_residual_J_m2']<1e-8 and state['diagnostics']['max_full_cycle_combined_residual_J_m2']<1e-5;print('EXACT_PRODUCTION_STATE_RECOVERED',sha(state_path),flush=True)
height=np.empty(normalized.shape,'f4');total=0.
for y in range(0,N,128):
 h=np.interp(np.asarray(normalized[y:y+128],dtype='f8')*g['source_tool_metric_max_m'],depths,state['retained_displacement_m']).astype('f4');height[y:y+128]=h;total+=float(h.sum(dtype='f8'))
mean=total/N/N;assert mean==b['mean_removed_retained_displacement_m'];height-=mean;del normalized;gc.collect()
rough_u16=np.empty((N,N),'u2');sigma=4.;halo=40
for y in range(0,N,128):
 a=max(0,y-halo);e=min(N,y+128+halo);h=height[a:e];res=h-gaussian_filter(h,sigma,mode='reflect');gy,gx=np.gradient(res,PIX);variance=gaussian_filter((gx*gx+gy*gy)*.5,sigma,mode='reflect');rough=np.power(.475**4+np.maximum(variance,0.),.25);block=np.clip(rough[y-a:min(y+128,N)-a],0,1);rough_u16[y:min(y+128,N)]=np.rint(block*65535).astype('u2')
rough_path=R/'maps/production_state_roughness_native4096.png';Image.fromarray(np.flipud(rough_u16)).save(rough_path);assert sha(rough_path)==b['sha256']['production_state_roughness_native4096.png'];del rough_u16,h,res,gy,gx,variance,rough,block;gc.collect()
enc=np.empty((N,N),'u2');maxerr=0.
for y in range(0,N,64):
 h=height[y:y+64].astype('f8');q=np.rint((.5+h/80e-6)*65535);assert q.min()>=0 and q.max()<=65535;enc[y:y+64]=q.astype('u2');maxerr=max(maxerr,float(abs((q/65535-.5)*80e-6-h).max()))
height_path=R/'maps/production_retained_height_native4096.png';Image.fromarray(np.flipud(enc)).save(height_path,compress_level=6);assert sha(height_path)==json.loads((R/'inputs/retained_height_audited.json').read_text())['sha256'];del enc;gc.collect()
# This is exactly the flat top's draft-scaled geometry grid in the preserved builder.
scale=1-np.tan(np.deg2rad(.8))*(.008-.004)/.04;axis=np.linspace(-.04,.04,513)*scale;lo,hi=axis[192],axis[320];fine_axis=np.linspace(lo,hi,1025);xx,yy=np.meshgrid(fine_axis,fine_axis);fine=sample_native(height,np.stack([xx.ravel(),yy.ravel()],axis=1)).reshape(1025,1025);coarse=fine[::8,::8].copy();np.savez_compressed(R/'state/interior_geometry_fields.npz',fine_axis_m=fine_axis,coarse_axis_m=axis[192:321],fine_height_m=fine,coarse_height_m=coarse);assert np.max(abs(fine_axis[::8]-axis[192:321]))<1e-17
report={'source_tool_sha256':sha(tool_path),'tool_exact_archive_bytes':True,'production_state_sha256':sha(state_path),'production_state_exact_archive_bytes':True,'full_cycle_ledger':state['diagnostics'],'native_roughness_sha256':sha(rough_path),'native_height_sha256':sha(height_path),'both_native_maps_exact_archive_bytes':True,'native_resolution':[N,N],'height_range_m':[float(height.min()),float(height.max())],'height_encoding_max_error_m':maxerr,'crop_extent_m':float(hi-lo),'crop_index_range_original_geometry':[192,320],'coarse_pitch_m':float(axis[1]-axis[0]),'fine_pitch_m':float(fine_axis[1]-fine_axis[0]),'coarse_geometry_grid':[129,129],'fine_geometry_grid':[1025,1025],'shared_samples_exact':True,'same_state_palette_and_roughness_for_first_control':True,'roughness_limitation':'Original Gaussian high-pass closure intentionally held for isolated geometry comparison; it is not yet the residual of the actual rendered basis. No full geometry/BRDF partition qualification.','model_limitations':g['limitations'],'elapsed_seconds':time.monotonic()-start,'peak_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024};(R/'receipts/native_recovery.json').write_text(json.dumps(report,indent=2));print('NATIVE_RECOVERY_READY',json.dumps({k:report[k] for k in ['crop_extent_m','coarse_pitch_m','fine_pitch_m','elapsed_seconds','peak_rss_mib']}),flush=True)
