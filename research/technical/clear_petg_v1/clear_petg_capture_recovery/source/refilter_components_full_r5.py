"""Isolated retained-pass experiment; no render and no live pipeline mutation.
Prepare exports only positive T and unchanged guides. Operator owns OIDN.
Import combines filtered positive T/R with unchanged signed T terms and delta64.
"""
import bpy,sys,json,hashlib,numpy as np
from pathlib import Path
HERE=Path(__file__).resolve().parent;R=HERE.parent;P=R.parents[1];RENDER=P.parent/'cloud-eevee';sys.path.insert(0,str(RENDER));sys.path.insert(0,str(HERE))
assert hashlib.sha256((RENDER/'component_image_bridge.py').read_bytes()).hexdigest()=='540b9b39df71576661a58bbba81ed3818fe1317f4224ed570ef03241b1cd6b5b','Unreviewed shared bridge revision'
import component_image_bridge as bridge
from png_precision_check import compare_opaque_pngs
from oidn_contract import validate_stats
from two_component_signed_math_r2 import combine,partition,tests
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
mode,runarg,outarg=sys.argv[sys.argv.index('--')+1:];assert mode in {'prepare','import'};run=Path(runarg);out=Path(outarg);oldout=run/'INTERNAL_PANEL_guides/components';old=json.loads((oldout/'export_receipt.json').read_text());native=json.loads((run/'INTERNAL_PANEL.png.json').read_text());source=Path(native['job']['source']);assert sha(source)==old['source_sha256'];assert old['identity_exact_float32']
for k,v in old['pass_hashes'].items():assert sha(old['input_paths'][k])==v
account=oldout/'linear_component_accounting.npz';assert sha(account)==old['component_state_sha256'];state=np.load(account);t=state['transmission_emission'];r=state['physical_remainder'];delta=state['roundoff_residual'];positive_t,signed_t=partition(t,r,delta);assert np.any(positive_t>0)
bpy.ops.wm.open_mainfile(filepath=str(source));s=bpy.context.scene
for k,v in native['display_settings'].items():setattr(s.view_settings,k,v)
s.render.image_settings.file_format='PNG';s.render.image_settings.color_depth='16';s.render.image_settings.color_mode='RGB';s.render.dither_intensity=0
helpers={str(p):sha(p) for p in [Path(__file__),HERE/'two_component_math.py',HERE/'two_component_signed_math_r2.py',RENDER/'component_image_bridge.py',RENDER/'png_precision_check.py',RENDER/'oidn_contract.py']}
filtered_r_path=oldout/'denoised_remainder.pfm';original_result=json.loads((oldout/'component_receipt.json').read_text());assert sha(filtered_r_path)==original_result['denoised_remainder_sha256'];filtered_r=bridge.read_pfm(filtered_r_path);assert np.all(filtered_r>=0)
def layer_energy(fields):
 cam=s.objects[native['job']['camera']];h,w=t.shape[:2];q=cam.matrix_world.to_3x3();right=np.array(q.col[0]);up=np.array(q.col[1]);direction=-np.array(q.col[2]);origin=np.array(cam.location)
 xx=((np.arange(w)+.5)/w-.5)*cam.data.ortho_scale;yy=((np.arange(h)+.5)/h-.5)*cam.data.ortho_scale;X,Y=np.meshgrid(xx,yy);distance=(.0012-origin[2]-X*right[2]-Y*up[2])/direction[2];world_y=origin[1]+X*right[1]+Y*up[1]+distance*direction[1];world_x=origin[0]+X*right[0]+Y*up[0]+distance*direction[0];phase=2*np.pi*world_y/.0002
 regions=[]
 for xc in [-.010,-.003,.004,.010]:
  mask=(abs(world_x-xc)<=.0005)&(abs(world_y)<=.010);assert mask.sum()>1000
  A=np.stack([np.ones(mask.sum()),X[mask],Y[mask],np.cos(phase[mask]),np.sin(phase[mask])],1);rows={}
  for name,a in fields.items():
   lum=a[mask].astype('f8')@np.array([.2126,.7152,.0722]);beta=np.linalg.lstsq(A,lum,rcond=None)[0];rows[name]={'layer_frequency_amplitude':float(np.linalg.norm(beta[-2:])),'phase_coefficients':beta[-2:].tolist(),'mean_linear_luminance':float(lum.mean()),'fit_residual_rms':float(np.sqrt(np.mean((A@beta-lum)**2)))}
  regions.append({'center_x_m':xc,'width_m':.001,'y_interval_m':[-.010,.010],'pixels':int(mask.sum()),'fields':rows})
 return {'physical_pitch_m':.0002,'method':'Four preregistered interior physical strips; linear-light5,000cycles/m demodulation with planar illumination terms. Excludes cut edges and start-stop seam, which require separate native visual review.','regions':regions}
if mode=='prepare':
 assert not out.exists(),'Fresh output directory required';out.mkdir(parents=True)
 beauty=bridge.read_exr(old['input_paths']['beauty']);identity=combine(positive_t,r,delta,signed_t);assert np.array_equal(identity,beauty)
 bridge.save(identity,out/'identity_raw_passthrough.png',s);roundtrip=compare_opaque_pngs(run/'INTERNAL_PANEL.png',out/'identity_raw_passthrough.png');assert roundtrip['bit_depth']==16 and roundtrip['max_error_native']<=1
 albedo=bridge.read_exr(old['input_paths']['albedo']);normal=bridge.read_exr(old['input_paths']['normal']);features={'beauty':bridge.stats(positive_t),'albedo':bridge.stats(albedo),'normal':bridge.stats(normal)};contract=validate_stats(features,require_nonzero_features=True);assert contract['valid'],contract
 for name,a in [('transmission',positive_t),('albedo',albedo),('normal',normal)]:bridge.write_pfm(out/(name+'.pfm'),a)
 bridge.save(filtered_r,out/'REVIEW_existing_filtered_remainder.png',s)
 receipt={'mode':'both-physical-components-filtered','input_run':str(run),'source':str(source),'source_sha256':sha(source),'helpers':helpers,'component_state_path':str(account),'component_state_sha256':sha(account),'pass_hashes':old['pass_hashes'],'pass_paths':old['input_paths'],'existing_filtered_remainder':str(filtered_r_path),'existing_filtered_remainder_sha256':sha(filtered_r_path),'filter_inputs':{n:sha(out/(n+'.pfm')) for n in ['transmission','albedo','normal']},'identity_float32_exact':True,'native16_identity_roundtrip':roundtrip,'feature_contract':contract,'unit_tests':tests(),'retained_transmission_claim':False,'delta64_untouched':True,'signed_transmission_terms_untouched':True,'clipping':False,'noise_or_detail_acceptance':False,'scope':'Filter positive transmission/emission separately. Reuse the previously filtered diffuse/glossy/environment/volume remainder. No source material or radiance rescaling.'}
 negative=signed_t<0;scale=np.maximum(np.maximum(abs(identity.astype('f8')),abs(t.astype('f8'))+r),np.finfo(np.float32).tiny)
 receipt['signed_transmission_accounting']={'negative_channel_count':int(negative.sum()),'negative_pixel_count':int(np.any(negative,axis=2).sum()),'min':float(signed_t.min()),'sum':float(signed_t.astype('f8').sum()),'relative_max':float(np.max(-signed_t.astype('f8')/scale)),'relative_guard':2.**-16,'original_negative_terms_preserved_exactly':True,'signed_term_bytes_sha256':hashlib.sha256(signed_t.tobytes()).hexdigest(),'partition_rule':'Tpositive=T where T>=0 else 0; Tsigned=T where T<0 else 0. Final=filteredTpositive + unchangedTsigned + filteredR + unchangedDelta64; no original term discarded.'}
 receipt['existing_filtered_remainder_format']='Hash-verified PFM, not a saved-EXR proxy'
 receipt['layer_energy']=layer_energy({'raw_transmission':t,'raw_remainder':r,'filtered_remainder':filtered_r})
 (out/'prepare_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print('TRANSMISSION_REPROCESS_PREPARED',str(out),flush=True)
else:
 receipt=json.loads((out/'prepare_receipt.json').read_text());assert receipt['helpers']==helpers and receipt['source_sha256']==sha(source) and receipt['component_state_sha256']==sha(account) and receipt['existing_filtered_remainder_sha256']==sha(filtered_r_path)
 for n,h in receipt['filter_inputs'].items():assert sha(out/(n+'.pfm'))==h
 filtered_t=bridge.read_pfm(out/'denoised_transmission.pfm');final=combine(filtered_t,filtered_r,delta,signed_t);dest=out/'TWO_COMPONENT_FILTERED_PANEL.png';assert not dest.exists()
 bridge.save(filtered_t,out/'REVIEW_filtered_transmission.png',s,out/'filtered_transmission_linear.exr');bridge.save(final,dest,s,out/'two_component_filtered_linear.exr')
 receipt.update(output=str(dest),output_sha256=sha(dest),output_linear_sha256=sha(out/'two_component_filtered_linear.exr'),filtered_transmission_pfm_sha256=sha(out/'denoised_transmission.pfm'),filtered_transmission_stats=bridge.stats(filtered_t),final_stats=bridge.stats(final),both_components_filtered=True,transmission_emission_untouched=False,selection_promoted=False)
 receipt['final_layer_energy']=layer_energy({'raw_combined':combine(positive_t,r,delta,signed_t),'filtered_transmission':filtered_t,'two_component_filtered':final})
 (out/'result_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print('TWO_COMPONENT_FILTERED_PROBE_READY',str(dest),flush=True)
assert sha(source)==old['source_sha256']
