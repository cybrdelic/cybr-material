"""Verify actual native-cell, height, and analytic normal readbacks.
Geometry qualification is inherited only through frozen bitwise equality checks.
"""
from pathlib import Path
import sys,json,hashlib
import numpy as np
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'source'))
from surface_contract import HEIGHT,TANGENTS,normal,native_field,base_and_direction
from analytic_shader import NumericOps,field_expressions
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

def angles(x,y):
    x=x/np.linalg.norm(x,axis=-1,keepdims=True);y=y/np.linalg.norm(y,axis=-1,keepdims=True)
    return np.degrees(np.arctan2(np.linalg.norm(np.cross(x,y),axis=-1),np.sum(x*y,axis=-1)))
def stats(x):
    return {'max':float(x.max()),'rms':float(np.sqrt(np.mean(x*x))),'p99':float(np.quantile(x,.99))} if len(x) else {'max':None,'rms':None,'p99':None}

def verify(actual_path,meta):
    extraction=json.loads(actual_path.with_suffix('.json').read_text());assert sha(actual_path)==extraction['npz_sha256'];a=np.load(actual_path);mask=np.any(a['q']!=0,axis=-1);assert mask.sum()>500
    q_local=(a['q'][mask].astype('f8')-.5)/meta['q_encoding_scale'];q=q_local+meta['q_encoding_origin_m'];q_input=q_local.astype('f4')+np.array(meta['q_shader_origin_actual_f32'],'f4');measured=a['normal_field'][mask].astype('f8')*2-1;length=np.linalg.norm(measured,axis=-1);assert np.all(length>0);measured/=length[:,None]
    direct=a['native_cell_height'][mask].astype('f8');cell_float=direct[:,:2]*4096;observed_cell=np.rint(cell_float).astype('i4');assert np.all((observed_cell>=0)&(observed_cell<=4095));cell_integer_error=float(np.max(np.abs(cell_float-observed_cell)));assert cell_integer_error<.01;measured_h=(direct[:,2]-.5)*80e-6
    chart=np.where(abs(q[:,2]-.004)<=abs(q[:,0]-.04),0,2);exact_normal=normal(q,chart);base=base_and_direction(q)[0];exact_h,_=native_field(base[:,:2]);raw=(base[:,:2]/.08+.5)*4096-.5;exact_cell=np.floor(np.clip(raw,0,4095)).astype('i4')
    oracle_normal=np.empty_like(q);observed_owner_normal=np.empty_like(q);oracle_cell=np.empty_like(observed_cell);oracle_h=np.empty(len(q));observed_owner_h=np.empty(len(q))
    for c in [0,2]:
        sel=chart==c;op=NumericOps(HEIGHT,'f4');normal_call=field_expressions(op,list(q_input[sel].T),TANGENTS[c,0],TANGENTS[c,1]);owner_call=field_expressions(op,list(q_input[sel].T),TANGENTS[c,0],TANGENTS[c,1],native_cell_override=observed_cell[sel].astype('f4'))
        oracle_normal[sel]=np.stack(normal_call['normal'],-1);observed_owner_normal[sel]=np.stack(owner_call['normal'],-1);oracle_cell[sel]=np.stack(normal_call['cell'],-1);oracle_h[sel]=normal_call['height'];observed_owner_h[sel]=owner_call['height']
    actual_exact_swap=np.any(observed_cell!=exact_cell,axis=-1);actual_float_swap=np.any(observed_cell!=oracle_cell,axis=-1);oracle_exact_swap=np.any(oracle_cell!=exact_cell,axis=-1);distance=np.min(np.abs(raw-np.rint(raw)),axis=-1)*(.08/4096)
    actual_exact=angles(measured,exact_normal);implementation=angles(measured,oracle_normal);condition=angles(oracle_normal,exact_normal);owner_implementation=angles(measured,observed_owner_normal);same_cell=~actual_exact_swap
    inherited=meta['inherited_position_qualification'];assert inherited['mesh_transform_q_attribute_bitwise_equal'] and inherited['position_pass'];assert sha(inherited['receipt'])==inherited['sha256']
    return {'data_origin':'synthetic_verifier_control' if extraction.get('synthetic') else 'actual_Cycles_EXRs','actual_input_hashes':extraction['input_hashes'],'hit_pixels':int(mask.sum()),'inherited_geometry_qualification':inherited,'observed_native_cell_integer_readback_max_error':cell_integer_error,'actual_vs_exact_coordinate_cell_mismatch_count':int(actual_exact_swap.sum()),'actual_vs_float32_oracle_cell_mismatch_count':int(actual_float_swap.sum()),'float32_oracle_vs_exact_cell_mismatch_count':int(oracle_exact_swap.sum()),'actual_cell_mismatch_knot_distance_max_m':float(distance[actual_exact_swap].max()) if actual_exact_swap.any() else 0.,'normal_implementation_vs_float32_oracle_degrees':stats(implementation),'normal_conditioning_float32_vs_exact_degrees':stats(condition),'normal_actual_vs_exact_degrees':stats(actual_exact),'normal_actual_same_cell_vs_exact_degrees':stats(actual_exact[same_cell]),'normal_implementation_with_observed_cell_degrees':stats(owner_implementation),'height_actual_vs_float32_oracle_m':stats(np.abs(measured_h-oracle_h)),'height_actual_vs_observed_cell_oracle_m':stats(np.abs(measured_h-observed_owner_h)),'height_actual_vs_exact_coordinate_m':stats(np.abs(measured_h-exact_h)),'actual_normal_field_lengths':[float(length.min()),float(length.max())],'normal_measurement_scope':meta['normal_capture_semantics'],'source_sha256':meta['source_sha256'],'shader_source_sha256':meta['shader_source_sha256'],'new_acceptance_gate_created':False,'whole_part_or_appearance_qualified':False,'no_automatic_promotion':True,'oracle_limits':['Coordinate AOV decoding has finite precision; float32 oracle reproduces declared operation sequence, not undocumented compiler fusion.','Observed-cell oracle deliberately uses the renderer-reported neighboring cell at the reconstructed coordinate; any tiny extrapolation is diagnostic only and never changes the material.','All actual-versus-exact worst cases are reported, including derivative jumps.']}

if __name__=='__main__':
    actual,out=map(Path,sys.argv[1:]);meta=json.loads((R/'receipts/shader_fixture_build.json').read_text());assert meta['shader_source_sha256']==sha(R/'source/analytic_shader.py');result=verify(actual,meta);out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
