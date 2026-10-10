"""Bounded insertion-only native correction; frozen keys remain byte-identical."""
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
UNIT=ROOT.parents[1];WORK=UNIT.parent
sys.path.insert(0,str(WORK))
sys.path.insert(0,str(UNIT/'revisions/r6_c2_bridge/source'))
from bridge_field import BridgeField
sys.path.insert(0,str(UNIT/'native-qualification'))
import native_curve_gate as native

def digest(path):return hashlib.file_digest(Path(path).open('rb'),'sha256').hexdigest()

def split_at(bez,t):
    first=(1-t)*bez[:,:-1]+t*bez[:,1:]
    second=(1-t)*first[:,:-1]+t*first[:,1:]
    mid=(1-t)*second[:,0]+t*second[:,1]
    return np.stack((bez[:,0],first[:,0],second[:,0],mid),axis=1),np.stack((mid,second[:,1],first[:,2],bez[:,-1]),axis=1)

def restrict(bez,lo,hi):
    if hi<1:bez=split_at(bez,hi)[0]
    if lo>0:bez=split_at(bez,lo/hi)[1]
    return bez

def changed_spans(body,old,key_ids):
    B=np.stack([native.native_beziers(p) for p in body])
    O=np.stack([native.native_beziers(p) for p in old])
    changed=[];mapping=[]
    for k in range(body.shape[1]-1):
        i,j=key_ids[k:k+2]
        same=i>=0 and j==i+1 and np.array_equal(B[:,k],O[:,i])
        changed.append(not same);mapping.append(i if same else -1)
    return B,O,np.asarray(changed),np.asarray(mapping)

def generate(state,field,insertions):
    stations=state['stations'];old=state['world'].astype('f4').transpose(1,0,2)
    points=[];new_s=[];ids=[];refs=[]
    for k,s in enumerate(stations):
        points.append(old[:,k]);new_s.append(s);ids.append(k);refs.append(state['references'][k])
        for fraction in insertions.get(k,[]):
            t=s+fraction*(stations[k+1]-s)
            xy=(1-fraction)*state['xy'][k]+fraction*state['xy'][k+1]
            reference=(1-fraction)*state['references'][k]+fraction*state['references'][k+1]
            delta=xy-reference;length=np.linalg.norm(delta,axis=1)
            xy=reference+delta*np.minimum(1,(30e-6-1e-9)/np.maximum(length,1e-30))[:,None]
            points.append(field.world(t,xy).astype('f4'));new_s.append(t);ids.append(-1);refs.append(reference)
    return np.stack(points,axis=1),np.asarray(new_s),np.asarray(ids),np.asarray(refs)

def sample_screen(B,radius):
    t=np.linspace(0,1,129)
    D=3*np.diff(B,axis=2);E=2*np.diff(D,axis=2)
    velocity=(1-t)[None,None,:,None]**2*D[:,:,0,None]+2*(t*(1-t))[None,None,:,None]*D[:,:,1,None]+t[None,None,:,None]**2*D[:,:,2,None]
    acceleration=(1-t)[None,None,:,None]*E[:,:,0,None]+t[None,None,:,None]*E[:,:,1,None]
    speed=np.linalg.norm(velocity,axis=-1)
    ratio=np.linalg.norm(np.cross(velocity,acceleration),axis=-1)/np.maximum(speed,1e-30)**3*radius
    arg=np.unravel_index(np.argmax(ratio),ratio.shape)
    return {'maximum_sampled_kappa_radius':float(ratio[arg]),'owner':int(arg[0]),
            'changed_span_index':int(arg[1]),'parameter':float(t[arg[2]]),
            'minimum_sampled_speed_m':float(speed.min()),'sampled_pass_is_not_qualification':True}

def run():
    start=time.monotonic()
    (ROOT/'arrays').mkdir(exist_ok=True);(ROOT/'receipts').mkdir(exist_ok=True)
    frozen=UNIT/'revisions/r9_physical_resume/arrays/loop_968_partial.npz'
    assert digest(frozen)=='ab087996ecf9c0bc26cff51520e4ab3ca26b40ff7cc19bae228d3de163839cee'
    state=np.load(frozen);data=np.load(WORK/'carpet-r5-packed-core/receipts/baseline_loops.npz')
    field=BridgeField(968,data);radius=np.float32(field.r)
    old=state['world'].astype('f4').transpose(1,0,2)
    wraps=np.load(UNIT/'revisions/r10_physical_length/arrays/trial_968_prefix_native.npz')['wraps']
    options=[{35:[.5]},{35:[.75]},{35:[.625]},{35:[.875]},
             {35:[.5,.75]},{35:[.5],36:[.5]},
             {35:[.5,.75],36:[.5]},
             {35:[.5,.75,.875],36:[.5,.75]}]
    rows=[];chosen=None
    for index,insertions in enumerate(options):
        body,stations,ids,refs=generate(state,field,insertions)
        B,O,changed,mapping=changed_spans(body,old,ids)
        sample=sample_screen(B[:,changed],float(radius))
        row={'candidate':index,'insertions':insertions,'added_keys':body.shape[1]-old.shape[1],
             'changed_native_spans':np.flatnonzero(changed).tolist(),'screen':sample,
             'all_original_keys_bitwise_preserved':np.array_equal(body[:,ids>=0],old),
             'continuous_local_curvature_pass':False}
        if sample['maximum_sampled_kappa_radius']<.999:
            config=native.GateConfig(wall_seconds=20,max_leaves=100000)
            budget=native.Budget(config);maximum=0;failure=None
            for owner in range(187):
                pad=256*native.EPS*max(1e-6,float(np.abs(body[owner]).max()))*19
                try:
                    proof,_=native.owner_curvature(B[owner,changed],float(radius),config,budget,pad,owner)
                    maximum=max(maximum,proof['kappa_radius_upper'])
                except native.Unproven as exc:
                    failure={'code':exc.code,**exc.details};break
            row['continuous_local_curvature_pass']=failure is None
            row['maximum_continuous_local_kappa_radius_upper']=maximum
            row['continuous_failure']=failure
            if failure is None:chosen=(body,stations,ids,refs,changed,mapping,insertions,row)
        rows.append(row);print(json.dumps(row),flush=True)
        if chosen is not None:break
        if time.monotonic()-start>30:break
    report={'scope':'Insertion-only correction on frozen48-span loop968 prefix; no solver or promotion.',
            'source_frozen_sha256':digest(frozen),'candidates':rows,'selected_candidate':None,
            'original_keys_and_prefix_endpoints_preserved':True,'production_loop_complete':False,
            'native_curve_qualified':False,'elapsed_seconds':time.monotonic()-start}
    if chosen is not None:
        body,stations,ids,refs,changed,mapping,insertions,row=chosen
        path=ROOT/'arrays/native_key_candidate_968.npz'
        np.savez_compressed(path,body=body,wraps=wraps,body_radius=radius,wrap_radius=np.float32(16e-6),
                            stations=stations,key_ids=ids,references=refs,changed_spans=changed,
                            unchanged_original_span_mapping=mapping,original_body=old)
        report.update(selected_candidate=row['candidate'],candidate_sha256=digest(path),
                      continuous_local_curvature_pass=True,
                      maximum_continuous_local_kappa_radius_upper=row['maximum_continuous_local_kappa_radius_upper'],
                      new_reference_provenance='Linear interpolation of original transverse reference disks; only inserted keys have new references. Original keys/references remain exact.',
                      inserted_world_provenance='Qualified deterministic C2 guide plus linear interpolation of frozen transverse coordinates, with1nm quantization inset inside30um disks.')
    (ROOT/'receipts/native_key_stage.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':run()
