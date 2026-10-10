"""Minimal disk-constrained deterministic refits of four existing native keys.

Only keys35..38 may change. This is a research candidate, never an export.
Sampled curvature chooses a diagnostic candidate; it cannot accept geometry.
"""
import json
import time
import numpy as np
from scipy.interpolate import PchipInterpolator, CubicHermiteSpline
from stage_keys import ROOT,UNIT,WORK,BridgeField,native,digest,changed_spans,sample_screen

def refit(state,field,method):
    old=state['world'].astype('f4').transpose(1,0,2)
    s=state['stations'];local=s[35:39];a,b=34,39
    t=(local-s[a])/(s[b]-s[a]);xy=state['xy'];refs=state['references'][35:39]
    if method=='linear_transverse':
        candidate=(1-t)[:,None,None]*xy[a]+t[:,None,None]*xy[b]
    else:
        if method=='world_chord':
            world=(1-t)[None,:,None]*old[:,a,None]+t[None,:,None]*old[:,b,None]
        elif method=='physical_tangent_hermite':
            d0=(old[:,35].astype(float)-old[:,33])/(s[35]-s[33])
            d1=(old[:,40].astype(float)-old[:,38])/(s[40]-s[38])
            spline=CubicHermiteSpline(s[[a,b]],old[:,[a,b]].transpose(1,0,2),np.stack((d0,d1)))
            world=spline(local).transpose(1,0,2)
        elif method=='pchip_transverse':
            candidate=PchipInterpolator(s[[33,34,39,40]],xy[[33,34,39,40]],axis=0)(local)
            world=None
        else:raise ValueError(method)
        if world is not None:
            candidate=[]
            for k,parameter in enumerate(local):
                center,e1,e2,_=field.basis(parameter)
                delta=world[:,k].astype(float)-center
                candidate.append(np.stack((delta@e1,delta@e2),axis=1))
            candidate=np.asarray(candidate)
    raw=np.linalg.norm(candidate-refs,axis=2)
    delta=candidate-refs
    candidate=refs+delta*np.minimum(1,(30e-6-1e-9)/np.maximum(raw,1e-30))[:,:,None]
    body=old.copy()
    for k,parameter in enumerate(local):body[:,35+k]=field.world(parameter,candidate[k]).astype('f4')
    return body,{'unclipped_max_reference_motion_m':float(raw.max()),'clipped_keys':int((raw>30e-6-1e-9).sum())}

def run():
    start=time.monotonic();frozen=UNIT/'revisions/r9_physical_resume/arrays/loop_968_partial.npz'
    assert digest(frozen)=='ab087996ecf9c0bc26cff51520e4ab3ca26b40ff7cc19bae228d3de163839cee'
    state=np.load(frozen);field=BridgeField(968,np.load(WORK/'carpet-r5-packed-core/receipts/baseline_loops.npz'))
    old=state['world'].astype('f4').transpose(1,0,2);ids=np.arange(49);rows=[];best=None
    for method in ('linear_transverse','world_chord','physical_tangent_hermite','pchip_transverse'):
        body,details=refit(state,field,method);B,O,changed,mapping=changed_spans(body,old,ids)
        screen=sample_screen(B[:,changed],float(np.float32(field.r)))
        row={'method':method,**details,'screen':screen,'changed_native_spans':np.flatnonzero(changed).tolist(),
             'changed_key_indices':[35,36,37,38],'outside_keys_preserved':np.array_equal(body[:,np.r_[0:35,39:49]],old[:,np.r_[0:35,39:49]])}
        rows.append(row);print(json.dumps(row),flush=True)
        score=screen['maximum_sampled_kappa_radius']
        if best is None or score<best[0]:best=(score,body,changed,mapping,row)
    score,body,changed,mapping,row=best
    path=ROOT/'arrays/native_key_candidate_968.npz'
    trial=np.load(UNIT/'revisions/r10_physical_length/arrays/trial_968_prefix_native.npz')
    np.savez_compressed(path,body=body,wraps=trial['wraps'],body_radius=np.float32(field.r),wrap_radius=trial['wrap_radius'],
                        stations=state['stations'],key_ids=ids,references=state['references'],changed_spans=changed,
                        unchanged_original_span_mapping=mapping,original_body=old)
    report={'scope':'Four local keys35..38 refitted separately from frozen48-span source; no solver expansion.',
            'candidates':rows,'diagnostic_candidate_method':row['method'],'diagnostic_candidate_sha256':digest(path),
            'selection_is_not_acceptance':True,'sampled_pass_is_not_qualification':True,
            'native_curve_qualified':False,'production_loop_complete':False,'selected_r5_promoted':False,
            'source_frozen_sha256':digest(frozen),'elapsed_seconds':time.monotonic()-start}
    (ROOT/'receipts/native_key_refit.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'candidate_method':row['method'],'screen_maximum_kappa_radius':score,'elapsed_seconds':report['elapsed_seconds']}),flush=True)

if __name__=='__main__':run()
