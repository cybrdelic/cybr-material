"""Eight-variable native-curvature refit. Every other fibre/key is fixed.

Sampled residuals steer one bounded diagnostic trial. Authoritative acceptance
is separate and requires continuous native curvature plus all geometry gates.
"""
import json
import time
import numpy as np
from scipy.optimize import least_squares
from stage_keys import ROOT,UNIT,WORK,BridgeField,native,digest,changed_spans,sample_screen

def run():
    start=time.monotonic();frozen=UNIT/'revisions/r9_physical_resume/arrays/loop_968_partial.npz'
    assert digest(frozen)=='ab087996ecf9c0bc26cff51520e4ab3ca26b40ff7cc19bae228d3de163839cee'
    state=np.load(frozen);field=BridgeField(968,np.load(WORK/'carpet-r5-packed-core/receipts/baseline_loops.npz'))
    old=state['world'].astype('f4').transpose(1,0,2);s=state['stations'];refs=state['references'][35:39,0]/1e-6
    initial=state['xy'][35:39,0]/1e-6;radius=float(np.float32(field.r));t=np.linspace(0,1,129)
    basis=[field.basis(v) for v in s[35:39]]
    def points(z):
        xy=z.reshape(4,2)*1e-6;p=old[0].astype(float).copy()
        for k,(c,e1,e2,_) in enumerate(basis):p[35+k]=c+xy[k,0]*e1+xy[k,1]*e2
        return p
    def ratios(p):
        B=native.native_beziers(p)[33:40];D=3*np.diff(B,axis=1);E=2*np.diff(D,axis=1)
        V=(1-t)[None,:,None]**2*D[:,0,None]+2*(t*(1-t))[None,:,None]*D[:,1,None]+t[None,:,None]**2*D[:,2,None]
        A=(1-t)[None,:,None]*E[:,0,None]+t[None,:,None]*E[:,1,None]
        speed=np.linalg.norm(V,axis=2)
        return np.linalg.norm(np.cross(V,A),axis=2)/np.maximum(speed,1e-30)**3*radius,speed
    best=[float('inf'),initial.ravel().copy()];calls=0
    class Stop(Exception):pass
    def objective(z):
        nonlocal calls
        calls+=1
        if calls>600 or time.monotonic()-start>25:raise Stop()
        ratio,speed=ratios(points(z));movement=np.linalg.norm(z.reshape(4,2)-refs,axis=1)
        raw=np.r_[np.maximum(ratio.ravel()-.75,0),20*np.maximum(movement-29.999,0),
                  .0001*(z-initial.ravel())]
        score=float(np.maximum(ratio.max()-1,0)+20*np.maximum(movement.max()-29.999,0))
        if score<best[0]:best[:]=[score,z.copy()]
        return raw
    before=ratios(old[0].astype(float))[0]
    initial=np.clip(initial,refs-29.998,refs+29.998)
    status='completed'
    try:
        result=least_squares(objective,initial.ravel(),bounds=((refs-29.999).ravel(),(refs+29.999).ravel()),
                             max_nfev=60,ftol=1e-10,xtol=1e-10,gtol=1e-10,x_scale='jac')
        status=result.message
    except Stop:status='bounded_solver_stop'
    xy=best[1].reshape(4,2);delta=xy-refs;length=np.linalg.norm(delta,axis=1)
    xy=refs+delta*np.minimum(1,29.999/np.maximum(length,1e-30))[:,None]
    body=old.copy();body[0]=points(xy.ravel()).astype('f4')
    ids=np.arange(49);B,O,changed,mapping=changed_spans(body,old,ids)
    screen=sample_screen(B[0:1,changed],radius)
    config=native.GateConfig(wall_seconds=15,max_leaves=100000,memory_limit_mib=900)
    budget=native.Budget(config);pad=256*native.EPS*max(1e-6,float(np.abs(body[0]).max()))*19
    try:
        proof,_=native.owner_curvature(B[0,changed],radius,config,budget,pad,0)
        local={'pass':True,**proof}
    except native.Unproven as exc:local={'pass':False,'code':exc.code,**exc.details}
    trial=np.load(UNIT/'revisions/r10_physical_length/arrays/trial_968_prefix_native.npz')
    path=ROOT/'arrays/native_key_candidate_968.npz'
    np.savez_compressed(path,body=body,wraps=trial['wraps'],body_radius=np.float32(field.r),wrap_radius=trial['wrap_radius'],
                        stations=s,key_ids=ids,references=state['references'],changed_spans=changed,
                        unchanged_original_span_mapping=mapping,original_body=old,changed_owners=np.array([0]),
                        allowed_key_indices=np.array([35,36,37,38]))
    report={'scope':'Owner0 only, existing keys35..38 only; eight transverse variables; no contact solver expansion.',
            'status':str(status),'optimizer_function_calls':calls,'elapsed_seconds':time.monotonic()-start,
            'original_affected_sampled_maximum_kappa_radius':float(before.max()),'actual_float32_screen':screen,
            'continuous_local_curvature':local,'candidate_sha256':digest(path),
            'changed_native_spans':np.flatnonzero(changed).tolist(),'other186_fibres_bitwise_preserved':np.array_equal(body[1:],old[1:]),
            'owner0_outside_keys_bitwise_preserved':np.array_equal(body[0,np.r_[0:35,39:49]],old[0,np.r_[0:35,39:49]]),
            'source_frozen_sha256':digest(frozen),'sampled_pass_is_not_qualification':True,
            'full_native_prefix_qualified':False,'accepted':False,'selected_r5_promoted':False}
    (ROOT/'receipts/owner0_refit.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(report),flush=True)

if __name__=='__main__':run()
