"""One unchanged original 43.63 micrometre span from qualified section 0."""
from pathlib import Path
import json,time,resource,hashlib
import numpy as np
from swept_feasibility import ROOT,R2,gradient_checks,SweptConstraints,solve
from advance_probe import Field,UNIT,check_section

def run():
    start=time.time();gradient=gradient_checks();print(json.dumps({'segment_gradient_validation':gradient}),flush=True)
    data=np.load(UNIT.parent/'carpet-r5-packed-core/receipts/baseline_loops.npz');field=Field(968,data);trace=np.load(UNIT/'arrays/constructor_trace_968.npz');initial=np.load(R2/'arrays/section_0_candidate.npy')
    ok,initial_check,a=check_section(field,0,initial,trace['xy'][0]);assert ok
    end=field.length/96;reference=initial+field.target(end)-field.target(0);c,e1,e2,k=field.basis(end);basis=np.stack((e1,e2),axis=1)
    wrap=np.stack([field.wraps(s) for s in np.linspace(0,end,5)]);oa=wrap[:-1].reshape(-1,3);ob=wrap[1:].reshape(-1,3)
    system=SweptConstraints(reference/field.R,k*field.R,field.rate*field.R,.045,16e-6/field.R,2e-6/field.R,30e-6/field.R,a/field.R,c/field.R,basis,oa/field.R,ob/field.R,np.full(len(oa),16e-6/field.R))
    # Directional check of the actual full endpoint residual Jacobian, away
    # from hinge boundaries. This includes body and fixed-wrap swept terms.
    x=reference.ravel()/field.R;values,jac=system.residual_jacobian(x);direction=np.random.default_rng(968).normal(size=len(x));direction/=np.linalg.norm(direction);h=1e-6
    plus=system.residual_jacobian(x+h*direction,False);minus=system.residual_jacobian(x-h*direction,False);numeric=(plus-minus)/(2*h);analytic=jac@direction;active=(values>1e-4)&(plus>0)&(minus>0);error=float(abs(numeric[active]-analytic[active]).max());assert error<1e-5
    gradient['actual_combined_active_rows_directional_error']=error;gradient['actual_combined_rows_tested']=int(active.sum());print(json.dumps({'combined_gradient_error':error,'rows_checked':int(active.sum())}),flush=True)
    candidate,report=solve(system,reference/field.R,.05e-6/field.R,seconds=30,max_evaluations=160)
    xy=candidate*field.R;ok,section_check,b=check_section(field,end,xy,reference);recovered=(b-c)@basis/field.R
    residual=system.residual_jacobian(recovered.ravel(),False);swept=system.assess_sweeps(recovered.ravel());accepted=ok and residual.max()*field.R<=.05e-6
    result={'scope':'One direct original near-root span, loop 968, same guide/seed/morphology and qualified initial section. No density, amplitude or count variant.','initial_section_check':initial_check,'segment_gradient_validation':gradient,'span_arclength_m':end,'solve':report,'float32_section_check':section_check,'float32_swept_check_normalized':swept,'float32_body_swept_margin_m':swept['body_swept_gap_lower_bound_normalized']*field.R,'float32_wrap_swept_margin_m':swept['wrap_swept_gap_lower_bound_normalized']*field.R,'float32_maximum_ordinary_distance_residual_m':float(residual.max()*field.R),'maximum_control_movement_m':float(np.linalg.norm(xy-reference,axis=1).max()),'physical_nonpenetration':bool(swept['body_swept_gap_lower_bound_normalized']>=0 and swept['wrap_swept_gap_lower_bound_normalized']>=0),'construction_margin_pass':bool(accepted),'all_swept_pairs_complete':True,'initial_section_changed':False,'native_curve_qualified':False,'seconds':time.time()-start,'peak_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'source').glob('*.py')}}
    np.savez_compressed(ROOT/'arrays/span_968_0.npz',initial_xy=initial,final_xy=xy,start_world=a,end_world=b,wraps=wrap,span_arclength=end)
    (ROOT/'receipts/swept_span_probe.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)

if __name__=='__main__':run()
