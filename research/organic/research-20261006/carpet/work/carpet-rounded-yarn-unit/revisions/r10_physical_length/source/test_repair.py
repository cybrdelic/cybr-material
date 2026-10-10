"""Meaningful window, complete coverage, derivative and Windows budget tests."""
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
from physical_length_window import (make_mask,pair_candidates,segment_closest,
                                   ROOT,UNIT,WORK,DISK_M)
from probe_repair import setup
from portable_resources import peak_rss_mib,accounting_receipt

def run():
    results={}
    stations=np.r_[np.arange(0.,10.),9.01,9.02,9.03]
    anchor,mask=make_mask(stations,[35,36],2.)
    assert stations[-1]-stations[anchor]>=2 and not mask[anchor].any()
    assert mask[-2:].all() and not mask[0].any()
    assert not mask[anchor+1:-2,np.r_[0:35,37:187]].any()
    other,_=make_mask(stations,[35,36],2.,physical_stations=stations*2)
    assert other>anchor, 'physical metric must affect selection when speed differs'
    results['physical_length_survives_subdivision_and_preserves_root']={'passed':True}
    rng=np.random.default_rng(81241)
    a=rng.uniform(-1,1,(300,3));b=a+rng.uniform(-.3,.3,(300,3))
    # A distant midpoint can still own an intersecting long segment.
    a[0]=[-1000.,0,0];b[0]=[.01,0,0];a[1]=[0,-.01,0];b[1]=[0,.01,0]
    affected=np.arange(300)<60;owner=np.arange(300)%187
    move=np.where(affected,.08,0.)
    i,j=pair_candidates(a,b,affected,owner,move,.03)
    got=set(zip(i.tolist(),j.tolist()));required=set()
    for first in np.flatnonzero(affected):
        second=np.arange(300);valid=(owner[first]!=owner[second])&(~affected[second]|(first<second))
        second=second[valid]
        d,_,_,_,e,_=segment_closest(np.broadcast_to(a[first],(len(second),3)),np.broadcast_to(b[first],(len(second),3)),a[second],b[second])
        for other in second[d-e<=.03+move[first]+move[second]]:required.add((int(first),int(other)))
    assert required<=got and (0,1) in got
    results['complete_disk_bounds_against_exhaustive_pairs_and_long_segment']={
        'passed':True,'exhaustive_required':len(required),'retained':len(got)}
    system,wraps,parameters,inputs=setup()
    assert system.anchor==34 and system.n==478 and len(system.stations)==50
    assert system.physical_window_width_m-system.arc_integration_error_m>=87.26166464481879e-6
    assert np.array_equal(system.refs*system.field.R,system.references[system.mask]) or np.max(abs(system.refs*system.field.R-system.references[system.mask]))<1e-18
    value,jac=system.evaluate(system.initial.ravel())
    worst=0.;rows_tested=0
    for _ in range(8):
        direction=rng.normal(size=system.n*2);direction/=np.linalg.norm(direction)
        # These adaptive spans can be less than1um long. A larger step changes
        # the closest segment feature; convergence was checked at five step sizes.
        h=2e-7
        plus=system.evaluate(system.initial.ravel()+h*direction,False)
        minus=system.evaluate(system.initial.ravel()-h*direction,False)
        stable=(value>1e-7)&(plus>1e-7)&(minus>1e-7)
        error=abs((plus-minus)/(2*h)-jac@direction)
        if float(error[stable].max(initial=0))>5e-5:
            ids=np.flatnonzero(stable)[np.argsort(error[stable])[-8:]]
            diagnostic=[]
            aa,bb=system.segments(system.initial.ravel())
            for row in ids:
                item={'row':int(row),'analytic':float((jac@direction)[row]),'finite':float(((plus-minus)/(2*h))[row]),'residual':float(value[row])}
                if system.n<=row<system.n+len(system.i):
                    k=row-system.n;i=system.i[k:k+1];j=system.j[k:k+1]
                    d,s,t,n,e,feature=segment_closest(aa[i],bb[i],aa[j],bb[j]);item.update(group='body',distance=float(d[0]),error=float(e[0]),s=float(s[0]),t=float(t[0]),feature=int(feature[0]))
                elif system.n+len(system.i)<=row<system.n+len(system.i)+len(system.wi):
                    k=row-system.n-len(system.i);i=system.wi[k:k+1];j=system.wj[k:k+1]
                    d,s,t,n,e,feature=segment_closest(aa[i],bb[i],system.wa[j],system.wb[j]);item.update(group='wrap',distance=float(d[0]),error=float(e[0]),s=float(s[0]),t=float(t[0]),feature=int(feature[0]))
                else:item['group']='disk_or_envelope'
                diagnostic.append(item)
            print(json.dumps({'gradient_diagnostic':diagnostic}),flush=True)
        worst=max(worst,float(error[stable].max(initial=0)))
        rows_tested+=int(stable.sum())
    assert rows_tested>100 and worst<5e-5,(rows_tested,worst)
    results['sparse_jacobian_finite_difference_on_real_window']={
        'passed':True,'directions':8,'stable_active_rows_tested':rows_tested,
        'normalized_finite_difference_step':h,'max_derivative_error':worst}
    actual,check=system.assess_actual(system.initial.ravel())
    assert np.array_equal(actual[~system.mask],system.original_world.astype('f4')[~system.mask])
    assert check['original30um_disks_pass']
    assert not check['native_curve_qualified']
    results['stored_float32_fixed_history_and_original_disk_identity']={'passed':True,
        'fixed_points':int((~system.mask).sum()),'max_movement_m':check['max_movement_from_original_reference_m'],
        'integrated_physical_window_width_m':system.physical_window_width_m,
        'arc_quadrature_error_estimate_m':system.arc_integration_error_m}
    # Cross-module collision regression: solver and gate must use their own kernels.
    import native_curve_gate as native
    assert native.projection_distance_lower.__module__=='_native_clearance_segments'
    before=peak_rss_mib();allocation=np.ones(8*1024**2,dtype=np.uint8);after=peak_rss_mib()
    del allocation
    assert after>=before and peak_rss_mib()>=after
    results['native_kernel_namespace_and_true_peak_memory']={'passed':True,**accounting_receipt()}
    source=ROOT/'arrays/trial_968_prefix_native.npz'
    receipt=ROOT/'receipts/tiny_wall_test.json'
    process=subprocess.run([sys.executable,str(UNIT/'native-qualification/native_curve_gate.py'),
                            str(source),'--output',str(receipt),'--identity-world','--wall-seconds','.01'],
                           capture_output=True,text=True,timeout=10)
    report=json.loads(receipt.read_text())
    assert process.returncode==2 and not report['qualified']
    assert report['status']=='unproven_cli_hard_wall_budget_exceeded'
    results['windows_cli_hard_wall_deadline']={'passed':True,'status':report['status']}
    report={'all_assertions_passed':True,'test_groups':len(results),'results':results,
            'scope':'Source and real-window structural tests only; no material promotion.'}
    (ROOT/'receipts/repair_tests.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':run()
