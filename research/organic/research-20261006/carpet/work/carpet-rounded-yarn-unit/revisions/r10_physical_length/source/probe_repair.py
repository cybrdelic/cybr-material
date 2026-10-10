"""One bounded trial from the last frozen48-span state, never span39."""
import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
from scipy.optimize import least_squares
from physical_length_window import (ROOT,UNIT,WORK,BridgeField,PhysicalLengthWindow,
                                   DISK_M,ROUNDING_INSET_M)
sys.path.insert(0,str(UNIT/'native-qualification'))
import native_curve_gate as native
from portable_resources import peak_rss_mib,accounting_receipt

def digest(path):return hashlib.file_digest(Path(path).open('rb'),'sha256').hexdigest()

def setup():
    provenance=json.loads((WORK.parent/'receipts/guide_regeneration.json').read_text())
    assert provenance['status']=='deterministic_guide_regenerated_and_restart_checks_passed'
    assert provenance['restart_checks']['accepted_float32_world_max_error_m']==0
    data=np.load(WORK/'carpet-r5-packed-core/receipts/baseline_loops.npz')
    field=BridgeField(968,data)
    frozen=UNIT/'revisions/r9_physical_resume/arrays/loop_968_partial.npz'
    rejected=UNIT/'revisions/r9_physical_resume/arrays/rejected_attempt_0025.npz'
    assert digest(frozen)=='ab087996ecf9c0bc26cff51520e4ab3ca26b40ff7cc19bae228d3de163839cee'
    assert digest(rejected)=='ecba003421ca383baffe9d447130c6a0de294f74e126c4e905fb964925eb51ba'
    state=np.load(frozen);last=np.load(rejected)
    assert len(state['stations'])==49
    assert np.array_equal(state['references'][-1],last['references'][0])
    assert np.array_equal(state['world'][-2],last['previous_fixed'])
    stations=np.r_[state['stations'],last['end_parameter']]
    xy=np.concatenate((state['xy'],last['candidate_xy'][1:]))
    xy[-2:]=last['candidate_xy']
    world=np.concatenate((state['world'],last['candidate_world'][1:]))
    world[-2:]=last['candidate_world']
    refs=np.concatenate((state['references'],last['references'][1:]))
    expected=state['xy'][-1]+field.target(stations[-1])-field.target(stations[-2])
    assert np.array_equal(expected,refs[-1]) or np.max(abs(expected-refs[-1]))<1e-17
    plan=json.loads((UNIT/'revisions/r9_physical_resume/receipts/physical_window_repair_plan.json').read_text())
    # Native fixed wrapper keys come from the qualified deterministic guide.
    # They are an explicit new representation for this trial, not recovery of
    # missing historical wrapper arrays. Full guide includes distant history.
    parameters=np.linspace(0,field.length,1025)
    wraps=field.wrapper_array(parameters).astype('f4')
    config=native.GateConfig(wall_seconds=60,max_leaves=200000)
    leaves=[];budget=native.Budget(config)
    for owner,keys in enumerate(wraps):
        part,_=native.owner_leaves(native.native_beziers(keys),owner,np.inf,config,budget,1e-12,200000)
        leaves.append(part)
    wrapper_leaves={key:np.concatenate([v[key] for v in leaves]) for key in ('a','b','error')}
    system=PhysicalLengthWindow(field,stations,xy,world,refs,plan['contact_connected_fibres'],
                               plan['physical_window_target_m'],wrapper_leaves)
    assert system.anchor==plan['fixed_anchor_index']
    assert system.n==plan['movable_point_count']
    assert system.affected.sum()==plan['affected_body_segments']
    assert abs(stations[-1]-stations[system.anchor]-plan['actual_window_width_m'])<1e-14
    return system,wraps,parameters,{'frozen_state_sha256':digest(frozen),'rejected_state_sha256':digest(rejected),
                                  'guide_regeneration_receipt_sha256':digest(WORK.parent/'receipts/guide_regeneration.json'),
                                  'wrapper_leaf_count':len(wrapper_leaves['a']),
                                  'wrapper_representation':'new deterministic1025-key float32 full-guide wrapper representation; not exact historical recovery'}

def run():
    parser=argparse.ArgumentParser()
    parser.add_argument('--seconds',type=float,default=60)
    parser.add_argument('--evaluations',type=int,default=80)
    parser.add_argument('--setup-only',action='store_true')
    args=parser.parse_args()
    start=time.monotonic();system,wraps,parameters,inputs=setup()
    (ROOT/'arrays').mkdir(exist_ok=True);(ROOT/'receipts').mkdir(exist_ok=True)
    setup_seconds=time.monotonic()-start
    x=system.initial.copy()
    delta=x-system.refs;norm=np.linalg.norm(delta,axis=1)
    x=system.refs+delta*np.minimum(1,(DISK_M-ROUNDING_INSET_M)/system.field.R/np.maximum(norm,1e-30))[:,None]
    best=x.ravel().copy();best_score=np.inf;count=0;cache=None;pair=None;termination='evaluation_limit'
    class Stop(Exception):pass
    def evaluate(z):
        nonlocal best,best_score,count,cache,pair,termination
        if cache is not None and np.array_equal(cache,z):return pair
        if time.monotonic()-start>args.seconds:termination='wall_budget';raise Stop
        if peak_rss_mib()>900:termination='memory_budget';raise Stop
        if count>=args.evaluations:raise Stop
        pair=system.evaluate(z);cache=z.copy();count+=1;score=float(pair[0].max())
        if score<best_score:best=z.copy();best_score=score
        if count%10==0:print(json.dumps({'evaluations':count,'best_max_residual_nm':best_score*system.field.R*1e9,'elapsed_seconds':time.monotonic()-start}),flush=True)
        if score*system.field.R<=.5e-9:termination='residual_target';raise Stop
        return pair
    if args.setup_only:termination='setup_only';system.evaluate(best,False)
    else:
        try:
            bound=DISK_M/system.field.R
            result=least_squares(lambda z:evaluate(z)[0],best,jac=lambda z:evaluate(z)[1],
                bounds=((system.refs-bound).ravel(),(system.refs+bound).ravel()),
                method='trf',tr_solver='lsmr',tr_options={'atol':1e-9,'btol':1e-9,'maxiter':150},
                x_scale='jac',max_nfev=args.evaluations,ftol=1e-12,xtol=1e-12,gtol=1e-12)
            termination=f'scipy_status_{result.status}'
        except Stop:pass
    actual,checks=system.assess_actual(best)
    accepted=all(checks[name] for name in ('positive_physical_clearance_pass','construction_reserve_pass',
                        'continuous_straight_containment_pass','original30um_disks_pass','fixed_world_bytes_preserved'))
    candidate=ROOT/'arrays/trial_968_prefix_native.npz'
    np.savez_compressed(candidate,body=actual.transpose(1,0,2),wraps=wraps,
                        body_radius=np.float32(system.field.r),wrap_radius=np.float32(16e-6),
                        stations=system.stations,wrapper_parameters=parameters,
                        references=system.references,movable_mask=system.mask,
                        original_partial_world=system.original_world.astype('f4'))
    np.savez_compressed(ROOT/'arrays/optimizer_state.npz',xy=best.reshape(system.n,2)*system.field.R,
                        original_references=system.references,mask=system.mask,stations=system.stations)
    report={'scope':'One selective physical-length-window trial from final48-span state, loop968 only.',
            'status':'straight_trial_passed_native_prefix_gate_pending' if accepted else 'straight_trial_unproven',
            'production_loop_complete':False,'initial_spans':48,'trial_spans':len(system.stations)-1,
            'progress_if_accepted_fraction':float(system.stations[-1]/system.field.length),
            'guide_length_m':float(system.field.length),'anchor_index':system.anchor,
            'window_width_m':float(system.stations[-1]-system.stations[system.anchor]),
            'window_parameter_width_m':float(system.stations[-1]-system.stations[system.anchor]),
            'window_integrated_physical_width_m':system.physical_window_width_m,
            'physical_arc_quadrature_error_estimate_m':system.arc_integration_error_m,
            'movable_points':system.n,'scalar_unknowns':2*system.n,
            'affected_body_segments':int(system.affected.sum()),'fixed_body_segments':int((~system.affected).sum()),
            'complete_body_pairs':len(system.i),'complete_wrapper_pairs':len(system.wi),
            'complete_fixed_history_included':True,'original_references_preserved':True,
            'roots_and_apex_unmodified':True,'same_fibre_self_clearance_qualified':False,
            'native_curve_qualified':False,'straight_trial_accepted':accepted,'actual_float32_checks':checks,
            'evaluations':count,'termination':termination,'setup_seconds':setup_seconds,
            'source_sha256':{name:digest(Path(__file__).with_name(name)) for name in
                             ('probe_repair.py','physical_length_window.py')},
            'elapsed_seconds':time.monotonic()-start,'rss_accounting':accounting_receipt(),
            'limits':{'solver_seconds':args.seconds,'solver_evaluations':args.evaluations,'memory_mib':900,
                      'original_movement_disks_m':DISK_M,'construction_reserve_m':1.95e-6,
                      'solver_rounding_inset_m':ROUNDING_INSET_M},'inputs':inputs,
            'candidate_sha256':digest(candidate),
            'limitations':['No complete body loop or three-loop assembly exists.',
                           'Straight-body sufficient conditions cannot qualify native Catmull-Rom fibres.',
                           'Fixed native wrapper keys are explicitly regenerated source-derived trial geometry.']}
    (ROOT/'receipts/repair_result.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':run()
