"""Exact saved failed span; only physical contact/containment authority changes."""
from pathlib import Path
import time,json,resource,hashlib
import numpy as np
from scipy.optimize import least_squares
from physical_window import ROOT,PhysicalWindow as Window,Field,UNIT,witness_search,span_envelope_upper
def run():
    begin=time.monotonic();deadline=begin+30;data=np.load(UNIT.parent/'carpet-r5-packed-core/receipts/baseline_loops.npz');field=Field(968,data);R6=ROOT.parent/'r6_c2_bridge';state=np.load(R6/'arrays/loop_968_partial.npz');original=json.loads((R6/'receipts/bounded_loop_968.json').read_text());expected=original['attempts'][-1]
    stations=list(state['stations']);xy=[v.copy() for v in state['xy']];positions=[v.copy() for v in state['world']];references=[v.copy() for v in state['references']];base_step=field.length/96;step=expected['end_arclength_m']-expected['start_arclength_m'];total_evaluations=0
    class Stop(Exception):pass
    start=stations[-1];end=min(start+step,field.length);refs=np.stack((references[-1],xy[-1]+field.target(end)-field.target(start)));current=np.stack((xy[-1],refs[1]))/field.R
    h_a=np.concatenate(positions[:-2]) if len(positions)>2 else np.empty((0,3));h_b=np.concatenate(positions[1:-1]) if len(positions)>2 else np.empty((0,3));h_owner=np.tile(np.arange(187),max(0,len(positions)-2))
    candidate_a=np.concatenate((positions[-2],field.world(start,xy[-1])));candidate_b=np.concatenate((field.world(start,xy[-1]),field.world(end,refs[1])))
    low=max(0.,stations[-2]-2*base_step);high=min(field.length,end+2*base_step);owner,witness,initial_wrap=witness_search(field,candidate_a,candidate_b,low,high);all_owner=owner.copy();all_witness=witness.copy();rounds=[];attempt_start=time.monotonic();attempt_evaluations=0;accepted=False
    for outer in range(3):
        if time.monotonic()>=deadline or attempt_evaluations>=160:break
        system=Window(field,start,end,refs,positions[-2],h_a,h_b,h_owner,all_owner,all_witness);best=current.ravel().copy();best_score=np.inf;cache=None;pair=None;round_evaluations=0
        def evaluate(z):
            nonlocal best,best_score,cache,pair,round_evaluations,total_evaluations,attempt_evaluations
            if cache is not None and np.array_equal(cache,z):return pair
            if time.monotonic()>=deadline or time.monotonic()-attempt_start>=30 or attempt_evaluations>=160:raise Stop
            pair=system.evaluate(z);cache=z.copy();round_evaluations+=1;total_evaluations+=1;attempt_evaluations+=1;score=float(pair[0].max())
            if score<best_score:best_score=score;best=z.copy()
            if score<=.02e-6/field.R:raise Stop
            return pair
        try:
            least_squares(lambda z:evaluate(z)[0],current.ravel(),jac=lambda z:evaluate(z)[1],bounds=((refs/field.R-30e-6/field.R).ravel(),(refs/field.R+30e-6/field.R).ravel()),method='trf',tr_solver='lsmr',tr_options={'atol':1e-10,'btol':1e-10,'maxiter':300},x_scale='jac',max_nfev=min(60,160-attempt_evaluations),ftol=1e-12,xtol=1e-12,gtol=1e-12)
        except Stop:pass
        current=best.reshape(2,187,2);nodes=[(p*field.R).astype('f4').astype('f8') for p in system.nodes(best)];recovered=np.stack([(nodes[k]/field.R-system.centers[k])@system.bases[k] for k in range(2)])
        values=system.full_evaluate(recovered.ravel(),False);body=system.body_margin(recovered.ravel());a,b=system.segments(recovered.ravel());new_owner,new_witness,wrap=witness_search(field,a*field.R,b*field.R,low,high)
        accepted=bool(values.max()*field.R<=.05e-6 and body>=1.95e-6 and wrap>=1.95e-6);rounds.append({'round':outer+1,'evaluations':round_evaluations,'max_distance_residual_m':float(values.max()*field.R),'body_swept_margin_m':body,'curved_wrapper_witness_margin_m':wrap,'body_variable_pairs_retained':len(system.i),'body_variable_pairs_full_check':len(system.all_i),'body_history_candidates':len(system.hi),'complete_body_broad_phase':True,'accepted':accepted})
        np.savez_compressed(ROOT/f'arrays/round_{outer+1}_optimizer.npz',candidate_xy=current*field.R,recovered_xy=recovered*field.R,candidate_world=np.asarray(nodes),references=refs,previous_fixed=positions[-2],history_a=h_a,history_b=h_b,history_owner=h_owner,wrapper_owner=system.wo,wrapper_points=system.wp*field.R,stations=stations,start_parameter=start,end_parameter=end)
        if accepted:break
        all_owner=np.r_[all_owner,new_owner];all_witness=np.concatenate((all_witness,new_witness))

    envelope=[span_envelope_upper(field,stations[-2],start,positions[-2],nodes[0]),span_envelope_upper(field,start,end,nodes[0],nodes[1])]
    diagnostics=[part.assess(recovered[k].ravel()) for k,part in enumerate(system.parts)]
    prior_metric=[part.assess(recovered[k].ravel()) for k,part in enumerate(system.metric_diagnostics)]
    rejected=np.load(R6/'arrays/last_rejected_window.npz');same_refs=np.array_equal(refs,rejected['references']);same_history=np.array_equal(h_a,rejected['history_a']) and np.array_equal(h_b,rejected['history_b']);assert same_refs and same_history
    result={'scope':'Only the exact final failed span. Correctly differentiated infinite-tangent section metric remains diagnostic, not hard contact authority. Actual swept body/history and wrapper witnesses remain authoritative.','accepted':accepted and all(x['physical_containment_pass'] for x in envelope),'references_bitwise_unchanged':same_refs,'history_bitwise_unchanged':same_history,'guide_phase_seed_targets_unchanged':True,'span_parameters_m':[start,end],'swept_target_m':2e-6,'minimum_accepted_swept_margin_m':1.95e-6,'movement_disk_m':30e-6,'physical_center_containment_limit_m':field.R-field.r,'former_internal_center_packing_limit_m':.91*field.R,'rounds':rounds,'physical_section_checks':diagnostics,'former_metric_diagnostics_on_same_output':prior_metric,'continuous_straight_span_envelope':envelope,'native_curve_qualified':False,'seconds':time.monotonic()-begin,'peak_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'source').glob('*.py')}}
    np.savez_compressed(ROOT/'arrays/physical_span_candidate.npz',xy=current*field.R,recovered_xy=recovered*field.R,world=np.asarray(nodes),references=refs,start_parameter=start,end_parameter=end)
    (ROOT/'receipts/physical_span_probe.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
if __name__=='__main__':run()
