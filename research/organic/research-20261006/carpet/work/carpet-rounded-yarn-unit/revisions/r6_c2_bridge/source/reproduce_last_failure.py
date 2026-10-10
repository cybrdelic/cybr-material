"""One authorized exact failed-window reproduction; accepted prefix is not rerun."""
from pathlib import Path
import time,json,resource,copy
import numpy as np
from scipy.optimize import least_squares
from continue_bridge_loop import Window,witness_search,Field,UNIT
ROOT=Path(__file__).resolve().parents[1]
def run():
    begin=time.monotonic();deadline=begin+30;data=np.load(UNIT.parent/'carpet-r5-packed-core/receipts/baseline_loops.npz');field=Field(968,data);state=np.load(ROOT/'arrays/loop_968_partial.npz');original=json.loads((ROOT/'receipts/bounded_loop_968.json').read_text());expected=original['attempts'][-1]
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
        if accepted:break
        all_owner=np.r_[all_owner,new_owner];all_witness=np.concatenate((all_witness,new_witness))
    np.savez_compressed(ROOT/'arrays/last_rejected_window.npz',candidate_xy=current*field.R,recovered_xy=recovered*field.R,candidate_world=np.asarray(nodes),references=refs,previous_fixed=positions[-2],history_a=h_a,history_b=h_b,history_owner=h_owner,wrapper_owner=system.wo,wrapper_points=system.wp*field.R,stations=stations,start_parameter=start,end_parameter=end)
    full=copy.copy(system);full.i=system.all_i;full.j=system.all_j;values=full.evaluate(recovered.ravel(),False);groups=[];offset=0
    definitions=[]
    for node,part in enumerate(system.parts):
        for name,count in zip(('body_metric','wrap_metric','envelope','movement'),part.counts):definitions.append((f'section{node}_{name}',count,node,part,name))
    definitions.extend([('swept_variable_pairs',len(full.i),None,None,None),('fixed_history_pairs',len(full.hi),None,None,None),('curved_wrapper_witness_bank',len(full.wo),None,None,None)])
    for label,count,node,part,name in definitions:
        v=values[offset:offset+count];local=int(v.argmax()) if count else None;maximum=float(v.max(initial=0)*field.R);witness={'constraint':label,'maximum_distance_residual_m':maximum,'rows_above_0_05um':int((v*field.R>.05e-6).sum()),'local_row':local}
        if count and node is not None:
            x=recovered[node];witness['section_parameter_m']=(start,end)[node]
            if name=='body_metric':
                i=int(part.i[local]);j=int(part.j[local]);witness.update(fibres=[i,j],metric_surface_margin_m=2e-6-maximum,actual_same_section_euclidean_margin_m=float((np.linalg.norm(x[i]-x[j])-.09)*field.R))
            elif name=='wrap_metric':
                i=int(part.wi[local]);j=int(part.wj[local]);witness.update(fibre=i,wrap=j,metric_surface_margin_m=2e-6-maximum,actual_same_section_euclidean_margin_m=float((np.linalg.norm(x[i]-part.wraps[j])-.045-16e-6/field.R)*field.R))
            elif name=='envelope':witness.update(fibre=local,radial_center_distance_m=float(np.linalg.norm(x[local])*field.R),permitted_center_radius_m=.91*field.R,physical_surface_radius_m=float((np.linalg.norm(x[local])+.045)*field.R),original_bundle_radius_m=field.R)
            elif name=='movement':witness.update(fibre=local,actual_movement_m=float(np.linalg.norm(x[local]-part.reference[local])*field.R),movement_cap_m=30e-6)
        groups.append(witness);offset+=count
    matches=all(abs(a['max_distance_residual_m']-b['max_distance_residual_m'])<1e-15 and abs(a['body_swept_margin_m']-b['body_swept_margin_m'])<1e-15 and abs(a['curved_wrapper_witness_margin_m']-b['curved_wrapper_witness_margin_m'])<1e-15 for a,b in zip(rounds,expected['rounds'])) and len(rounds)==len(expected['rounds'])
    assert matches
    result={'scope':'Only the final0.681732micrometre rejected attempt reproduced with identical references and160total evaluations; accepted prefix not replayed.','scalar_receipt_reproduced_to_1e_15m':matches,'rounds':rounds,'constraint_groups':groups,'maximum_group':max(groups,key=lambda x:x['maximum_distance_residual_m']),'seconds':time.monotonic()-begin,'peak_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024}
    (ROOT/'receipts/rejected_window_constraints.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
if __name__=='__main__':run()
