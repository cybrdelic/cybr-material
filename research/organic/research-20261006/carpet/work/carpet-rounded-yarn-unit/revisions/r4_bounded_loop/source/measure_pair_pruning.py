"""Read-only equivalence and cost probe on identical saved window states."""
import copy,json,time,resource
import numpy as np
from certified_segments import ROOT,segment_closest
from advance_loop import Window,witness_search,Field,UNIT

def run():
    start=time.monotonic();data=np.load(UNIT.parent/'carpet-r5-packed-core/receipts/baseline_loops.npz');field=Field(968,data);saved=np.load(ROOT/'arrays/loop_968_partial.npz');stations=saved['stations'];xy=saved['xy'];world=saved['world'];references=saved['references'];rows=[]
    for end_index in (2,13,26):
        before=end_index-2;start_s=stations[end_index-1];end_s=stations[end_index];refs=references[end_index-1:end_index+1]
        history_a=world[:before].reshape(-1,3);history_b=world[1:before+1].reshape(-1,3);history_owner=np.tile(np.arange(187),before)
        a=np.concatenate((world[before],world[end_index-1]));b=np.concatenate((world[end_index-1],world[end_index]));low=max(0.,stations[before]-2*field.length/96);high=min(field.length,end_s+2*field.length/96)
        owner,witness,_=witness_search(field,a,b,low,high);full=Window(field,start_s,end_s,refs,world[before],history_a,history_b,history_owner,owner,witness)
        a_ref,b_ref=full.segments(full.refs.ravel());i=full.i;j=full.j;dist,_,_,_,error,_=segment_closest(a_ref[i],b_ref[i],a_ref[j],b_ref[j]);lower=dist-error
        # sqrt(2)*movement covers every point in the optimizer's coordinate box,
        # including infeasible iterates outside the circular movement constraint.
        displacement=np.sqrt(2)*full.move;threshold=.09+full.margin+2*displacement+1e-12;keep=lower<=threshold
        pruned=copy.copy(full);pruned.i=i[keep];pruned.j=j[keep]
        current=np.stack([(world[end_index-1+k]/field.R-full.centers[k])@full.bases[k] for k in range(2)]).ravel()
        base_rows=sum(sum(part.counts) for part in full.parts);full_value,full_jac=full.evaluate(current);selection=np.r_[np.arange(base_rows),base_rows+np.flatnonzero(keep),np.arange(base_rows+len(i),len(full_value))]
        compact_value,compact_jac=pruned.evaluate(current);difference=full_jac[selection]-compact_jac;max_jac=float(abs(difference.data).max()) if difference.nnz else 0
        assert np.array_equal(full_value[selection],compact_value) and max_jac==0 and np.all(full_value[base_rows+np.flatnonzero(~keep)]==0)
        # Check omitted rows on identical random optimizer-box states as well.
        random=np.random.default_rng(968+end_index);random_max=0.
        for _ in range(3):
            sample=full.refs.ravel()+random.uniform(-full.move,full.move,748);values=full.evaluate(sample,False);random_max=max(random_max,float(values[base_rows+np.flatnonzero(~keep)].max(initial=0)))
        assert random_max==0
        elapsed={}
        for label,system in (('unculled',full),('culled',pruned)):
            times=[]
            for _ in range(4):
                tick=time.perf_counter();system.evaluate(current);times.append(time.perf_counter()-tick)
            elapsed[label]=float(np.median(times))
        row={'end_section':end_index,'saved_arclength_m':float(end_s),'unculled_variable_segment_pairs':len(i),'retained_variable_segment_pairs':int(keep.sum()),'removed_fraction':float((~keep).mean()),'section_rows_not_pruned':base_rows,'fixed_history_candidates_unchanged':len(full.hi),'all_retained_residuals_bitwise_equal':True,'retained_jacobian_max_difference':max_jac,'maximum_omitted_residual_saved_or_random_box_states':random_max,'timing_median_seconds':elapsed,'residual_jacobian_speedup':elapsed['unculled']/elapsed['culled'],'movement_bound_used_normalized':displacement,'movement_bound_used_m':displacement*field.R,'proof':'For each segment, every point moves at most the largest endpoint displacement. Segment-set distance changes by at most the sum of both segment bounds. The conservative reference lower bound exceeds target plus both bounds for every omitted pair. sqrt(2) times the original30micrometre cap covers the full optimizer box; physical movement acceptance remains30micrometres.'}
        rows.append(row);print(json.dumps(row),flush=True)
    report={'scope':'No construction rerun. Conservative candidate reduction on three identical saved window states; morphology and solver budgets unchanged.','rows':rows,'full_final_check_required':True,'loop_or_crop_runtime_not_measured':True,'seconds':time.monotonic()-start,'peak_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024}
    (ROOT/'receipts/pair_pruning_probe.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:v for k,v in report.items() if k!='rows'}))

if __name__=='__main__':run()
