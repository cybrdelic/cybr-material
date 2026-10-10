"""Bounded moving two-section construction with complete body-history checks."""
import time,json,resource
import numpy as np
from scipy.sparse import coo_matrix,block_diag,vstack
from scipy.spatial import cKDTree
from scipy.optimize import least_squares
from certified_segments import ROOT,R3,segment_closest
from joint_span import wrappers
from swept_feasibility import SectionConstraints
from advance_probe import Field,UNIT

def witness_search(field,a,b,lo_s,hi_s):
    count=len(a);owner=np.repeat(np.arange(count),24);which=np.tile(np.repeat(np.arange(6),4),count);part=np.tile(np.arange(4),count*6)
    lo=lo_s+part*(hi_s-lo_s)/4;hi=lo_s+(part+1)*(hi_s-lo_s)/4;original_lo=lo.copy();original_hi=hi.copy();u=b[owner]-a[owner];length=np.sum(u*u,axis=1)
    def evaluate(s):
        p=wrappers(field,s,which);v=np.clip(np.sum((p-a[owner])*u,axis=1)/np.maximum(length,1e-30),0,1);delta=a[owner]+v[:,None]*u-p;return np.sum(delta*delta,axis=1)
    ratio=(np.sqrt(5)-1)/2
    for _ in range(36):
        x1=hi-ratio*(hi-lo);x2=lo+ratio*(hi-lo);left=evaluate(x1)<=evaluate(x2);hi=np.where(left,x2,hi);lo=np.where(left,lo,x1)
    samples=np.stack((original_lo,original_hi,lo,hi));values=np.stack([evaluate(s) for s in samples]);choice=np.argmin(values,axis=0);ids=np.arange(len(owner));parameters=samples[choice,ids];distance=np.sqrt(values[choice,ids]);choose=np.argmin(distance.reshape(count,6,4),axis=2).ravel();ids=np.arange(count*6)*4+choose
    return owner[ids],wrappers(field,parameters[ids],which[ids]),float(distance[ids].min()-field.r-16e-6)

class Window:
    def __init__(self,field,start,end,refs,fixed_start,history_a,history_b,history_owner,wo,wp):
        self.field=field;self.refs=refs/field.R;self.margin=2e-6/field.R;self.move=30e-6/field.R;self.parts=[];self.centers=[];self.bases=[];self.fixed=fixed_start/field.R
        for s,ref in zip((start,end),refs):
            c,e1,e2,k=field.basis(s);self.centers.append(c/field.R);self.bases.append(np.stack((e1,e2),axis=1));self.parts.append(SectionConstraints(ref/field.R,k*field.R,field.rate*field.R,.045,16e-6/field.R,self.margin,movement=self.move))
        self.left_ids=np.r_[np.full(187,-1),np.arange(187)];self.right_ids=np.r_[np.arange(187),187+np.arange(187)];self.owners=np.tile(np.arange(187),2)
        i,j=np.triu_indices(374,1);valid=self.owners[i]!=self.owners[j];self.i=i[valid];self.j=j[valid]
        self.ha=history_a/field.R;self.hb=history_b/field.R;self.howner=history_owner;self.wo=wo;self.wp=wp/field.R
        a,b=self.segments(self.refs.ravel());self.hi=np.empty(0,dtype=int);self.hj=np.empty(0,dtype=int)
        if len(self.ha):
            mid=(a+b)*.5;hm=(self.ha+self.hb)*.5;half=np.linalg.norm(b-a,axis=1)*.5;hh=np.linalg.norm(self.hb-self.ha,axis=1)*.5
            # Both moving endpoints can shift by the original movement cap.
            radii=.09+self.margin+half+hh.max()+2*self.move+2*.05e-6/field.R
            lists=cKDTree(hm).query_ball_point(mid,radii,workers=1);i=np.repeat(np.arange(374),[len(x) for x in lists]);j=np.concatenate([np.asarray(x,dtype=int) for x in lists]) if len(i) else np.empty(0,dtype=int);keep=self.owners[i]!=history_owner[j];self.hi=i[keep];self.hj=j[keep]
    def nodes(self,flat):
        x=flat.reshape(2,187,2);return [self.centers[k]+x[k]@self.bases[k].T for k in range(2)]
    def segments(self,flat):
        a,b=self.nodes(flat);return np.concatenate((self.fixed,a)),np.concatenate((a,b))
    def evaluate(self,flat,jacobian=True):
        x=flat.reshape(2,187,2);base=[self.parts[k].residual_jacobian(x[k].ravel(),jacobian) for k in range(2)]
        if jacobian:residual=np.r_[base[0][0],base[1][0]];matrix=block_diag((base[0][1],base[1][1]),format='csr')
        else:residual=np.r_[base[0],base[1]]
        a,b=self.segments(flat);blocks=[];jacblocks=[]
        def add_block(values,terms):
            nonlocal blocks,jacblocks
            blocks.append(np.maximum(values,0))
            if not jacobian:return
            rows=[];columns=[];entries=[];active=values>0
            for ids,world_gradient in terms:
                live=ids>=0
                if not live.any():continue
                ids0=ids[live];basis=np.where((ids0<187)[:,None,None],self.bases[0],self.bases[1]);gradient=np.einsum('ij,ijk->ik',world_gradient[live],basis)
                for k in range(2):rows.append(np.flatnonzero(live));columns.append(2*ids0+k);entries.append(gradient[:,k]*active[live])
            extra=coo_matrix((np.concatenate(entries),(np.concatenate(rows),np.concatenate(columns))),shape=(len(values),748)).tocsr();extra.eliminate_zeros();jacblocks.append(extra)
        d,s,t,n,error,_=segment_closest(a[self.i],b[self.i],a[self.j],b[self.j]);raw=.09+self.margin-(d-error)
        add_block(raw,[(self.left_ids[self.i],-(1-s)[:,None]*n),(self.right_ids[self.i],-s[:,None]*n),(self.left_ids[self.j],(1-t)[:,None]*n),(self.right_ids[self.j],t[:,None]*n)])
        if len(self.hi):
            d,s,t,n,error,_=segment_closest(a[self.hi],b[self.hi],self.ha[self.hj],self.hb[self.hj]);raw=.09+self.margin-(d-error);add_block(raw,[(self.left_ids[self.hi],-(1-s)[:,None]*n),(self.right_ids[self.hi],-s[:,None]*n)])
        u=b[self.wo]-a[self.wo];square=np.sum(u*u,axis=1);parameter=np.clip(np.sum((self.wp-a[self.wo])*u,axis=1)/np.maximum(square,1e-30),0,1);delta=a[self.wo]+parameter[:,None]*u-self.wp;length=np.linalg.norm(delta,axis=1);n=delta/np.maximum(length[:,None],1e-30);raw=.045+16e-6/self.field.R+self.margin-length
        add_block(raw,[(self.left_ids[self.wo],-(1-parameter)[:,None]*n),(self.right_ids[self.wo],-parameter[:,None]*n)])
        result=np.concatenate([residual]+blocks)
        return (result,vstack([matrix]+jacblocks,format='csr')) if jacobian else result
    def body_margin(self,flat):
        a,b=self.segments(flat);d,s,t,n,e,_=segment_closest(a[self.i],b[self.i],a[self.j],b[self.j]);minimum=float((d-e-.09).min())
        if len(self.hi):
            d,s,t,n,e,_=segment_closest(a[self.hi],b[self.hi],self.ha[self.hj],self.hb[self.hj]);minimum=min(minimum,float((d-e-.09).min()))
        return minimum*self.field.R

def run():
    begin=time.monotonic();deadline=begin+120;data=np.load(UNIT.parent/'carpet-r5-packed-core/receipts/baseline_loops.npz');field=Field(968,data);proof=np.load(R3/'arrays/joint_span_968.npz');base_step=field.length/96
    stations=[0.,base_step];xy=[proof['initial_xy'].copy(),proof['end_xy'].copy()];positions=[proof['start_world'].copy(),proof['end_world'].copy()];references=[proof['reference0'].copy(),proof['reference1'].copy()];attempts=[];step=base_step;status='running';total_evaluations=0
    class Stop(Exception):pass
    while stations[-1]<field.length-1e-12 and time.monotonic()<deadline:
        if len(stations)>=385:status='384_span_limit';break
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
            values=system.evaluate(recovered.ravel(),False);body=system.body_margin(recovered.ravel());a,b=system.segments(recovered.ravel());new_owner,new_witness,wrap=witness_search(field,a*field.R,b*field.R,low,high)
            accepted=bool(values.max()*field.R<=.05e-6 and body>=1.95e-6 and wrap>=1.95e-6);rounds.append({'round':outer+1,'evaluations':round_evaluations,'max_distance_residual_m':float(values.max()*field.R),'body_swept_margin_m':body,'curved_wrapper_witness_margin_m':wrap,'body_variable_pairs':len(system.i),'body_history_candidates':len(system.hi),'complete_body_broad_phase':True,'accepted':accepted})
            if accepted:break
            all_owner=np.r_[all_owner,new_owner];all_witness=np.concatenate((all_witness,new_witness))
        row={'start_arclength_m':start,'end_arclength_m':end,'span_um':(end-start)*1e6,'accepted':accepted,'rounds':rounds,'evaluations':attempt_evaluations,'seconds':time.monotonic()-attempt_start};attempts.append(row)
        if accepted:
            xy[-1]=current[0]*field.R;xy.append(current[1]*field.R);positions[-1]=nodes[0];positions.append(nodes[1]);stations.append(end);references.append(refs[1]);step=min(base_step,step*2)
        else:
            step*=.5
            if step<base_step/64-1e-14:status='minimum_span_failed';break
        print(json.dumps({'accepted_spans':len(stations)-1,'length_fraction':stations[-1]/field.length,'attempt_seconds':row['seconds'],'evaluations':attempt_evaluations,'accepted':accepted,'span_um':row['span_um']}),flush=True)
        np.savez_compressed(ROOT/'arrays/loop_968_partial.npz',stations=stations,xy=xy,world=positions,references=references)
        partial={'status':'running','attempts':attempts,'accepted_arclength_m':stations[-1],'accepted_length_fraction':stations[-1]/field.length,'seconds':time.monotonic()-begin};(ROOT/'receipts/loop_968_progress.json').write_text(json.dumps(partial,indent=2))
    if stations[-1]>=field.length-1e-12:status='piecewise_linear_body_passed_native_gate_pending'
    elif status=='running':status='120_second_bound'
    result={'scope':'Loop 968 bounded overlapping-window advance; guide/morphology unchanged; no scene.','status':status,'accepted_spans':len(stations)-1,'accepted_arclength_m':stations[-1],'guide_length_m':field.length,'accepted_length_fraction':stations[-1]/field.length,'evaluations':total_evaluations,'attempts':attempts,'seconds':time.monotonic()-begin,'peak_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'continuous_wrapper_certificate':False,'same_fibre_self_clearance_qualified':False,'native_curve_qualified':False,'limits':{'total_seconds':120,'per_attempt_seconds':30,'per_attempt_evaluations':160,'minimum_accepted_margin_m':1.95e-6,'movement_reference_disks_m':30e-6}}
    (ROOT/'receipts/bounded_loop_968.json').write_text(json.dumps(result,indent=2));print(json.dumps({k:v for k,v in result.items() if k!='attempts'}))

if __name__=='__main__':run()
