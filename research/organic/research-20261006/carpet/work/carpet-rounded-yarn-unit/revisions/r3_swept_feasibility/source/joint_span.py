"""Joint first two sections with true-wrapper witness constraints."""
import time,json,portable_resources
import numpy as np
from scipy.sparse import coo_matrix,block_diag,vstack
from scipy.optimize import least_squares
from swept_feasibility import ROOT,R2,SectionConstraints,segment_closest
from advance_probe import Field,UNIT

def wrappers(field,s,which):
    c=field.spline(s);v=field.spline(s,1);t=v/np.linalg.norm(v,axis=1)[:,None];inside=np.cross(field.normal,t);inside/=np.linalg.norm(inside,axis=1)[:,None];a=inside*np.cos(field.alpha)+field.normal*np.sin(field.alpha);b=np.cross(t,a)
    theta=field.phase0+field.rate*s+which*np.pi/3
    return c+.99*field.R*(np.cos(theta)[:,None]*a+np.sin(theta)[:,None]*b)

def wrapper_witnesses(field,a,b,end):
    owner=np.repeat(np.arange(187),24);which=np.tile(np.repeat(np.arange(6),4),187);part=np.tile(np.arange(4),187*6);lo=part*end/4;hi=(part+1)*end/4;u=b[owner]-a[owner];square=np.sum(u*u,axis=1)
    def evaluate(s):
        p=wrappers(field,s,which);parameter=np.clip(np.sum((p-a[owner])*u,axis=1)/np.maximum(square,1e-30),0,1);d=a[owner]+parameter[:,None]*u-p
        return np.sum(d*d,axis=1)
    original_lo=lo.copy();original_hi=hi.copy();ratio=(np.sqrt(5)-1)/2;x1=hi-ratio*(hi-lo);x2=lo+ratio*(hi-lo);f1=evaluate(x1);f2=evaluate(x2)
    for _ in range(42):
        left=f1<=f2;hi=np.where(left,x2,hi);lo=np.where(left,lo,x1);x1=hi-ratio*(hi-lo);x2=lo+ratio*(hi-lo);f1=evaluate(x1);f2=evaluate(x2)
    positions=np.stack((original_lo,original_hi,x1,x2));values=np.stack([evaluate(x) for x in positions]);best=np.argmin(values,axis=0);ids=np.arange(len(owner));s=positions[best,ids];distance=np.sqrt(values[best,ids]);choice=np.argmin(distance.reshape(187,6,4),axis=2).ravel();ids=np.arange(187*6)*4+choice
    return owner[ids],wrappers(field,s[ids],which[ids]),{'minimum_true_wrapper_margin_m':float(distance[ids].min()-field.r-16e-6),'closest_parameters_m':s[ids],'wrapper_indices':which[ids],'method':'Actual transported wrapper field evaluated directly; vectorized golden search on four parameter intervals, including every interval endpoint. This is a construction witness search, not the final native cubic certificate.'}

class JointSystem:
    def __init__(self,field,end,reference0,reference1,owner,witness):
        self.field=field;self.end=end;self.reference=np.concatenate((reference0,reference1))/field.R;self.margin=2e-6/field.R;self.movement=30e-6/field.R;self.parts=[];self.centers=[];self.bases=[]
        for s,reference in ((0,reference0),(end,reference1)):
            c,e1,e2,k=field.basis(s);self.centers.append(c/field.R);self.bases.append(np.stack((e1,e2),axis=1));self.parts.append(SectionConstraints(reference/field.R,k*field.R,field.rate*field.R,.045,16e-6/field.R,self.margin,movement=self.movement))
        self.i,self.j=np.triu_indices(187,1);self.owner=owner;self.witness=witness/field.R
    def world(self,x):
        x=x.reshape(2,187,2);return [self.centers[k]+x[k]@self.bases[k].T for k in range(2)]
    def residual_jacobian(self,flat,jacobian=True):
        x=flat.reshape(2,187,2);pieces=[self.parts[k].residual_jacobian(x[k].ravel(),jacobian) for k in range(2)]
        if jacobian:values=np.r_[pieces[0][0],pieces[1][0]];base=block_diag((pieces[0][1],pieces[1][1]),format='csr')
        else:values=np.r_[pieces[0],pieces[1]]
        a,b=self.world(flat);distance,s,t,n,error,_=segment_closest(a[self.i],b[self.i],a[self.j],b[self.j]);raw=.09+self.margin-(distance-error)
        u=b[self.owner]-a[self.owner];square=np.sum(u*u,axis=1);parameter=np.clip(np.sum((self.witness-a[self.owner])*u,axis=1)/np.maximum(square,1e-30),0,1);delta=a[self.owner]+parameter[:,None]*u-self.witness;length=np.linalg.norm(delta,axis=1);normal=delta/np.maximum(length[:,None],1e-30);raw_wrap=.045+16e-6/self.field.R+self.margin-length
        result=np.r_[values,np.maximum(raw,0),np.maximum(raw_wrap,0)]
        if not jacobian:return result
        rows=[];columns=[];entries=[]
        constraints=[(self.i,-(1-s)[:,None]*(n@self.bases[0]),0,raw>0),(187+self.i,-s[:,None]*(n@self.bases[1]),0,raw>0),(self.j,(1-t)[:,None]*(n@self.bases[0]),0,raw>0),(187+self.j,t[:,None]*(n@self.bases[1]),0,raw>0),(self.owner,-(1-parameter)[:,None]*(normal@self.bases[0]),len(raw),raw_wrap>0),(187+self.owner,-parameter[:,None]*(normal@self.bases[1]),len(raw),raw_wrap>0)]
        for ids,gradient,offset,active in constraints:
            for k in range(2):rows.append(np.arange(len(ids))+offset);columns.append(2*ids+k);entries.append(gradient[:,k]*active)
        extra=coo_matrix((np.concatenate(entries),(np.concatenate(rows),np.concatenate(columns))),shape=(len(raw)+len(raw_wrap),748)).tocsr();extra.eliminate_zeros();return result,vstack((base,extra),format='csr')

def run():
    start=time.monotonic();data=np.load(UNIT.parent/'carpet-r5-packed-core/receipts/baseline_loops.npz');field=Field(968,data);trace=np.load(UNIT/'arrays/constructor_trace_968.npz');qualified=np.load(R2/'arrays/section_0_candidate.npy');end=field.length/96
    reference0=trace['xy'][0];reference1=qualified+field.target(end)-field.target(0);current=np.concatenate((qualified,reference1))/field.R;refs=np.concatenate((reference0,reference1))/field.R
    initial_a=field.world(0,qualified);initial_b=field.world(end,reference1);owner,witness,initial=wrapper_witnesses(field,initial_a,initial_b,end);all_owner=owner.copy();all_witness=witness.copy();rounds=[];total_evaluations=0;status='unstarted'
    class Stop(Exception):pass
    for outer in range(4):
        if time.monotonic()-start>=30 or total_evaluations>=160:status='budget_exhausted';break
        system=JointSystem(field,end,reference0,reference1,all_owner,all_witness);best=current.ravel().copy();best_score=np.inf;cache=None;pair=None;round_evaluations=0
        def evaluate(z):
            nonlocal best,best_score,cache,pair,round_evaluations,total_evaluations
            if cache is not None and np.array_equal(cache,z):return pair
            if time.monotonic()-start>=30 or total_evaluations>=160:raise Stop
            pair=system.residual_jacobian(z);cache=z.copy();round_evaluations+=1;total_evaluations+=1;score=float(pair[0].max())
            if score<best_score:best_score=score;best=z.copy()
            if score<=.02e-6/field.R:raise Stop
            return pair
        try:
            least_squares(lambda z:evaluate(z)[0],current.ravel(),jac=lambda z:evaluate(z)[1],bounds=((refs-30e-6/field.R).ravel(),(refs+30e-6/field.R).ravel()),method='trf',tr_solver='lsmr',tr_options={'atol':1e-10,'btol':1e-10,'maxiter':300},x_scale='jac',max_nfev=min(60,160-total_evaluations),ftol=1e-12,xtol=1e-12,gtol=1e-12)
        except Stop:pass
        current=best.reshape(374,2);a,b=system.world(current);a=(a*field.R).astype('f4').astype('f8');b=(b*field.R).astype('f4').astype('f8');new_owner,new_witness,actual=wrapper_witnesses(field,a,b,end)
        recovered=np.concatenate(((a/field.R-system.centers[0])@system.bases[0],(b/field.R-system.centers[1])@system.bases[1]));final=system.residual_jacobian(recovered.ravel(),False);distance,_,_,_,error,_=segment_closest(a[system.i],b[system.i],a[system.j],b[system.j]);body_margin=float((distance-error-2*field.r).min());section=[system.parts[k].assess(recovered[k*187:(k+1)*187].ravel()) for k in range(2)]
        accepted=bool(final.max()*field.R<=.05e-6 and body_margin>=1.95e-6 and actual['minimum_true_wrapper_margin_m']>=1.95e-6)
        row={'round':outer+1,'evaluations':round_evaluations,'maximum_constrained_distance_residual_m':float(final.max()*field.R),'body_swept_margin_m':body_margin,'actual_curved_wrapper_margin_m':actual['minimum_true_wrapper_margin_m'],'section_reports_normalized':section,'accepted':accepted};rounds.append(row);print(json.dumps(row),flush=True)
        if accepted:status='qualified_straight_span';break
        all_owner=np.r_[all_owner,new_owner];all_witness=np.concatenate((all_witness,new_witness));status='witness_round_failed'
    np.savez_compressed(ROOT/'arrays/joint_span_968.npz',initial_xy=current[:187]*field.R,end_xy=current[187:]*field.R,start_world=a,end_world=b,reference0=reference0,reference1=reference1,end_arclength=end)
    result={'scope':'Joint first two section states, unchanged guide and morphology. Both remain constrained to their original 30 micrometre reference disks.','status':status,'span_arclength_m':end,'span_passed':status=='qualified_straight_span','initial_metric_only_proof_retained':True,'true_wrapper_witness_search':initial['method'],'initial_true_wrapper_gap_m':initial['minimum_true_wrapper_margin_m'],'rounds':rounds,'solves_per_accepted_span':len(rounds) if status=='qualified_straight_span' else None,'total_evaluations':total_evaluations,'seconds':time.monotonic()-start,'peak_rss_mib':portable_resources.getrusage(portable_resources.RUSAGE_SELF).ru_maxrss/1024,'limits':{'seconds':30,'total_evaluations':160,'witness_rounds':4,'movement_per_section_m':30e-6,'minimum_final_margin_m':1.95e-6},'native_curve_passed':False,'three_loop_feasibility_established':False}
    (ROOT/'receipts/joint_span_probe.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)

if __name__=='__main__':run()
