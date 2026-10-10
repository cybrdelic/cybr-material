"""Fail-closed section advancement with complete swept-segment searches."""
from pathlib import Path
import sys,json,time,portable_resources as resource
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.spatial import cKDTree
from section_feasibility import SectionConstraints,solve_section
ROOT=Path(__file__).resolve().parents[1];UNIT=ROOT.parents[1]
sys.path.insert(0,str(UNIT/'source'))
import construct_unit as old

MINIMUM=1.95e-6
RESIDUAL=0.05e-6

def exact_segments(a,b,c,d):
    """Exact finite-segment minimum: interior stationary pair or four edges."""
    u=b-a;v=d-c;w=a-c;aa=np.sum(u*u,axis=-1);bb=np.sum(u*v,axis=-1);cc=np.sum(v*v,axis=-1);dd=np.sum(u*w,axis=-1);ee=np.sum(v*w,axis=-1);det=aa*cc-bb*bb
    candidates=[]
    for point in (a,b):
        t=np.clip(np.sum((point-c)*v,axis=-1)/np.maximum(cc,1e-30),0,1);candidates.append(np.sum((point-c-t[...,None]*v)**2,axis=-1))
    for point in (c,d):
        s=np.clip(np.sum((point-a)*u,axis=-1)/np.maximum(aa,1e-30),0,1);candidates.append(np.sum((point-a-s[...,None]*u)**2,axis=-1))
    s=np.divide(bb*ee-cc*dd,det,out=np.zeros_like(det),where=det>1e-30);t=np.divide(aa*ee-bb*dd,det,out=np.zeros_like(det),where=det>1e-30)
    interior=(det>1e-30)&(s>=0)&(s<=1)&(t>=0)&(t<=1);distance=np.sum((w+s[...,None]*u-t[...,None]*v)**2,axis=-1);candidates.append(np.where(interior,distance,np.inf))
    return np.sqrt(np.minimum.reduce(candidates))

def complete_cross(a,b,owner,c,d,other_owner,threshold,exclude_same):
    if not len(c):return np.inf,0
    middle=(a+b)*.5;other=(c+d)*.5;half=np.linalg.norm(b-a,axis=1)*.5;other_half=np.linalg.norm(d-c,axis=1)*.5
    tree=cKDTree(other);lists=tree.query_ball_point(middle,threshold+half+other_half.max()+1e-12,workers=1)
    first=np.repeat(np.arange(len(a)),[len(x) for x in lists]);second=np.concatenate([np.asarray(x,dtype='i8') for x in lists]) if len(first) else np.empty(0,dtype='i8')
    if exclude_same:
        valid=owner[first]!=other_owner[second];first=first[valid];second=second[valid]
    if not len(first):return np.inf,0
    return float(exact_segments(a[first],b[first],c[second],d[second]).min()),len(first)

class Field:
    def __init__(self,index,data):
        self.index=index;self.data=data;self.R=float(data['radius'][index].mean());self.r=.045*self.R
        self.guide=np.load(UNIT/f'arrays/rounded_unit_{index}.npz')['guide'];arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(self.guide,axis=0),axis=1))];self.length=arc[-1];self.spline=CubicSpline(arc,self.guide,axis=0)
        c=data['center'][index];axis=old.unit(c[-1]-c[0]);self.normal=old.unit(np.cross(axis,c[24]-(c[0]+c[-1])*.5));t=old.unit(self.spline(0,1));inside=old.unit(np.cross(self.normal,t));a0=data['frame'][index,0,0];self.alpha=np.arctan2(a0@self.normal,a0@inside)
        self.phase0=float(data['phase'][index,0]);self.rate=float((data['phase'][index,-1]-self.phase0)/self.length)
        rng=np.random.default_rng(241016+index*104729);self.raw_seed=old.migration.old.PACK*.13*self.R+rng.normal(0,19e-6,(187,2));self.amplitude=rng.uniform(30e-6,60e-6,(187,2));self.wave=rng.uniform(.0006,.0012,(187,2));self.offset=rng.uniform(0,2*np.pi,(187,2))
    def basis(self,s):
        c=self.spline(s);velocity=self.spline(s,1);t=old.unit(velocity);a=old.unit(np.cross(self.normal,t))*np.cos(self.alpha)+self.normal*np.sin(self.alpha);b=old.unit(np.cross(t,a));phase=self.phase0+self.rate*s;e1=np.cos(phase)*a+np.sin(phase)*b;e2=-np.sin(phase)*a+np.cos(phase)*b
        accel=self.spline(s,2);bend=(accel-t*(accel@t))/(velocity@velocity);k=np.array([bend@e1,bend@e2]);return c,e1,e2,k
    def target(self,s):return np.sin(np.pi*s/self.length)**2*self.amplitude*(np.sin(2*np.pi*s/self.wave+self.offset)-np.sin(self.offset))*.5
    def world(self,s,xy):
        c,e1,e2,k=self.basis(s);return c+xy[:,0,None]*e1+xy[:,1,None]*e2
    def wraps(self,s):
        angle=np.arange(6)*np.pi/3;return self.world(s,.99*self.R*np.stack((np.cos(angle),np.sin(angle)),axis=1))

def check_section(field,s,xy,reference):
    c,e1,e2,k=field.basis(s);world=field.world(s,xy).astype('f4').astype('f8');recovered=np.stack(((world-c)@e1,(world-c)@e2),axis=1)/field.R
    system=SectionConstraints(reference/field.R,k*field.R,field.rate*field.R,.045,16e-6/field.R,2e-6/field.R,movement=30e-6/field.R)
    result=system.assess(recovered.ravel());violations=result['constraint_violations_normalized'];accepted=all(v*field.R<=RESIDUAL for v in violations.values())
    return accepted,{'accepted':accepted,'distance_residuals_normalized':violations,'distance_residuals_m':{key:value*field.R for key,value in violations.items()},'minimum_body_distance_margin_m':result['minimum_body_metric_gap_normalized']*field.R,'minimum_wrap_distance_margin_m':result['minimum_wrap_metric_gap_normalized']*field.R,'float32_world_positions_checked':True,'residual_kind':'ordinary distance, not squared distance'},world

def swept(field,s0,s1,a,b,history_a,history_b,history_owner,wrap_a,wrap_b):
    i,j=np.triu_indices(187,1);minimum=float(exact_segments(a[i],b[i],a[j],b[j]).min())-2*field.r
    body_candidates=len(i)
    if history_a:
        value,count=complete_cross(a,b,np.arange(187),np.concatenate(history_a),np.concatenate(history_b),np.concatenate(history_owner),2*field.r+MINIMUM,True);minimum=min(minimum,value-2*field.r);body_candidates+=count
    ss=np.linspace(s0,s1,5);ww=np.stack([field.wraps(s) for s in ss]);ca=ww[:-1].reshape(-1,3);cb=ww[1:].reshape(-1,3)
    all_a=np.concatenate(wrap_a+[ca]);all_b=np.concatenate(wrap_b+[cb]);value,count=complete_cross(a,b,np.arange(187),all_a,all_b,np.zeros(len(all_a),dtype=int),field.r+16e-6+MINIMUM,False);wrap_gap=value-field.r-16e-6
    return {'accepted':bool(minimum>=MINIMUM and wrap_gap>=MINIMUM),'body_minimum_margin_m':minimum,'wrap_minimum_margin_m':wrap_gap,'body_candidate_pairs':body_candidates,'wrap_candidate_pairs':count,'broad_phase_complete':True,'broad_phase':'Variable midpoint radius = clearance threshold + both segment half lengths; full query_ball_point result, then exact finite distances.','wrapper_chord_subdivisions':4},ca,cb

def run():
    start=time.time();data=np.load(UNIT.parent/'carpet-r5-packed-core/receipts/baseline_loops.npz');field=Field(968,data);saved=np.load(UNIT/'arrays/constructor_trace_968.npz');reference=saved['xy'][0]
    initial=np.load(ROOT/'arrays/section_0_candidate.npy');ok,initial_check,world=check_section(field,0,initial,reference);assert ok
    history_a=[];history_b=[];history_owner=[];wrap_a=[];wrap_b=[];stations=[0.];coordinates=[initial];positions=[world];attempts=[];status='running';base_step=field.length/96;step=base_step
    while stations[-1]<field.length-1e-12:
        if time.time()-start>120:status='global_120_second_cap';break
        if len(stations)>=385:status='384_segment_cap';break
        s0=stations[-1];s1=min(s0+step,field.length);trial=coordinates[-1]+field.target(s1)-field.target(s0);_,_,_,k=field.basis(s1)
        fixed,solve=solve_section(trial/field.R,k*field.R,field.rate*field.R,.045,16e-6/field.R,2e-6/field.R,30e-6/field.R,seconds=min(10,120-(time.time()-start)),max_evaluations=160,tolerance=RESIDUAL/field.R)
        xy=fixed*field.R;ok,check,new_world=check_section(field,s1,xy,trial);sweep=None;ca=cb=None
        if ok:sweep,ca,cb=swept(field,s0,s1,positions[-1],new_world,history_a,history_b,history_owner,wrap_a,wrap_b)
        accepted=ok and sweep['accepted'];row={'start_s_m':s0,'end_s_m':s1,'step_um':(s1-s0)*1e6,'accepted':accepted,'solve':solve,'section_check':check,'swept_check':sweep,'maximum_transverse_step_um':float(np.linalg.norm(xy-coordinates[-1],axis=1).max()*1e6)};attempts.append(row)
        print(json.dumps({'attempt':len(attempts),'step_um':row['step_um'],'accepted':accepted,'section_ok':ok,'swept_check':sweep},default=lambda x:x.item()),flush=True)
        if accepted:
            history_a.append(positions[-1]);history_b.append(new_world);history_owner.append(np.arange(187));wrap_a.append(ca);wrap_b.append(cb);stations.append(s1);coordinates.append(xy);positions.append(new_world);step=min(base_step,step*2)
        else:
            step*=.5
            if step<base_step/64-1e-14:status='minimum_step_reached_without_clear_sweep';break
        if len(attempts)%10==0:print(json.dumps({'attempts':len(attempts),'accepted_sections':len(stations),'progress':stations[-1]/field.length,'last_step_um':row['step_um'],'last_acceptance':accepted}),flush=True)
    if stations[-1]>=field.length-1e-12:status='straight_advancement_passed_native_gate_pending'
    np.savez_compressed(ROOT/'arrays/advance_968_trace.npz',stations=stations,xy=coordinates,world=positions)
    result={'scope':'Loop 968 only; unchanged morphology, authoritative float32 section gates and adaptive swept checks. No 16-pass 3D correction, scene or interpolation.','status':status,'straight_advancement_passed':stations[-1]>=field.length-1e-12,'native_curve_qualified':False,'initial_section_check':initial_check,'limits':{'global_seconds':120,'maximum_segments':384,'minimum_step_m':base_step/64,'body_wrap_clearance_target_m':2e-6,'minimum_accepted_clearance_m':MINIMUM,'maximum_distance_residual_m':RESIDUAL,'section_movement_cap_m':30e-6},'accepted_sections':len(stations),'attempts':attempts,'completed_length_fraction':stations[-1]/field.length,'seconds':time.time()-start,'peak_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024}
    (ROOT/'receipts/advance_968_probe.json').write_text(json.dumps(result,indent=2,default=lambda x:x.item()));print(json.dumps({key:value for key,value in result.items() if key!='attempts'},default=lambda x:x.item()))

if __name__=='__main__':run()
