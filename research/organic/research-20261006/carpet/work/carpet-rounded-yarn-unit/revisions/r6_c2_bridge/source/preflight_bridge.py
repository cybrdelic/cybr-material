"""Geometry-only preflight before any contact solve or checkpoint reuse."""
import json,time,resource
import numpy as np
from scipy.interpolate import CubicSpline
from bridge_field import ROOT,UNIT,OriginalField,BridgeField

def curve_measures(p,s,radius):
    v=np.gradient(p,s,axis=1);a=np.gradient(v,s,axis=1);speed=np.linalg.norm(v,axis=-1);curvature=np.linalg.norm(np.cross(v,a),axis=-1)/np.maximum(speed,1e-30)**3
    return {'maximum_sampled_kappa_radius':float(curvature[:,3:-3].max()*radius),'minimum_parameter_speed':float(speed[:,3:-3].min())}

def run():
    started=time.time();data=np.load(UNIT.parent/'carpet-r5-packed-core/receipts/baseline_loops.npz');old=OriginalField(968,data);new=BridgeField(968,data);span=new.half_width;rows=[]
    for join in new.joins:
        samples=np.linspace(join-2*span,join+2*span,4001);p=new.spline(samples);before=old.spline(samples);v=new.spline(samples,1);a=new.spline(samples,2);speed=np.linalg.norm(v,axis=-1);kappa=np.linalg.norm(np.cross(v,a),axis=-1)/speed**3
        bounds=[]
        for boundary in (join-span,join+span):bounds.append({'parameter_m':float(boundary),'position_difference_m':float(np.linalg.norm(new.spline(boundary)-old.spline(boundary))),'derivative_difference':float(np.linalg.norm(new.spline(boundary,1)-old.spline(boundary,1))),'second_derivative_difference_per_m':float(np.linalg.norm(new.spline(boundary,2)-old.spline(boundary,2)))})
        wrappers=[]
        for count in (2001,4001):
            s=np.linspace(join-2*span,join+2*span,count);wrappers.append({'samples':count,**curve_measures(new.wrapper_array(s),s,16e-6)})
        rows.append({'join_parameter_m':float(join),'support_half_width_m':span,'support_full_width_m':2*span,'maximum_same_parameter_guide_displacement_m':float(np.linalg.norm(p-before,axis=1).max()),'speed_range':[float(speed.min()),float(speed.max())],'maximum_sampled_guide_kappa_bundle_radius':float(kappa.max()*new.R),'boundary_matching':bounds,'actual_wrapper_field_curvature_convergence':wrappers,'physical_twist_rate_range_rad_per_m':[float(new.physical_twist_rate(samples).min()),float(new.physical_twist_rate(samples).max())]})
    restart=np.load(ROOT.parent/'r4_bounded_loop/restart/span_26/state.npz');prefix=restart['stations'];checks={f'guide_derivative_{order}_bitwise_equal':bool(np.array_equal(old.spline(prefix,order),new.spline(prefix,order))) for order in range(3)};checks['phase_map_bitwise_equal']=bool(np.array_equal(old.phase0+old.rate*prefix,new.phase0+new.rate*prefix));checks['world_positions_from_same_xy_bitwise_equal']=all(np.array_equal(old.world(float(s),xy),new.world(float(s),xy)) for s,xy in zip(prefix,restart['xy']))
    special=np.array([0.,new.length*.5,new.length]);root=data['center'][968];apex_parameter=float(np.sum(np.linalg.norm(np.diff(new.guide[:1025],axis=0),axis=1)));special[1]=apex_parameter
    checks['roots_and_apex_equal_original_evaluated_guide']=bool(np.array_equal(new.spline(special),old.spline(special)))
    # This is the actual saved pre-correction morphology field transported by
    # the new guide, not a claim about the future solved native fibres.
    trace=np.load(UNIT/'arrays/constructor_trace_968.npz');xy_spline=CubicSpline(trace['s'],trace['xy'],axis=0);body=[]
    for join in new.joins:
        tests=[]
        for count in (1001,2001):
            s=np.linspace(join-2*span,join+2*span,count);c,e1,e2=new.component_arrays(s);xy=xy_spline(s);p=(c[:,None]+xy[:,:,0,None]*e1[:,None]+xy[:,:,1,None]*e2[:,None]).transpose(1,0,2);tests.append({'samples':count,**curve_measures(p,s,field_radius:=.045*new.R)})
        body.append({'join_parameter_m':float(join),'saved_transverse_morphology_field_tests':tests})
    all_s=np.linspace(0,new.length,20001);physical_length=float(np.trapezoid(np.linalg.norm(new.spline(all_s,1),axis=1),all_s));old_length=float(np.trapezoid(np.linalg.norm(old.spline(all_s,1),axis=1),all_s));phase_turns=float(new.rate*new.length/(2*np.pi))
    result={'scope':'One span on each side of both original joins. Geometry-only C2 bridge preflight; no contact solve, source scene or native qualification.','support_choice':'half-width original field.length/96, with the old parameter and phase maps unchanged. No40micrometre alternative was built or tested.','join_checks':rows,'prefix_checks':checks,'saved_body_field_curvature':body,'physical_length_m':physical_length,'physical_length_change_m':physical_length-old_length,'parameter_length_unchanged_m':new.length,'phase_turns_unchanged':phase_turns,'all_prefix_checks_pass':all(checks.values()),'sampled_wrapper_nonfold_pass':all(x['actual_wrapper_field_curvature_convergence'][-1]['maximum_sampled_kappa_radius']<1 for x in rows),'sampled_saved_body_field_nonfold_pass':all(x['saved_transverse_morphology_field_tests'][-1]['maximum_sampled_kappa_radius']<1 for x in body),'limitations':['Sampled true-field curvature convergence is a preflight diagnostic, not a continuous native certificate.','The saved transverse morphology field has not been repaired for contact. Future solved fibres must still pass the actual native curve gates.','The source parameter is intentionally preserved; physical twist rate is phase derivative divided by current guide speed.'],'seconds':time.time()-started,'peak_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024}
    (ROOT/'receipts/bridge_preflight.json').write_text(json.dumps(result,indent=2));np.savez_compressed(ROOT/'arrays/guide_bridge_coefficients.npz',supports=np.array([(lo,hi) for lo,hi,c in new.spline.parts]),coefficients=np.array([c for lo,hi,c in new.spline.parts]));print(json.dumps(result))

if __name__=='__main__':run()
