"""Physical-length tangents and explicit inverse native control conversion.

The Cycles exporter has no knot/tangent channel: it copies float32 keys and the
kernel constructs uniform Catmull-Rom controls. Reparameterizing a polynomial
does not change curvature. Altering these tangents changes the centreline.
"""
import json
import math
import time
import numpy as np
from scipy.integrate import quad
from stage_keys import ROOT,UNIT,WORK,BridgeField,native,digest

def hermite_bezier(p1,p2,d1,d2,length):
    return np.stack((p1,p1+length*d1/3,p2-length*d2/3,p2))

def inverse_native(B):
    p1,p2=B[0],B[3]
    return np.stack((p2-6*(B[1]-p1),p1,p2,p1+6*(p2-B[2])))

def run():
    start=time.monotonic();frozen=UNIT/'revisions/r9_physical_resume/arrays/loop_968_partial.npz'
    state=np.load(frozen);old=state['world'].astype('f4').transpose(1,0,2);p=old[0].astype(float)
    field=BridgeField(968,np.load(WORK/'carpet-r5-packed-core/receipts/baseline_loops.npz'))
    s=state['stations'];arc=[];error=[]
    for lo,hi in zip(s[34:39],s[35:40]):
        value,e=quad(lambda x:float(np.linalg.norm(field.spline(x,1))),lo,hi,epsabs=1e-15,epsrel=1e-11)
        arc.append(value);error.append(e)
    h=np.asarray(arc);left=(p[36]-p[35])/h[1];middle=(p[37]-p[36])/h[2];right=(p[38]-p[37])/h[3]
    d1=(h[2]*left+h[1]*middle)/(h[1]+h[2]);d2=(h[3]*middle+h[2]*right)/(h[2]+h[3])
    sec1=(p[37]-p[35])/(h[1]+h[2]);sec2=(p[38]-p[36])/(h[2]+h[3])
    methods={'quadratic_consistent_physical_tangents':(d1,d2),
             'central_physical_secants':(sec1,sec2),'local_chord_tangents':(middle,middle)}
    rows=[];arrays={};native_B=native.native_beziers(p)[36]
    for name,(a,b) in methods.items():
        target=hermite_bezier(p[36],p[37],a,b,h[2]);new=inverse_native(target).astype('f4')
        actual=native.native_beziers(new)[1:2]
        config=native.GateConfig(wall_seconds=10,max_leaves=100000)
        pad=256*native.EPS*max(1e-6,float(np.abs(new).max()))*19
        try:
            proof,_=native.owner_curvature(actual,float(np.float32(field.r)),config,native.Budget(config),pad,0)
            proof={'pass':True,**proof}
        except native.Unproven as exc:proof={'pass':False,'code':exc.code,**exc.details}
        delta=new.astype(float)-p[35:39]
        plane=[];disk=[]
        for k in range(4):
            center,e1,e2,_=field.basis(s[35+k]);tangent=np.cross(e1,e2)
            reference=field.world(s[35+k],state['references'][35+k,0:1])[0]
            plane.append(float((new[k].astype(float)-center)@tangent))
            disk.append(float(np.linalg.norm(new[k].astype(float)-reference)))
        rows.append({'method':name,'actual_float32_central_native_span_curvature':proof,
                     'control_conversion_max_roundoff_m':float(np.abs(actual[0]-target).max()),
                     'stored_inner_endpoints_bitwise_preserved':np.array_equal(new[1:3],old[0,36:38]),
                     'required_outer_key35_and38_world_shift_m':[float(np.linalg.norm(delta[0])),float(np.linalg.norm(delta[3]))],
                     'original_section_plane_axial_residuals_m':plane,'distance_from_original_disk_centres_m':disk,
                     'original_transverse_disk_plane_tolerance_m':1e-9,
                     'original_transverse_disks_pass':max(abs(v) for v in plane)<=1e-9 and max(disk)<=30e-6+1e-12,
                     'native_controls_equal_existing_curve':np.array_equal(actual[0],native_B),
                     'curve_change_control_hull_upper_m':float(np.linalg.norm(actual[0]-native_B,axis=1).max())})
        arrays[name+'_four_native_keys']=new;arrays[name+'_target_bezier']=target;arrays[name+'_actual_bezier']=actual[0]
        print(json.dumps(rows[-1]),flush=True)
    # Preserve the physical tangent candidate for all-prefix revalidation.
    # The central span is certified; original transverse planes still reject it.
    body=old.copy();body[0,35:39]=arrays['quadratic_consistent_physical_tangents_four_native_keys']
    wraps=np.load(UNIT/'revisions/r10_physical_length/arrays/trial_968_prefix_native.npz')
    np.savez_compressed(ROOT/'arrays/physical_tangent_conversion.npz',**arrays,physical_lengths=h,
                        original_four_keys=old[0,35:39],original_native_bezier=native_B)
    np.savez_compressed(ROOT/'arrays/native_tangent_candidate_968.npz',body=body,original_body=old,
                        wraps=wraps['wraps'],body_radius=np.float32(field.r),wrap_radius=wraps['wrap_radius'],
                        stations=s,references=state['references'],key_ids=np.arange(49),
                        changed_owners=np.array([0]),allowed_key_indices=np.array([35,36,37,38]))
    report={'scope':'One native span36 of owner0, unchanged centre keys36/37; separate tangent/interpolation conversion experiment.',
            'source_frozen_sha256':digest(frozen),'inspected_exporter':{
                'local_builder':'carpet-r5-packed-core/source/build_packed_panel.py:68..75',
                'local_builder_sha256':digest(WORK/'carpet-r5-packed-core/source/build_packed_panel.py'),
                'pinned_exporter_evidence':'carpet-self-clearance-audit/cycles_representation_receipt.json',
                'pinned_blender_commit':'32f5fdce0a0a571a39198edf8d7b5d719bf59f2c',
                'native_point_position_export':'Raw point positions; no authored stations, derivative or tangent attributes copied to Cycles hair.',
                'kernel_controls':'B=[P1,P1+(P2-P0)/6,P2-(P3-P1)/6,P2]',
                'uniform_knot_spacing_assumption_confirmed':True},
            'station35_through38_steps_m':np.diff(s[35:39]).tolist(),
            'physical_lengths34_through39_m':h.tolist(),'arc_quadrature_error_m':float(sum(error)),
            'largest_to_smallest_step35_through38':float(h[1:4].max()/h[1:4].min()),
            'inverse_native_conversion':'P0=P2-6(B1-P1);P3=P1+6(P2-B2). Outer keys cannot stay fixed when these tangents change.',
            'fixed_key_derivative_conflict':{'required_native_start_tangent':(.5*(p[37]-p[35])).tolist(),
                'physical_quadratic_start_tangent':(h[2]*d1).tolist(),
                'start_tangent_difference_m':float(np.linalg.norm(h[2]*d1-.5*(p[37]-p[35]))),
                'required_native_end_tangent':(.5*(p[38]-p[36])).tolist(),
                'physical_quadratic_end_tangent':(h[2]*d2).tolist()},
            'methods':rows,'all_prefix_candidate_method':'quadratic_consistent_physical_tangents',
            'isolated_local_span_pass_is_not_full_prefix_qualification':True,
            'pure_reparameterization_changes_curvature':False,'new_tangents_change_physical_centreline':True,
            'original_keys_outside35_through38_and_other_fibres_preserved':True,
            'C1_native_parameter_tangents_at_joins_preserved_by_uniform_basis':True,
            'C2_seams_not_claimed':'The native hair contract is piecewise C2 and C1 at keys; second derivative continuity is not automatic. The separate guide C2 polynomial is untouched.',
            'accepted':False,'selected_r5_promoted':False,'elapsed_seconds':time.monotonic()-start}
    (ROOT/'receipts/physical_tangent_conversion.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'physical_tangent_probe_complete':True,'local_span_passes':[x['method'] for x in rows if x['actual_float32_central_native_span_curvature']['pass']]}),flush=True)

if __name__=='__main__':run()
