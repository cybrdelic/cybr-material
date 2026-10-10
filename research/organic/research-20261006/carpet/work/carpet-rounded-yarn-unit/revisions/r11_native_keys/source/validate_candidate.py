"""Read-only actual float32 key validation, including fixed-history contacts.

Partial affected-interfibre gate does not assert same-fibre clearance or qualify
the unchanged history. The complete native prefix gate remains authoritative.
"""
import argparse
from fractions import Fraction as F
import json
import math
import sys
import time
import numpy as np
from scipy.spatial import cKDTree
from stage_keys import ROOT,UNIT,WORK,BridgeField,native,digest,changed_spans,restrict
sys.path.insert(0,str(UNIT/'revisions/r10_physical_length/source'))
from physical_length_window import guide_pieces
from diagnose_retained_span import add,mul,evaluate,rational_report
from acceptance import decision

def elevate(B,degree):
    while B.shape[1]-1<degree:
        n=B.shape[1]-1;C=np.empty((len(B),n+2,3));C[:,0]=B[:,0];C[:,-1]=B[:,-1]
        for i in range(1,n+1):
            alpha=i/(n+1);C[:,i]=alpha*B[:,i-1]+(1-alpha)*B[:,i]
        B=C
    return B

def containment(B,stations,field,changed):
    total=-math.inf;local=-math.inf;worst=None
    for span,(lo,hi) in enumerate(zip(stations[:-1],stations[1:])):
        for fractions,G in guide_pieces(field,lo,hi):
            part=elevate(restrict(B[:,span],float(fractions[0]),float(fractions[-1])),len(G)-1)
            delta=part-G[None];bound=np.linalg.norm(delta,axis=2).max(axis=1)+float(np.float32(field.r))+1e-12
            k=int(np.argmax(bound));v=float(bound[k]);total=max(total,v)
            if changed[span] and bound[0]>local:local=float(bound[0]);worst={'owner':0,'native_span':span,'surface_envelope_upper_m':local}
    return {'whole_prefix_surface_envelope_upper_m':total,'affected_owner0_surface_envelope_upper_m':local,
            'bundle_radius_m':field.R,'whole_prefix_containment_pass':total<=field.R,
            'affected_owner0_containment_pass':local<=field.R,'worst_affected_bound':worst,
            'method':'Piecewise exact C2-guide polynomial subtraction and degree-elevated native Bezier convex-hull bound, with1pm margin.',
            'failed_upper_bound_is_unproven_not_a_penetration_witness':True}

def exact_witness(keys,radius,spans):
    B=native.native_beziers(keys)[spans];t=np.linspace(0,1,129)
    D=3*np.diff(B,axis=1);E=2*np.diff(D,axis=1)
    V=(1-t)[None,:,None]**2*D[:,0,None]+2*(t*(1-t))[None,:,None]*D[:,1,None]+t[None,:,None]**2*D[:,2,None]
    A=(1-t)[None,:,None]*E[:,0,None]+t[None,:,None]*E[:,1,None]
    speed=np.linalg.norm(V,axis=2);ratio=np.linalg.norm(np.cross(V,A),axis=2)/np.maximum(speed,1e-30)**3*radius
    i,j=np.unravel_index(np.argmax(ratio),ratio.shape);span=int(spans[i]);u=F(int(j),128)
    indices=np.clip(np.arange(span-1,span+3),0,len(keys)-1)
    p0,p1,p2,p3=[[F(float(v)) for v in p] for p in keys[indices]]
    power=[p1,[(b-a)/2 for a,b in zip(p0,p2)],
           [a-F(5,2)*b+2*c-d/2 for a,b,c,d in zip(p0,p1,p2,p3)],
           [-a/2+F(3,2)*b-F(3,2)*c+d/2 for a,b,c,d in zip(p0,p1,p2,p3)]]
    D=[[power[1][a],2*power[2][a],3*power[3][a]] for a in range(3)]
    E=[[2*power[2][a],6*power[3][a]] for a in range(3)]
    S=[F(0)];T=[F(0)]
    for a in range(3):S=add(S,mul(D[a],D[a]))
    for a,b in ((1,2),(2,0),(0,1)):
        cross=add(mul(D[a],E[b]),[-v for v in mul(D[b],E[a])]);T=add(T,mul(cross,cross))
    r=F(float(radius));H=add(mul(mul(S,S),S),[-r*r*v for v in T]);su,tu,hu=evaluate(S,u),evaluate(T,u),evaluate(H,u)
    return {'owner':0,'native_span':span,'exact_parameter':str(u),'exact_speed_squared':str(su),
            'exact_nonfold_polynomial':str(hu),'exact_curvature_radius_squared':str(r*r*tu/su**3) if su>0 else None,
            'strict_nonfold_violation_proved':su>0 and hu<0,'kappa_radius':math.sqrt(float(r*r*tu/su**3)) if su>0 else None,
            'exact_float32_native_power_coefficients':rational_report(power),'exact_radius':str(r)}

def affected_leaf_mask(owners,spans,changed):
    # Wrapper native spans can exceed the body's key count. Index the body's
    # mask only after selecting the changed owner, never via eager boolean &.
    result=np.zeros(len(owners),dtype=bool);selected=owners==0
    result[selected]=changed[spans[selected]]
    return result

def contact_gate(body,wraps,radius,wrap_radius,changed):
    config=native.GateConfig(wall_seconds=50,max_leaves=250000,max_candidate_pairs=2000000,memory_limit_mib=900)
    budget=native.Budget(config);groups=[];count=0
    for owner in range(193):
        keys=body[owner] if owner<187 else wraps[owner-187]
        pad=256*native.EPS*max(1e-6,float(np.abs(keys).max()))*19
        leaves,_=native.owner_leaves(native.native_beziers(keys),owner,math.inf,config,budget,pad,config.max_leaves-count)
        count+=len(leaves['a']);groups.append(leaves)
    leaves={k:np.concatenate([g[k] for g in groups]) for k in groups[0]};del groups
    a,b,error,owners=(leaves[k] for k in ('a','b','error','owner'))
    radii=np.where(owners<187,radius,wrap_radius);mid=(a+b)*.5;extent=np.linalg.norm(b-a,axis=1)*.5+error+radii
    affected=affected_leaf_mask(owners,leaves['native_span'],changed);ids=np.flatnonzero(affected)
    tree=cKDTree(mid);total=0;tested=0;bad=0;reserve_bad=0;minimum=math.inf;failures=[]
    pad=config.distance_margin_m+512*native.EPS*max(1e-6,float(np.abs(mid).max()))
    for first in range(0,len(ids),16):
        budget.check();sel=ids[first:first+16];horizon=extent[sel]+extent.max()+native.INTERFIBRE_GAP_M+pad
        sizes=tree.query_ball_point(mid[sel],horizon,workers=1,return_length=True)
        if sizes.max()>config.max_neighbors_per_query:raise native.Unproven('neighbor_budget_exceeded')
        if total+int(sizes.sum())>config.max_candidate_pairs:raise native.Unproven('affected_contact_candidate_budget_exceeded',count=total)
        lists=tree.query_ball_point(mid[sel],horizon,workers=1)
        ii=np.repeat(sel,[len(x) for x in lists]);jj=np.concatenate(lists).astype(int);total+=len(jj)
        keep=owners[ii]!=owners[jj];ii=ii[keep];jj=jj[keep]
        for offset in range(0,len(ii),8192):
            budget.check();i=ii[offset:offset+8192];j=jj[offset:offset+8192]
            separation=np.maximum(np.maximum(np.minimum(a[j],b[j])-np.maximum(a[i],b[i]),np.minimum(a[i],b[i])-np.maximum(a[j],b[j])),0)
            far=np.max(separation,axis=1)>radii[i]+radii[j]+error[i]+error[j]+native.INTERFIBRE_GAP_M+pad
            i=i[~far];j=j[~far]
            if not len(i):continue
            lower=native.projection_distance_lower(a[i],b[i],a[j],b[j])-error[i]-error[j]-radii[i]-radii[j]-pad
            tested+=len(i);minimum=min(minimum,float(lower.min()));bad+=int((lower<=0).sum());reserve_bad+=int((lower<native.INTERFIBRE_GAP_M).sum())
            for k in np.flatnonzero(lower<=0):
                if len(failures)>=16:break
                failures.append({'first':native._leaf_label(leaves,int(i[k])),'second':native._leaf_label(leaves,int(j[k])),
                                 'continuous_surface_clearance_lower_m':float(lower[k])})
    return {'complete_affected_interfibre_search':True,'leaf_count_all_193_owners':count,'affected_owner0_leaf_count':len(ids),
            'sphere_candidates':total,'projection_pairs_tested':tested,'nonpositive_bound_pairs':bad,'construction_reserve_unproven_pairs':reserve_bad,
            'physical_clearance_target_m':0,'construction_reserve_m':native.INTERFIBRE_GAP_M,
            'surface_clearance_lower_m':min(native.INTERFIBRE_GAP_M,minimum),
            'positive_affected_physical_clearance_pass':bad==0,'affected_construction_reserve_pass':reserve_bad==0,
            'includes_all_fixed_history_and_all_six_wrappers':True,'same_owner_nonlocal_clearance_qualified':False,
            'unchanged_owner_pair_clearance_qualified':False,'failure_examples':failures,'elapsed_seconds':time.monotonic()-budget.start}

def run(name):
    start=time.monotonic();path=ROOT/'arrays'/name;values=np.load(path)
    body,old,wraps=values['body'],values['original_body'],values['wraps'];stations=values['stations'];refs=values['references']
    radius=float(values['body_radius']);wrap_radius=float(values['wrap_radius'])
    field=BridgeField(968,np.load(WORK/'carpet-r5-packed-core/receipts/baseline_loops.npz'))
    B,O,changed,mapping=changed_spans(body,old,np.arange(49));allowed=np.zeros((187,49),bool);allowed[0,35:39]=True
    history={'other186_fibres_bitwise_preserved':np.array_equal(body[1:],old[1:]),
             'all_keys_outside_local_mask_bitwise_preserved':np.array_equal(body[~allowed],old[~allowed]),
             'all_native_spans_outside33_through39_exactly_preserved':np.array_equal(B[:,np.r_[0:33,40:48]],O[:,np.r_[0:33,40:48]]),
             'changed_native_spans':np.flatnonzero(changed).tolist(),'native_curves_remain_uniform_Catmull_Rom_C1':True,
             'guide_C2_polynomial_unchanged':True,'roots_and_prefix_endpoints_bitwise_preserved':np.array_equal(body[:,[0,-1]],old[:,[0,-1]])}
    motion=[];axial=[];original_motion=[]
    for k,s in enumerate(stations):
        center,e1,e2,_=field.basis(s);tangent=np.cross(e1,e2);reference=field.world(s,refs[k])
        motion.append(np.linalg.norm(body[:,k].astype(float)-reference,axis=1));original_motion.append(np.linalg.norm(old[:,k].astype(float)-reference,axis=1))
        axial.append(np.abs((body[:,k].astype(float)-center)@tangent))
    motion=np.asarray(motion).T;axial=np.asarray(axial).T;original_motion=np.asarray(original_motion).T
    disks={'original_radius_m':30e-6,'section_plane_rounding_tolerance_m':1e-9,
           'affected_maximum_distance_from_original_reference_m':float(motion[allowed].max()),
           'affected_maximum_axial_plane_residual_m':float(axial[allowed].max()),
           'affected_transverse_disks_pass':bool(motion[allowed].max()<=30e-6+1e-12 and axial[allowed].max()<=1e-9),
           'whole_prefix_maximum_reference_distance_m':float(motion.max()),'baseline_maximum_reference_distance_m':float(original_motion.max()),
           'whole_prefix_original30um_disks_pass':bool(motion.max()<=30e-6+1e-12 and axial.max()<=1e-9),
           'inherited_disk_overshoots_are_not_silently_accepted':True}
    print(json.dumps({'stage':'independent_containment_and_contact','candidate':name}),flush=True)
    envelope=containment(B,stations,field,changed)
    try:contacts=contact_gate(body,wraps,radius,wrap_radius,changed)
    except native.Unproven as exc:contacts={'complete_affected_interfibre_search':False,'positive_affected_physical_clearance_pass':False,
                                          'affected_construction_reserve_pass':False,'code':exc.code,**exc.details}
    print(json.dumps({'stage':'complete_native_prefix_curvature_and_clearance','candidate':name}),flush=True)
    native_report=native.qualify_npz(path,identity_world=True,config=native.GateConfig(wall_seconds=30,max_leaves=250000,memory_limit_mib=900))
    witness=exact_witness(body[0],radius,np.flatnonzero(changed))
    gates={'strict_nonfold_pass':native_report['strict_nonfold_pass'],
           'positive_physical_clearance_pass':native_report['physical_nonpenetration_pass'],
           'same_fibre_self_clearance_pass':native_report['same_fibre_self_clearance_pass'],
           'continuous_containment_pass':envelope['whole_prefix_containment_pass'],
           'original30um_transverse_disks_pass':disks['whole_prefix_original30um_disks_pass'],
           'roots_preserved':history['roots_and_prefix_endpoints_bitwise_preserved'],
           'unchanged_history_preserved':all(history[k] for k in ('other186_fibres_bitwise_preserved','all_keys_outside_local_mask_bitwise_preserved','all_native_spans_outside33_through39_exactly_preserved')),
           'candidate_hash_matches_gate':native_report.get('input_npz_sha256')==digest(path),
           'construction_reserve_pass':native_report['construction_gap_preserved_pass'],
           'production_loop_complete':False,'actual_scene_contract_verified':False,'three_loop_qualified':False,'matched_visual_review_pass':False}
    report={'candidate_filename':name,'candidate_sha256':digest(path),'actual_float32_required':body.dtype==np.dtype('f4'),
            'scope':'Partial loop968; full-prefix native gate plus separate changed-owner contact/containment/history constraints.',
            'history':history,'disks':disks,'containment':envelope,'affected_interfibre_contacts':contacts,
            'native_prefix_gate':native_report,'exact_new_candidate_curvature_witness':witness,
            'gates':gates,'decision':decision(gates),'elapsed_seconds':time.monotonic()-start,
            'positive_physical_clearance_separate_from_1_95um_reserve':True,'no_render_or_scene_export':True}
    output=ROOT/'receipts'/(path.stem+'_validation.json');output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'candidate':name,'decision':report['decision'],'native_status':native_report['status'],'contact_result':contacts.get('positive_affected_physical_clearance_pass'),
                      'elapsed_seconds':report['elapsed_seconds']}),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('candidate');run(parser.parse_args().candidate)
