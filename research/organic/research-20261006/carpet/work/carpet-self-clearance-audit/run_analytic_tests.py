"""Small analytic source-qualification tests; no scene/render/constructor work."""
from pathlib import Path
import hashlib
import json
import numpy as np
from self_clearance_reference import (Leaves,assess_self,join,line,circular_arc,
                                      cubic_leaves,projection_distance_lower)

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent


def run():
    results={}
    r=15e-6
    straight=line([0,0,0],[.004,0,0],144)
    results['straight_4mm']=assess_self(straight,0,r)
    assert results['straight_4mm']['pass']
    assert results['straight_4mm']['distant_pair_count']==0

    # A circle has reach exactly R; closed periodic seam is one parameter point.
    R=15.4e-6
    circle=circular_arc([0,0,0],R,0,2*np.pi,128)
    results['circle_radius_15p4um_tube_15um']=assess_self(circle,1/R,r,closed=True)
    assert results['circle_radius_15p4um_tube_15um']['pass']
    results['circle_focal_equality']=assess_self(circle,1/R,R,closed=True)
    assert not results['circle_focal_equality']['pass']

    # C1,1 line--semicircle--line hairpin. Its two arms have only 0.8 um
    # surface clearance, legitimately below the 1.95 um inter-fibre gap.
    H=120e-6
    hairpin=join(line([-H,0,0],[0,0,0],32),
                 circular_arc([0,R,0],R,-np.pi/2,np.pi/2,64),
                 line([0,2*R,0],[-H,2*R,0],32))
    results['hairpin_physical_gap_0p8um']=assess_self(hairpin,1/R,r)
    assert results['hairpin_physical_gap_0p8um']['pass']
    assert 0 < results['hairpin_physical_gap_0p8um']['minimum_tested_surface_clearance_lower_m'] < 1.95e-6

    # Near-return: wide gentle turn but the two open endpoints have separation
    # 20 um and opposite horizontal tangents. Their normal discs intersect at
    # (0,0,0), 10 um from each center, despite K*r=0.15 everywhere.
    R=100e-6;H=1.84e-3;d=20e-6
    lower_control=np.array([[0,-d/2,0],[H/3,-d/2,0],
                            [2*H/3,-R,0],[H,-R,0]])
    lower=cubic_leaves(lower_control,64)
    upper_control=lower_control.copy();upper_control[:,1]*=-1
    upper=cubic_leaves(upper_control[::-1],64)
    near=join(lower,circular_arc([H,0,0],R,-np.pi/2,np.pi/2,64),upper)
    # x'=H on each arm; |y''| <=6*(R-d/2), so kappa<=6*(R-d/2)/H^2.
    K=max(1/R,6*(R-d/2)/H**2)
    results['near_return_actual_endpoint_disc_overlap']=assess_self(near,K,r)
    results['near_return_actual_endpoint_disc_overlap']['analytic_witness']={
        'point_m':[0,0,0], 'normal_offsets_m':[d/2,d/2],
        'tube_radius_m':r, 'normal_dot_tangent':0.,
        'total_length_control_polygon_upper_m':float(near.arc_upper.sum())}
    assert not results['near_return_actual_endpoint_disc_overlap']['pass']
    assert results['near_return_actual_endpoint_disc_overlap']['kappa_radius']<1

    # Deliberately giant unresolved leaf: never discard intra-leaf contacts.
    results['long_leaf_refinement_required']=assess_self(
        Leaves(np.array([[0.,0,0]]),np.array([[0.,1e-6,0]]),
               np.array([1e-6]),np.array([1e-3])),1/(100e-6),r)
    assert results['long_leaf_refinement_required']['status']=='unproven_refine_long_leaves'

    # Recheck authoritative r4 lower bounds against analytic crossing distances.
    crossing=[]
    for angle in (2e-7,1e-6,1e-5,np.pi/2):
        for scale in (1e-9,1e-6,1.,1e6):
            a=np.zeros(3);b=np.array([1.,0,0]);v=np.array([np.cos(angle),np.sin(angle),0])
            c=.37*b-.61*v;dd=c+v
            lower=float(projection_distance_lower(a*scale,b*scale,c*scale,dd*scale))
            assert lower<=1e-14*scale
            crossing.append({'angle_rad':angle,'scale':scale,'lower_bound':lower,
                             'analytic_distance':0.})
    results['projection_gate_crossings']=crossing
    source_paths=[HERE/'self_clearance_reference.py',HERE/'run_analytic_tests.py',
                  ROOT/'carpet-rounded-yarn-unit/revisions/r4_bounded_loop/source/certified_segments.py']
    receipt={'scope':'Analytic source qualification only; no production fibre arrays, scenes, rendering, or source tuning.',
             'all_assertions_passed':True,'results':results,
             'limits':['Regularity, curvature and continuous leaf bounds are assumptions of assess_self.',
                       'These synthetic cases provide analytic bounds; no production native field is qualified.',
                       'Floating implementation is an engineering check, not a formal interval arithmetic proof.',
                       'A failed sufficient condition means unproven unless a separate collision witness is supplied.'],
             'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}}
    (HERE/'analytic_tests_receipt.json').write_text(json.dumps(receipt,indent=2))
    print(json.dumps({'all_assertions_passed':True,
                      'cases':{k:v['status'] for k,v in results.items() if isinstance(v,dict)},
                      'crossing_regressions':len(crossing)},indent=2))


if __name__=='__main__':run()
