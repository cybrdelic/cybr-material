"""Small algebra, actual converted-span and acceptance regressions; no solver."""
from fractions import Fraction as F
import json
import time
import numpy as np
from stage_keys import ROOT,UNIT,native,digest
from parameterization_probe import inverse_native,hermite_bezier
from validate_candidate import affected_leaf_mask
from acceptance import LOCAL_REQUIRED,decision

def run():
    start=time.monotonic();groups=[]
    # Arbitrary cubic controls, not a line: exact inverse native control algebra.
    B=[[F(1,8),F(2,7),F(-3,11)],[F(2,9),F(3,8),F(1,5)],
       [F(-1,6),F(4,9),F(2,3)],[F(1,2),F(-2,5),F(3,7)]]
    p1,p2=B[0],B[3]
    p0=[b-6*(c-a) for a,b,c in zip(p1,p2,B[1])]
    p3=[a+6*(b-c) for a,b,c in zip(p1,p2,B[2])]
    recovered=[p1,[a+(b-c)/6 for a,b,c in zip(p1,p2,p0)],
               [a-(b-c)/6 for a,b,c in zip(p2,p3,p1)],p2]
    assert recovered==B;groups.append('exact arbitrary cubic inverse native conversion')
    # Hermite physical derivatives are multiplied by the span length once.
    B=np.array([[0.,0,0],[1,2,0],[2,3,1],[4,5,2]])*1e-5
    h=7e-6;d1=3*(B[1]-B[0])/h;d2=3*(B[3]-B[2])/h
    assert np.allclose(hermite_bezier(B[0],B[3],d1,d2,h),B,rtol=1e-15,atol=1e-20)
    assert np.allclose(native.native_beziers(inverse_native(B))[1],B,rtol=1e-15,atol=1e-20)
    groups.append('physical derivative scaling and actual native conversion')
    # S_v=S_u a^2,T_v=T_u a^6,H_v=H_u a^6. Strict sign is invariant.
    S=F(7,9);T=F(8,5);r=F(3,4);H=S**3-r*r*T;assert H<0
    for a in (F(1,16),F(3,2),F(-7,3)):
        assert (S*a*a)**3-r*r*(T*a**6)==H*a**6<0
    groups.append('exact regular reparameterization nonfold sign invariance')
    values=np.load(ROOT/'arrays/physical_tangent_conversion.npz');keys=values['quadratic_consistent_physical_tangents_four_native_keys']
    assert keys.dtype==np.dtype('f4')
    actual=native.native_beziers(keys)[1:2]
    assert np.array_equal(actual[0],values['quadratic_consistent_physical_tangents_actual_bezier'])
    receipt=json.loads((ROOT/'receipts/physical_tangent_conversion.json').read_text())
    radius=float(F(995775,68719476736));config=native.GateConfig(wall_seconds=3,max_leaves=10000)
    pad=256*native.EPS*max(1e-6,float(np.abs(keys).max()))*19
    proof,_=native.owner_curvature(actual,radius,config,native.Budget(config),pad,0)
    assert proof['kappa_radius_upper']<.151 and proof['native_parameter_speed_lower_m']>0
    assert not receipt['methods'][0]['original_transverse_disks_pass']
    groups.append('actual float32 physical tangent span certified but disk planes rejected')
    candidate=np.load(ROOT/'arrays/native_tangent_candidate_968.npz');body=candidate['body'];old=candidate['original_body']
    assert np.array_equal(body[1:],old[1:]) and np.array_equal(body[:,[0,-1]],old[:,[0,-1]])
    B=native.native_beziers(body[0]);O=native.native_beziers(old[0])
    assert np.array_equal(B[np.r_[0:33,40:48]],O[np.r_[0:33,40:48]])
    assert np.max(np.abs(3*(B[:-1,3]-B[:-1,2])-3*(B[1:,1]-B[1:,0])))<2e-18
    groups.append('unchanged other fibres, roots, outside native cubics and C1 joins')
    full={key:True for key in LOCAL_REQUIRED}
    full.update(construction_reserve_pass=True,production_loop_complete=True,actual_scene_contract_verified=True,
                three_loop_qualified=True,matched_visual_review_pass=True)
    assert decision(full)['promotion_allowed']
    for key in LOCAL_REQUIRED:
        missing=full.copy();missing[key]=False
        assert not decision(missing)['production_export_allowed'] and not decision(missing)['promotion_allowed']
    for key in ('construction_reserve_pass','production_loop_complete','actual_scene_contract_verified'):
        missing=full.copy();missing[key]=False;assert not decision(missing)['production_export_allowed']
    missing=full.copy();missing['matched_visual_review_pass']=False
    assert not decision(missing)['promotion_allowed']
    assert not decision({'strict_nonfold_pass':True})['local_repair_accepted']
    groups.append('export and promotion fail closed on native curvature, every geometry gate and visual review')
    mask=np.zeros(48,bool);mask[36]=True
    assert affected_leaf_mask(np.array([0,1,187,192]),np.array([36,36,900,1023]),mask).tolist()==[True,False,False,False]
    groups.append('wrapper span indices cannot index body affected mask')
    output={'groups_passed':len(groups),'groups':groups,'elapsed_seconds':time.monotonic()-start,
            'local_kappa_radius_upper':proof['kappa_radius_upper'],
            'source_hashes':{p.name:digest(p) for p in Path(__file__).parent.glob('*.py')}}
    (ROOT/'receipts/native_key_tests.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output),flush=True)

if __name__=='__main__':
    from pathlib import Path
    run()
