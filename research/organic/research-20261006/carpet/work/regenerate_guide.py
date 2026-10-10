"""Guide-only deterministic regeneration. Never claims missing NPZ recovery."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np

ROOT = Path(__file__).resolve().parent
UNIT = ROOT/'carpet-rounded-yarn-unit'
sys.path.insert(0,str(ROOT/'carpet-crown-feasibility/source'))
from measure_guides_v2 import guide

def digest(path):
    return hashlib.file_digest(Path(path).open('rb'),'sha256').hexdigest()

def run():
    baseline=ROOT/'carpet-r5-packed-core/receipts/baseline_loops.npz'
    data=np.load(baseline)
    directory=UNIT/'arrays'
    directory.mkdir(exist_ok=True)
    rows=[]
    for index in (968,955,356):
        c=data['center'][index]
        dense,description=guide(c,float(data['radius'][index].mean()))
        assert dense.shape==(2049,3) and dense.dtype==np.float64
        assert np.array_equal(dense[[0,-1]],c[[0,-1]])
        assert np.array_equal(dense[1024],c[24])
        path=directory/f'rounded_unit_{index}.npz'
        if path.exists():
            assert np.array_equal(np.load(path)['guide'],dense), 'refuse overwrite of other guide'
        else:
            np.savez_compressed(path,guide=dense)
        rows.append({'loop':index,'guide_array_sha256':hashlib.sha256(dense.tobytes()).hexdigest(),
                     'npz_sha256':digest(path),'guide_length_m':float(np.linalg.norm(np.diff(dense,axis=0),axis=1).sum()),
                     'roots_exact':True,'apex_exact':True,**description})
    sys.path.insert(0,str(UNIT/'revisions/r6_c2_bridge/source'))
    from bridge_field import BridgeField
    field=BridgeField(968,data)
    recorded=json.loads((UNIT/'revisions/r6_c2_bridge/receipts/bridge_preflight.json').read_text())
    saved=np.load(UNIT/'revisions/r6_c2_bridge/arrays/guide_bridge_coefficients.npz')
    supports=np.array([[lo,hi] for lo,hi,_ in field.spline.parts])
    coefficients=np.array([c for _,_,c in field.spline.parts])
    length_error=abs(field.length-recorded['parameter_length_unchanged_m'])
    support_error=float(np.max(abs(supports-saved['supports'])))
    coefficient_error=float(np.max(abs(coefficients-saved['coefficients'])))
    # Different NumPy/libm builds may differ in float64 by a few ulps. These
    # sub-picometre checks precede a direct check of stored accepted coordinates.
    assert length_error<=1e-14 and support_error<=1e-14 and coefficient_error<=1e-12
    state=np.load(UNIT/'revisions/r9_physical_resume/arrays/loop_968_partial.npz')
    predicted=np.stack([field.world(s,xy) for s,xy in zip(state['stations'],state['xy'])])
    stored=state['world']
    error=float(np.max(np.linalg.norm(predicted-stored,axis=-1)))
    rounding=float(np.max(np.linalg.norm(predicted.astype('f4').astype('f8')-stored,axis=-1)))
    # Historical transverse states precede float32 round-trip projection; they
    # can differ by one coordinate ulp. This explicit tolerance is 1 nanometre,
    # ~0.0035% of diameter and far below the 50nm historical residual allowance.
    assert error<=1e-9 and rounding<=1e-9
    report={'status':'deterministic_guide_regenerated_and_restart_checks_passed',
            'exact_missing_archive_recovery':False,'guide_only_npzs':True,
            'constructor_trace_recovered':False,'baseline_sha256':digest(baseline),
            'guide_source_sha256':digest(ROOT/'carpet-crown-feasibility/source/measure_guides_v2.py'),
            'rows':rows,'restart_checks':{'loop':968,'accepted_spans':len(state['stations'])-1,
                'length_error_m':float(length_error),'bridge_support_max_error_m':support_error,
                'bridge_coefficient_max_error_m':coefficient_error,
                'accepted_world_max_error_m':error,'accepted_float32_world_max_error_m':rounding,
                'accepted_world_tolerance_m':1e-9,'checked_points':int(stored.shape[0]*stored.shape[1])},
            'limitations':['955 and356 guides lack recorded accepted-prefix validation; no production fibres created.',
                           'Original missing NPZ container hashes are not recovered; original body/wrap arrays and constructor trace are absent.']}
    output=ROOT.parent/'receipts/guide_regeneration.json'
    output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':run()
