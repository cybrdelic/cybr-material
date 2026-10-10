"""Only saved section 0 and worst section 54 of the unchanged loop 968."""
from pathlib import Path
import json,time,resource,hashlib
import numpy as np
from section_feasibility import validate_gradients,solve_section
ROOT=Path(__file__).resolve().parents[1];UNIT=ROOT.parents[1]

def run():
    start=time.time();gradient=validate_gradients();print(json.dumps({'gradient_validation':gradient}),flush=True)
    saved=np.load(UNIT/'arrays/constructor_trace_968.npz');data=np.load(UNIT.parent/'carpet-r5-packed-core/receipts/baseline_loops.npz');R=float(data['radius'][968].mean());rows=[]
    for section in (0,54):
        original=saved['xy'][section];reference=original/R;k=saved['curvature'][section]*R;w=float(saved['rate'])*R
        fixed,report=solve_section(reference,k,w,.045,16e-6/R,2e-6/R,30e-6/R)
        report.update({'section':section,'radius_m':R,'input_sha256':hashlib.sha256(original.tobytes()).hexdigest(),'maximum_violation_um':report['final']['maximum_constraint_violation_normalized']*R*1e6,'maximum_movement_um':report['final']['maximum_movement_normalized']*R*1e6,'minimum_body_gap_um':report['final']['minimum_body_metric_gap_normalized']*R*1e6,'minimum_wrap_gap_um':report['final']['minimum_wrap_metric_gap_normalized']*R*1e6,'input_is_saved_failed_section':'Section 0 is the saved seed/first bad section after the former fixed-pass projection; section 54 is the saved worst section.'})
        np.save(ROOT/f'arrays/section_{section}_candidate.npy',fixed*R);rows.append(report);print(json.dumps(report),flush=True)
        result={'scope':'Two unchanged failed section inputs only. Morphology fixed. No 3D correction, advancement, interpolation, scene or manufacture claim.','gradient_validation':gradient,'rows':rows,'terminal':len(rows)==2,'both_pass':len(rows)==2 and all(row['accepted'] for row in rows),'seconds':time.time()-start,'peak_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024}
        (ROOT/'receipts/section_feasibility_probe.json').write_text(json.dumps(result,indent=2))

if __name__=='__main__':run()
