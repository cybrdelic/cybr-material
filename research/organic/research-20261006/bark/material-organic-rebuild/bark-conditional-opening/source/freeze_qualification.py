"""Freeze the completed numerical comparisons; no image or morphology branch."""
import hashlib,json,py_compile
from pathlib import Path
import numpy as np
from opening_model import ROOT
from checkpoint_io import read_checkpoint

labels=['central_r2','radial_half_r2','increment_half','penalty_low','penalty_high','quadrature5']
states=[]
for label in labels:
    path=ROOT/'state'/label/'restart.npz';a,m=read_checkpoint(path)
    side=json.loads((path.parent/'restart.json').read_text())
    sha=hashlib.sha256(path.read_bytes()).hexdigest();assert side['state_sha256']==sha
    assert m['status']=='completed' and m['trace'][-1]['load']==1.
    trace=m['trace'];d=np.array([r['new_dissipation_from_initial_J'] for r in trace]);assert np.all(np.diff(d)>=-1e-13)
    assert np.all(a['maximum']>=a['initial_maximum']) and np.all(a['damage']>=0) and np.all(a['damage']<=1)
    assert np.all(a['initial_maximum'][~a['starter']]==0) and np.all(a['initial_maximum'][a['starter']]>0)
    assert trace[0]['stored_energy_J']==0. and trace[0]['new_dissipation_from_initial_J']==0.
    for row in trace:
        assert row['minimum_physical_J']>0
        assert row['maximum_free_force_N']<5e-5
        assert np.all(np.array(row['maximum_tangential_strain_by_cohort'])<np.array([.0597,.045,.0444])+1e-12)
        assert abs(row['reference_mass_kg']-m['reference_mass_kg'])<1e-15
    process_records=[m]+[json.loads(p.read_text()) for p in (ROOT/'receipts/processes').glob(label+'_*.json')]
    states.append(dict(label=label,state_sha256=sha,accepted_stations=len(trace)-1,load_increment=m['parameters']['initial_increment'],
        radial_pitch_m=m['parameters']['radial_pitch'],interface_quadrature_order=m['parameters']['quadrature'],
        reference_mass_kg=m['reference_mass_kg'],prepared_phi_J=m['prepared_slit_phi_J'],
        maximum_recorded_process_rss_MiB=max(p['peak_rss_MiB'] for p in process_records),
        maximum_recorded_process_wall_s=max(p['wall_s_this_process'] for p in process_records),last_process_wall_s=m['wall_s_this_process'],
        final=trace[-1],source_hashes=m['parameters']['source_hashes']))
comparisons={name:json.loads((ROOT/'receipts'/f'compare_{name}.json').read_text()) for name in ['radial_refinement','load_increment','penalty_low','penalty_high','quadrature']}
assert all(all(r['engineering_comparison_gates'].values()) for r in comparisons.values())
preferred=json.loads((ROOT/'receipts/increment_half_inspection.json').read_text())
stability=json.loads((ROOT/'receipts/increment_half_stability.json').read_text())
assert preferred['active_zone_minimum_six_intervals_met'] and stability['status']=='PASS'
assert preferred['fully_separated_new_area_m2']==0.
snapshot=ROOT/'checkpoints/central_solved_20261006';sm=json.loads((snapshot/'manifest.json').read_text())
assert hashlib.sha256((snapshot/'manifest.json').read_bytes()).hexdigest()=='b81fe7bed3710cc51cef359cc96d202986a7fa7ca1e4246a74cf615f93895130'
assert all(hashlib.sha256((snapshot/p).read_bytes()).hexdigest()==h for p,h in sm['files'].items())
source_index=[];seen=set()
for state in states:
    for original,expected in state['source_hashes'].items():
        if (original,expected) in seen:continue
        seen.add((original,expected));candidates=[Path(original)]+list(ROOT.rglob(Path(original).name))
        found=next((p for p in candidates if p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest()==expected),None)
        assert found is not None,(original,expected)
        source_index.append(dict(original=original,sha256=expected,preserved_path=str(found.relative_to(ROOT))))
report=dict(status='NUMERICAL_COMPARISONS_PASS_INTERNAL_CONTROL_ONLY',states=states,comparisons=comparisons,
            exact_recorded_source_index=source_index,
            preferred_state='increment_half',preferred_inspection=preferred,preferred_local_stability=stability,
            first_central_snapshot_unchanged=True,gray_coupon_created=False,image_requested=False,appearance_promotion=False,
            physical_Gc_sensitivity_bracket_run=False,
            remaining_scope_limits=['Conditional processed-cork-informed parameters, not measured virgin-cork calibration',
                'One prescribed path and prepared slit; no prediction of bark crack pattern or complete growth',
                'Fixed inner substrate and 12 mm axial repetition limit admissible modes',
                'Radial pitch refined with fixed circumferential mesh and axial period; no claim of all-direction mesh convergence',
                'Positive numerical local stability estimate does not prove a global minimum',
                'Damage is stiffness loss; no newly fully separated ligament extension occurred',
                'No evidence of visual improvement over the retained baseline; no material scene promotion'])
(ROOT/'receipts/final_qualification.json').write_text(json.dumps(report,indent=2))
summary=[]
for name,c in comparisons.items():
    summary.append(f"- {name}: outer-edge gap change {100*c['maximum_gap_relative_change']:.6g}%; work change {100*c['total_work_relative_change']:.6g}%; peak conjugate change {100*c['maximum_work_conjugate_relative_change']:.6g}%.")
gap=comparisons['load_increment']['maximum_gap_m'][1]
text='# Frozen internal numerical qualification\n\nThe predeclared numerical comparisons pass for this conditional single-path control. No image, gray coupon or material promotion is produced.\n\n'+'\n'.join(summary)+'\n\n'
text+=f"The preferred 0.5 mm / 0.025-load-increment state has a {gap*1000:.6f} mm maximum outer-edge gap at the common material nodes. Its work residual is {preferred['relative_work_error']*100:.6f}%, and its active zone spans {preferred['active_zone_point_span_m']*1000:.6f} mm over {preferred['active_radial_intervals']} radial intervals. New fully separated extension is zero. The smallest estimated generalized local eigenvalue is {stability['smallest_generalized_eigenvalue']:.6g}.\n\n"
text+='Each quadrature rule and physical mesh starts from its own analytical prepared slit; no evolved history is transferred between material-point identities. The lower/higher penalty runs keep T0 and Gc fixed. The physical Gc sensitivity proposal is not run under the later narrowed scope.\n\nAll source/parameter assumptions and scope limits remain explicit in README.md and the final JSON receipt. The first protected central snapshot and historical failed attempts are preserved. The strongest visual baseline is untouched.\n'
(ROOT/'QUALIFICATION.md').write_text(text)
for p in (ROOT/'source').glob('*.py'):py_compile.compile(str(p),doraise=True)
files={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(ROOT.rglob('*')) if p.is_file() and '__pycache__' not in str(p) and p!=ROOT/'receipts/final_manifest.json'}
manifest=dict(status='FROZEN_INTERNAL_NUMERICAL_CONTROL',files=files,image_or_scene_promoted=False,first_central_snapshot_unchanged=True)
path=ROOT/'receipts/final_manifest.json';path.write_text(json.dumps(manifest,indent=2))
print(json.dumps(dict(manifest=str(path),manifest_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),files=len(files),status=report['status']),indent=2))
