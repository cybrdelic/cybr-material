import json
import numpy as np
from periodic_sector import ROOT

a=json.loads((ROOT/'receipts/sector48.json').read_text());b=json.loads((ROOT/'receipts/sector96.json').read_text())
assert a['status']==b['status']=='completed'
fa=a['trace'][-1];fb=b['trace'][-1]
result=dict(status='PASS',common_load=1.,sector_angles_rad=[a['sector_angle_rad'],b['sector_angle_rad']],
    normalized_energy_relative_difference=abs(fa['stored_energy_J']-fb['stored_energy_J']/2)/fa['stored_energy_J'],
    normalized_work_relative_difference=abs(fa['eigenstrain_work_J']-fb['eigenstrain_work_J']/2)/fa['eigenstrain_work_J'],
    normalized_mass_relative_difference=abs(a['reference_mass_kg']-b['reference_mass_kg']/2)/a['reference_mass_kg'],
    maximum_tangential_strain_by_cohort=fa['maximum_tangential_strain_by_cohort'],
    physical_orientation_minimum=fa['minimum_physical_J'],
    eigenstrain_work_J=fa['eigenstrain_work_J'],stored_energy_J=fa['stored_energy_J'],
    relative_work_ledger_error=abs(fa['energy_work_error_J'])/fa['stored_energy_J'],
    boundary_interpretation='The periodic sector eliminates the fixed/free circumferential edge. Axial repetition eliminates the top/bottom cut. Rigid inner attachment remains an explicit idealization.',
    angular_scope='A local representative sector with angle width/radius. It is not a complete closed 2 pi trunk mesh. A later whole-trunk construction must select an integer tiling count and rerun this control at that commensurate angle.',
    stress_law_scope='Full original load is within the reported tangential endpoints and assumed 6% radial/axial scope. Native law and guards unchanged.')
for key in ('normalized_energy_relative_difference','normalized_work_relative_difference','normalized_mass_relative_difference'):assert result[key]<1e-10
(ROOT/'receipts/expanded_sector_comparison.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
