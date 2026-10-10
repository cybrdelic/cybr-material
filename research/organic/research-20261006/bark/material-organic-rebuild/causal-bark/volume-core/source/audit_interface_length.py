"""Minimal exact bar analogue for the ring's intrinsic cohesive compliance.
No plot, tuning, or fracture morphology. Fixed K versus mesh-scaled K is tested
against the series spring solution under a specified traction.
"""
from pathlib import Path
import json
R=Path(__file__).resolve().parents[1]
E=2.8e6;L=.064;eta=.01;unrelated_height=.008;sigma=2e4
Kfixed=E/(eta*unrelated_height)
rows=[]
for n in (4,8,16,32,64):
 h=L/n;bulk=L/E;fixed=(n-1)/Kfixed;scaled=(n-1)/(E/(eta*h));rows.append({'cells':n,'bulk_compliance_m_per_Pa':bulk,'fixed_K_relative_added_compliance':fixed/bulk,'normal_cell_length_K_relative_added_compliance':scaled/bulk,'fixed_K_extension_m':sigma*(bulk+fixed),'scaled_K_extension_m':sigma*(bulk+scaled)})
out={'finding':'ring_interfaces.py lines 14-15 use axial span to set every fracture-face stiffness, including circumferential-normal faces. Thus spatial refinement changes artificial interfacial compliance, unless it is separately taken to zero. This is a qualification gap rather than proof of wrong finite-K constitutive behavior.','minimal_repro':'N equal bars in series with N-1 intrinsic cohesive springs: C=L/E+(N-1)/K.','parameters_SI':{'E_Pa':E,'L_m':L,'unrelated_axial_height_m':unrelated_height,'compliance_fraction':eta,'traction_Pa':sigma},'results':rows,'conclusion':'A mesh-objective penalty limit needs normal-length scaling and an independent compliance refinement. Alternatively a measured physical interface stiffness must be assigned to actual anatomical interfaces; arbitrary cell boundaries cannot create additional physical layers.','source_fix_status':'Not yet changed. Keep existing ring experiments labeled finite-K diagnostic fixtures.'}
(R/'receipts/intrinsic_interface_compliance_audit.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
