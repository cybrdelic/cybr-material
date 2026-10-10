"""Integrate native traction, without a bulk solve or material parameter choice."""
import json,time
import numpy as np
from single_seam import SingleSeam,ConditionalFracture,ROOT

start=time.monotonic();seam=SingleSeam();fixture=ConditionalFracture(1.,1.,1000.,'Only an independent numerical operator-work integration; never a cork material input.',evidence_status='synthetic unit test only')
bonds=seam.bind(fixture);area=float(seam.geometry['area'].sum())
full=np.zeros_like(seam.base.X);full[seam.duplicated_nodes,0]=-.5;full[len(seam.control.base.X):,0]=.5
counts=np.asarray(seam.T.power(2).sum(axis=0)).ravel()
direction=(np.asarray(seam.T.T@full.ravel())/counts).reshape(seam.X.shape)
opening=np.r_[np.linspace(0,fixture.delta_0,51),np.linspace(fixture.delta_0,fixture.delta_f,951)[1:]]
traction=[]
for delta in opening:
    damage=bonds._law(np.full(bonds.N,delta))
    E,g,actual=seam.interface(delta*direction,bonds,damage)
    assert np.max(abs(actual-delta))<1e-12
    traction.append(float(np.sum(g*direction))/area)
integrated=float(np.sum(.5*(np.asarray(traction)[1:]+np.asarray(traction)[:-1])*np.diff(opening)));relative=abs(integrated-fixture.Gc_J_m2)/fixture.Gc_J_m2
assert relative<1e-10
assert np.array_equal(bonds.maximum_opening,np.zeros(bonds.N)) and np.array_equal(bonds.damage,np.zeros(bonds.N))
clone=seam.bind(fixture);clone.restore(np.full(clone.N,fixture.delta_f),np.ones(clone.N))
dissipation=clone.dissipation()/area
assert abs(dissipation-fixture.Gc_J_m2)<1e-12
result=dict(status='PASS',sample_count=len(opening),native_traction_integral_J_m2=integrated,declared_unit_fixture_Gc_J_m2=fixture.Gc_J_m2,
            relative_error=relative,independent_native_history_dissipation_J_m2=dissipation,accepted_history_untouched=True,
            no_bulk_evaluation_or_crack_evolution=True,physical_contract_bound=False,wall_s=time.monotonic()-start)
(ROOT/'receipts/cohesive_work_integral.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
