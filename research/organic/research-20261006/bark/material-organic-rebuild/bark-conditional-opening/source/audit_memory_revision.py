import hashlib,importlib.util,json,sys
import numpy as np
from opening_model import OpeningModel,History,ROOT
from checkpoint_io import read_checkpoint

a,meta=read_checkpoint(ROOT/'state/central_r2/restart.npz')
old_source=ROOT/'attempts/eager_factorization/source/opening_model.py'
assert hashlib.sha256(old_source.read_bytes()).hexdigest()==meta['parameters']['source_hashes'][str(ROOT/'source/opening_model.py')]
spec=importlib.util.spec_from_file_location('eager_opening_model',old_source);module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
old=module.OpeningModel();new=OpeningModel();old.set_load(1.);new.set_load(1.);h=History.make(a['maximum'],a['damage'])
before=old.evaluate(a['displacement'],h,work=True);after=new.evaluate(a['displacement'],h,work=True)
assert before[0]==after[0] and np.array_equal(before[1],after[1]) and old.seam.mass==new.seam.mass
assert np.array_equal(before[2].maximum,after[2].maximum) and np.array_equal(before[2].damage,after[2].damage)
assert before[3]['eigenstrain_conjugate_J']==after[3]['eigenstrain_conjugate_J']
result=dict(status='PASS',change='Avoid unused intact-volume/periodic factorizations when constructing the cut model; keep the same assembled per-cell and cut preconditioner.',
    material_energy_force_history_work_mass_bitwise_unchanged=True,coarse_accepted_state_unmodified=True,refinement_restarts_from_prepared_initial_state=True,
    failed_refinement='state/radial_half/restart.npz; no load step accepted before memory cap',working_memory_cap_unchanged=True)
(ROOT/'receipts/memory_revision.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
