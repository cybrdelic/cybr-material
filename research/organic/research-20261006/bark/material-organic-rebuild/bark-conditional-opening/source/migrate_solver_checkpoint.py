"""Preserve identical physical state across an optimizer/report-only revision."""
import hashlib,importlib.util,json,shutil,sys
from pathlib import Path
import numpy as np
from opening_model import OpeningModel,History,ROOT
from checkpoint_io import read_checkpoint,write_checkpoint

source=ROOT/'state/central_r1/restart.npz';arrays,meta=read_checkpoint(source)
old_source=ROOT/'attempts/steepest_descent/source/opening_model.py'
assert hashlib.sha256(old_source.read_bytes()).hexdigest()==meta['parameters']['source_hashes'][str(ROOT/'source/opening_model.py')]
spec=importlib.util.spec_from_file_location('opening_prior_optimizer',old_source);module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
old=module.OpeningModel();new=OpeningModel();load=meta['trace'][-1]['load'];old.set_load(load);new.set_load(load)
assert np.array_equal(old.seam.cell_map,new.seam.cell_map) and np.array_equal(old.interface_reference,new.interface_reference)
assert old.seam.mass==new.seam.mass==meta['reference_mass_kg']
history=History.make(arrays['maximum'],arrays['damage']);a=old.evaluate(arrays['displacement'],history,work=True);b=new.evaluate(arrays['displacement'],history,work=True)
assert a[0]==b[0] and np.array_equal(a[1],b[1]) and np.array_equal(a[2].maximum,b[2].maximum) and np.array_equal(a[2].damage,b[2].damage)
assert a[3]['eigenstrain_conjugate_J']==b[3]['eigenstrain_conjugate_J']
dependencies=[ROOT/'source'/name for name in ('opening_model.py','volume_coupon.py','periodic_sector.py','single_seam.py','checkpoint_io.py','run_opening.py')]+list((ROOT/'config').glob('*.json'))
meta['parameters']['label']='central_r2';meta['parameters']['source_hashes']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in dependencies}
for row in meta['trace']:
    row['fracture_dissipation_J']=row['new_dissipation_from_initial_J'];row.pop('damage',None)
meta['status']='paused';meta['error']=None
meta['events'].append(dict(event='optimizer_and_diagnostic_revision',load=load,physical_state_changed=False,history_quadrature_identity_changed=False,
    explanation='Replace slow preconditioned steepest descent with preconditioned L-BFGS on the same incremental potential. Correct inherited bulk-only fracture label. Energy, full gradient, history and work conjugate compare bit-for-bit.'))
dest=ROOT/'state/central_r2';dest.mkdir(exist_ok=True)
meta['state_sha256']=write_checkpoint(dest/'restart.npz',arrays,meta);(dest/'restart.json').write_text(json.dumps(meta,indent=2))
shutil.copyfile(ROOT/'state/central_r1/previous_accepted.npz',dest/'previous_accepted.npz')
result=dict(status='PASS',source_checkpoint_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),destination_checkpoint_sha256=meta['state_sha256'],
    load=load,energy_gradient_history_work_bitwise_equal=True,reference_mass_kg=new.seam.mass,mass_bitwise_equal=True,accepted_history_transferred_without_interpolation=True,
    identical_mesh_and_quadrature=True,no_damage_reset=True,no_physical_energy_changed=True,prior_source_preserved=str(old_source))
(ROOT/'receipts/optimizer_revision.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
