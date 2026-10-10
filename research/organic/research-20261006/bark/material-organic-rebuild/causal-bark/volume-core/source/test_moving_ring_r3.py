from pathlib import Path
import json,hashlib,numpy as np
from ring_system_r3 import RingSystem
R=Path(__file__).resolve().parents[1];m=RingSystem(4,0.);rng=np.random.default_rng(1292)
E,g,_,_=m.evaluate();v=rng.normal(0,.007,m.q.shape);ac=np.zeros_like(v);ac[~m.free_mask]=rng.normal(0,14.,ac[~m.free_mask].shape);alpha=231.
a,res=m.moving_acceleration_reaction(g,v,ac,alpha)
free_res=float(np.linalg.norm(res[m.free_mask]));total_scale=float(np.linalg.norm(g)+np.linalg.norm(m.mass_product(a))+alpha*np.linalg.norm(m.mass_product(v)))
assert free_res/total_scale<1e-13
Edot=float(np.sum(v*(g+m.mass_product(a))));heat=alpha*float(np.sum(v*m.mass_product(v)));power=float(np.sum(v[~m.free_mask]*res[~m.free_mask]));err=abs(Edot+heat-power)/max(abs(Edot),abs(power),heat,1e-20);power_operand_scale=float(np.sum(np.abs(v*g))+np.sum(np.abs(v*m.mass_product(a)))+heat+abs(power));assert abs(Edot+heat-power)<100*np.finfo(float).eps*power_operand_scale
# Omit the formerly unnecessary moving-boundary term: reproduce the wrong
# equation on this same non-diagonal consistent mass and prescribed acceleration.
wrong=m.acceleration(g)-alpha*np.where(m.free_mask,v,0.)+ac
wrong_res=m.mass_product(wrong)+alpha*m.mass_product(v)+g
wrong_error=float(np.linalg.norm(wrong_res[m.free_mask]));assert wrong_error>free_res*1e6
# Zero prescribed motion reduces exactly to the fixed-boundary evolution.
z=np.zeros_like(v);vf=np.where(m.free_mask,v,0.);a0,r0=m.moving_acceleration_reaction(g,vf,z,alpha);assert np.allclose(a0,m.acceleration(g)-alpha*vf,rtol=1e-13,atol=1e-11)
# Full kinetic energy must include free/fixed coupling, not just M_ff.
K=m.kinetic(v);Kparts=m.kinetic(np.where(m.free_mask,v,0.))+m.kinetic(np.where(m.free_mask,0.,v));cross=K-Kparts;assert abs(cross)>1e-14
out={'all_pass':True,'free_force_balance_relative':free_res/total_scale,'full_power_identity_relative_to_net_power':err,'full_power_identity_relative_to_operand_scale':abs(Edot+heat-power)/power_operand_scale,'omitted_coupling_force_residual_N':wrong_error,'correct_free_force_residual_N':free_res,'kinetic_cross_term_J':cross,'source_sha256':hashlib.sha256(Path(__file__).with_name('ring_system_r3.py').read_bytes()).hexdigest(),'scope':'Algebraic force/power identities with true ring consistent mass. Time integrator and material calibration have separate gates.'}
(R/'receipts/moving_ring_r3_tests.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
