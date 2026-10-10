import json,resource,time
import numpy as np
from periodic_sector import PeriodicSector,rotation_y,ROOT

start=time.monotonic();m=PeriodicSector();checks=[]
def check(name,value,limit):
    ok=bool(value<=limit);checks.append(dict(name=name,value=float(value),limit=float(limit),passed=ok))
    if not ok:raise AssertionError(checks[-1])
check('periodic_reference_reconstruction_m',m.reference_reconstruction_error,2e-12)
check('reference_mass_kg',abs(m.mass-m.base.weights.sum()*m.config['reference_density_kg_m3']),1e-15)
rng=np.random.default_rng(12813);v=rng.normal(size=m.X.shape)
check('condensed_consistent_kinetic_relative',abs(m.kinetic(v)-m.base.kinetic(m.expand(v)))/m.kinetic(v),1e-13)
v[:]=[0.,.7,0.]
check('axial_translation_kinetic_kg_m2_s2',abs(m.kinetic(v)-.5*m.mass*.7**2),1e-15)
R=rotation_y(.37);v=m.X@(R-np.eye(3)).T
E,g,r=m.evaluate(v)
check('finite_stem_rigid_rotation_energy_J',abs(E),1e-22)
check('finite_stem_rigid_rotation_free_force_N',r['maximum_free_force_N'],1e-9)
m.set_load(.3);zero=np.zeros_like(m.X)
E0,g0,r0=m.evaluate(zero,work=True);E1,g1,r1=m.evaluate(v,work=True)
check('prestrained_rigid_rotation_energy_relative',abs(E1-E0)/E0,1e-12)
check('prestrained_rigid_rotation_force_relative',np.linalg.norm(g1-g0@R.T)/np.linalg.norm(g0),1e-11)
check('prestrained_rigid_rotation_work_relative',abs(r1['eigenstrain_conjugate_J']-r0['eigenstrain_conjugate_J'])/abs(r0['eigenstrain_conjugate_J']),1e-12)
v=rng.normal(size=m.X.shape)*1e-6;v[m.fixed_nodes]=0
E,g,r=m.evaluate(v,work=True);errors=[]
for k in range(8):
    d=rng.normal(size=v.shape);d[m.fixed_nodes]=0;d/=np.linalg.norm(d);h=1e-8
    numeric=(m.evaluate(v+h*d)[0]-m.evaluate(v-h*d)[0])/(2*h);analytic=float(np.sum(g*d))
    errors.append(abs(numeric-analytic)/max(abs(numeric),abs(analytic),1e-8))
check('rotation_condensed_force_derivative_relative',max(errors),2e-5)
h=1e-5;s=m.load;m.set_load(s+h);plus=m.evaluate(v)[0];m.set_load(s-h);minus=m.evaluate(v)[0];m.set_load(s)
check('eigenstrain_work_derivative_relative',abs((plus-minus)/(2*h)-r['eigenstrain_conjugate_J'])/abs(r['eigenstrain_conjugate_J']),1e-8)
q=m.base.X+m.expand(v);canonical=m.X+v
expected=np.einsum('nij,nj->ni',m.rotations,canonical[m.groups]);expected[m.upper_nodes,1]+=m.config['height_m']
check('periodic_deformed_seam_continuity_m',np.max(np.linalg.norm(q-expected,axis=1)),2e-12)
m.set_load(.05);begin=time.monotonic();solution,pilot=m.equilibrate(zero);pilot_s=time.monotonic()-begin
check('pilot_relative_residual',pilot['relative_energy_residual'],m.config['equilibrium_relative_energy_norm_tolerance'])
check('pilot_force_N',pilot['maximum_free_force_N'],m.config['equilibrium_maximum_force_N'])
result=dict(status='PASS',checks=checks,sector_angle_rad=m.alpha,full_nodes=len(m.base.X),periodic_nodes=len(m.X),
            cells=len(m.base.cells),reference_mass_kg=m.mass,pilot_s=pilot_s,pilot=pilot,
            axial_assumption='Repeating 48 mm axial height with zero imposed mean axial extension; not a finite free-ended coupon',
            rigid_rotation_scope='Finite rotation about the stem axis, compatible with the rotational periodic chart. General material-point objectivity is inherited from the 57-check native port; no arbitrary observer rotation is newly claimed here.',
            wall_s=time.monotonic()-start,peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
(ROOT/'receipts/periodic_audit.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
