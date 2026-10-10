from pathlib import Path
import json,time,resource,hashlib
import numpy as np
from model import Coupon,DomainError,ROOT

start=time.monotonic(); m=Coupon(); checks=[]
def check(name,value,limit):
    passed=bool(value<=limit)
    checks.append(dict(name=name,value=float(value),limit=float(limit),passed=passed))
    if not passed: raise AssertionError(checks[-1])

u=np.zeros_like(m.X); E,g,r=m.evaluate(u,work=True)
check('initial_energy_J',abs(E),1e-20)
check('initial_gradient_N',np.max(np.abs(g)),1e-14)
check('initial_work_conjugate_J',abs(r['eigenstrain_conjugate_J']),1e-20)
check('state_topology_mapping_m',m.mapping_error,1e-14)
check('condensed_mass_relative_error',abs(m.mass-sum(c.mass for c in m.cells))/m.mass,1e-13)
v=np.broadcast_to([.3,-.2,.1],m.X.shape).copy()
check('condensed_kinetic_constant_velocity',abs(m.kinetic(v)-.5*m.mass*.14),1e-14)
rng=np.random.default_rng(8831)
v=rng.normal(size=m.X.shape)
cellv=v[m.map]; Mc=np.array([c.mass_scalar for c in m.cells])
check('condensed_kinetic_random_velocity_relative',abs(m.kinetic(v)-.5*np.einsum('cia,cij,cja->',cellv,Mc,cellv))/m.kinetic(v),1e-13)
check('proper_reference_material_frames',np.max(np.abs(np.linalg.det(m.axes)-1)),1e-12)
m.set_load(.4)
u=rng.normal(size=m.X.shape)*1e-6; u[m.fixed_nodes]=0
E,g,r=m.evaluate(u,work=True)
gradient_errors=[]
for k in range(12):
    d=rng.normal(size=u.shape);d[m.fixed_nodes]=0;d/=np.linalg.norm(d)
    h=1e-8
    numeric=(m.evaluate(u+h*d)[0]-m.evaluate(u-h*d)[0])/(2*h)
    analytic=float(np.sum(g*d))
    gradient_errors.append(abs(numeric-analytic)/max(abs(numeric),abs(analytic),1e-8))
check('global_directional_gradient_max_relative',max(gradient_errors),2e-5)
s=m.load; h=1e-5
m.set_load(s+h); plus=m.evaluate(u)[0]
m.set_load(s-h); minus=m.evaluate(u)[0]
m.set_load(s)
check('imposed_eigenstrain_work_derivative_relative',abs((plus-minus)/(2*h)-r['eigenstrain_conjugate_J'])/abs(r['eigenstrain_conjugate_J']),1e-8)
check('reference_volume_fixed_across_load',max(abs(op.reference_volume.sum()-m.weights[l*m.T:(l+1)*m.T].sum()) for l,op in enumerate(m.operators)),1e-18)
check('mass_fixed_across_load',abs(sum(op.mass.sum() for op in m.operators)-m.mass),1e-15)
# These operators see the elastic determinant; physical density must not use it.
zero=np.zeros_like(m.X); _,_,rz=m.evaluate(zero,work=True)
check('physical_volume_unmoved_eigenstrain_relative',abs(rz['physical_volume_m3']-m.volume)/m.volume,1e-13)
check('physical_density_unmoved_eigenstrain_relative',abs(rz['mean_current_density_kg_m3']-m.config['reference_density_kg_m3'])/m.config['reference_density_kg_m3'],1e-13)
# Guard rejects, rather than silently evaluating beyond, the measured endpoint.
m.set_load(1.2)
guard=False
try:m.evaluate(zero)
except DomainError:guard=True
check('measured_endpoint_guard',not guard,0)
m.set_load(.05)
pilot_started=time.monotonic()
for i in range(10):m.evaluate(zero)
evaluation_s=(time.monotonic()-pilot_started)/10
pilot_started=time.monotonic(); solution,pilot=m.equilibrate(zero); pilot_s=time.monotonic()-pilot_started
check('pilot_relative_residual',pilot['relative_energy_residual'],m.config['equilibrium_relative_energy_norm_tolerance'])
check('pilot_maximum_free_force',pilot['maximum_free_force_N'],m.config['equilibrium_maximum_force_N'])
result=dict(status='PASS',checks=checks,build_s=m.build_s,nodes=len(m.X),cells=len(m.cells),quadrature_points=int(m.weights.size),
            reference_mass_kg=m.mass,reference_volume_m3=m.volume,evaluation_s=evaluation_s,pilot_s=pilot_s,pilot=pilot,
            wall_s=time.monotonic()-start,peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
            estimated_32_station_force_wall_s=32*pilot['line_search_evaluations']*evaluation_s,
            source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(ROOT/'receipts/audit_pilot.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
