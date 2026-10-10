import hashlib,json,resource,time
import numpy as np
from single_seam import SingleSeam,ConditionalFracture,ROOT
from periodic_sector import rotation_y

start=time.monotonic();seam=SingleSeam();checks=[]
def check(name,value,limit):
    ok=bool(value<=limit);checks.append(dict(name=name,value=float(value),limit=float(limit),passed=ok))
    if not ok:raise AssertionError(checks[-1])
b=seam.base;control=seam.control
check('split_topology_preserves_reference_nodes_m',np.max(abs(b.X[seam.cell_map]-b.Xcell)),1e-14)
check('split_topology_preserves_reference_mass_kg',abs(seam.mass-control.mass),1e-15)
check('single_seam_area_m2',abs(seam.geometry['area'].sum()-sum(b.config['cohort_thickness_m'])*b.config['height_m']),1e-14)
check('nonseam_nodes_unchanged',np.count_nonzero(seam.cell_map[~np.isclose(b.Xcell[...,0],0,atol=1e-12)]!=control.base.map[~np.isclose(b.Xcell[...,0],0,atol=1e-12)]),0)
rng=np.random.default_rng(67313);v=rng.normal(size=seam.X.shape)
check('split_rotation_condensed_kinetic_relative',abs(seam.kinetic(v)-b.kinetic(seam.expand(v)))/seam.kinetic(v),1e-13)
# Material values below are an operator-only unit fixture, never a cork scenario.
fixture=ConditionalFracture(1.,1.,1000.,'Operator-only unit fixture. These values must not enter a material solve.',evidence_status='synthetic unit test only')
bonds=seam.bind(fixture);u=np.zeros_like(seam.X);E,g,opening=seam.interface(u,bonds)
check('initial_interface_energy_J',abs(E),1e-30);check('initial_interface_force_N',np.max(abs(g)),1e-25)
local=b.Xcell[seam.geometry['pairs']].reshape(-1,36,3);gap=np.einsum('qi,qia->qa',seam.geometry['J'],local-local[:,:1])
check('initial_cohesive_reference_gap_m',np.max(np.linalg.norm(gap,axis=1)),1e-14)
v=rng.normal(size=seam.X.shape)*1e-6
E,g,opening=seam.interface(v,bonds);errors=[]
for k in range(8):
    d=rng.normal(size=v.shape);d/=np.linalg.norm(d);h=1e-8
    numerical=(seam.interface(v+h*d,bonds)[0]-seam.interface(v-h*d,bonds)[0])/(2*h);analytic=float(np.sum(g*d))
    errors.append(abs(numerical-analytic)/max(abs(numerical),abs(analytic),1e-14))
check('single_surface_condensed_interface_gradient_relative',max(errors),2e-5)
k0=bonds.maximum_opening.copy();d0=bonds.damage.copy()
proposed,damage=seam.propose_history(bonds,k0,np.full(bonds.N,fixture.delta_0*1.5))
seam.interface(v,bonds,damage)
check('rejected_trial_keeps_accepted_maximum',np.max(abs(bonds.maximum_opening-k0)),0)
check('rejected_trial_keeps_accepted_damage',np.max(abs(bonds.damage-d0)),0)
check('proposal_is_distinct_irreversible_history',not (np.all(proposed>k0) and np.all(damage>d0)),0)
check('fixture_integrated_bilinear_energy',abs(fixture.report()['bilinear_integrated_work_J_m2']-fixture.Gc_J_m2),1e-14)
# Reproduce the intact qualified state with zero gap in the split topology.
state=np.load(ROOT/'state/sector48/restart.npz');full=state['full_displacement'];full_split=np.vstack((full,full[seam.duplicated_nodes]))
counts=np.asarray(seam.T.power(2).sum(axis=0)).ravel()
u=(np.asarray(seam.T.T@full_split.ravel())/counts).reshape(seam.X.shape)
b.set_load(1.);Eb,gb,r=seam.bulk(u)
reference=json.loads((ROOT/'state/sector48/restart.json').read_text())['trace'][-1]
check('intact_state_energy_unchanged_by_split_relative',abs(Eb-reference['stored_energy_J'])/Eb,1e-13)
Ei,gi,opening=seam.interface(u,bonds)
check('intact_state_single_surface_opening_m',np.max(opening),1e-14)
result=dict(status='PASS',checks=checks,cohesive_faces=len(seam.faces)*3,cohesive_quadrature_points=bonds.N,
            duplicate_full_nodes=len(seam.duplicated_nodes),full_nodes=len(b.X),periodic_nodes=len(seam.X),
            interface_area_m2=float(seam.geometry['area'].sum()),reference_mass_kg=seam.mass,
            fixture=fixture.report(),physical_fracture_contract_bound=False,crack_evolution_solved=False,
            topology='Only x=0 longitudinal radial plane duplicated. Every other element face and all cohort interfaces remain welded. Circumferential and axial periodicity retained.',
            wall_s=time.monotonic()-start,peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
np.savez_compressed(ROOT/'state/single_seam_topology.npz',reference_nodes=b.X,node_map=seam.cell_map,periodic_groups=seam.groups,rotations=seam.rotations,
                    **{f'interface_{key}':value for key,value in seam.geometry.items()})
(ROOT/'receipts/single_seam_audit.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
