import json,time,resource
import numpy as np
from opening_model import OpeningModel,History,ROOT
from volume_coupon import Coupon

start=time.monotonic();m=OpeningModel();checks=[]
def check(name,value,limit):
    ok=bool(value<=limit);checks.append(dict(name=name,value=float(value),limit=float(limit),passed=ok))
    if not ok:raise AssertionError(checks[-1])
zero=np.zeros_like(m.X);h0=m.initial_history
E,g,h,r,o=m.evaluate(zero,h0,work=True)
check('initial_stored_energy_J',abs(r['stored_energy_J']),1e-22)
check('initial_new_dissipation_J',abs(r['new_dissipation_from_initial_J']),1e-18)
check('initial_incremental_potential_J',abs(E),1e-22)
check('prepared_slit_area_m2',abs(m.area[m.starter].sum()-.001*.012),1e-15)
check('prepared_slit_phi_J',abs(m.initial_phi-m.contract.Gc_J_m2*.001*.012),1e-14)
check('cut_topology_mass_kg',abs(m.seam.mass-m.control.mass),1e-15)
config=m.base.config.copy();config.pop('radial_pitch_m')
unsplit=Coupon(.006,precondition=False,config_override=config)
check('radial_subdivision_mass_relative',abs(unsplit.mass-m.seam.mass)/m.seam.mass,1e-12)
rng=np.random.default_rng(7443);v=rng.normal(size=m.X.shape)
check('cut_periodic_consistent_kinetic_relative',abs(m.seam.kinetic(v)-m.base.kinetic(m.seam.expand(v)))/m.seam.kinetic(v),1e-13)
clone=m.seam.bind(m.contract);phi_errors=[]
for k in (0.,.5*m.contract.delta_0,m.contract.delta_0,1.3*m.contract.delta_0,.2*m.contract.delta_f,m.contract.delta_f,1.1*m.contract.delta_f):
    ks=np.full(clone.N,k);clone.restore(ks,clone._law(ks))
    phi_errors.append(abs(clone.dissipation()-float(m.area@m.phi(ks)))/max(float(m.area.sum())*m.contract.Gc_J_m2,1e-20))
check('analytic_Phi_matches_native_dissipation_relative',max(phi_errors),1e-13)
# A small admissible gap field activates damage without leaving the bulk domain.
full=np.zeros_like(m.base.X);ids=m.seam.duplicated_nodes;rad=m.base.X[ids,2];inner=m.base.config['outer_radius_m']-sum(m.base.config['cohort_thickness_m'])
profile=(rad-inner)/sum(m.base.config['cohort_thickness_m'])
full[ids,0]=-.5*profile;full[len(m.control.base.X):,0]=.5*profile
counts=np.asarray(m.seam.T.power(2).sum(axis=0)).ravel();direction=(np.asarray(m.seam.T.T@full.ravel())/counts).reshape(m.X.shape)
direction[m.fixed_nodes]=0
m.set_load(.1);u=direction*1e-5
cases=[('active',h0),('unloading',History.make(np.where(m.starter,m.contract.delta_f,2e-5),m.bonds._law(np.where(m.starter,m.contract.delta_f,2e-5)))),
       ('saturated_unloading',History.make(np.full(m.bonds.N,m.contract.delta_f),np.ones(m.bonds.N)))]
results=[]
for name,accepted in cases:
    oldk=accepted.maximum.copy();oldd=accepted.damage.copy();E,g,h,r,o=m.evaluate(u,accepted)
    errors=[]
    for k in range(8):
        d=rng.normal(size=u.shape);d[m.fixed_nodes]=0;d/=np.linalg.norm(d);eps=1e-9
        numeric=(m.evaluate(u+eps*d,accepted)[0]-m.evaluate(u-eps*d,accepted)[0])/(2*eps);analytic=float(np.sum(g*d))
        errors.append(abs(numeric-analytic)/max(abs(numeric),abs(analytic),1e-6))
    check(name+'_full_incremental_gradient_relative',max(errors),1e-4)
    check(name+'_accepted_history_unchanged',max(np.max(abs(oldk-accepted.maximum)),np.max(abs(oldd-accepted.damage))),0)
    results.append(dict(case=name,gradient_error=max(errors),maximum_trial_damage=float(h.damage[~m.starter].max())))
# Explicit counterexample: stored energy alone has the wrong active derivative.
eps=1e-9;E,g,h,r,o=m.evaluate(u,h0)
rp=m.evaluate(u+eps*direction,h0)[3];rm=m.evaluate(u-eps*direction,h0)[3]
stored_numeric=(rp['stored_energy_J']-rm['stored_energy_J'])/(2*eps);full_analytic=float(np.sum(g*direction))
missing_phi_error=abs(stored_numeric-full_analytic)/max(abs(full_analytic),1e-12)
assert missing_phi_error>.1
# Scalar reference covers junctions, saturation and irreversible unloading.
scalar=[];T=m.contract.peak_Pa;K=m.contract.penalty_Pa_per_m;d0=m.contract.delta_0;df=m.contract.delta_f
for label,delta,kold in [('elastic',.5*d0,0),('onset',d0,0),('active',.3*df,0),('complete',df,0),('saturated',1.1*df,0),('unload',.1*df,.5*df)]:
    def value(x):
        kk=max(kold,x);dd=float(m.bonds._law(np.full(m.bonds.N,kk))[0])
        return .5*(1-dd)*K*x*x+float(m.phi(np.array([kk]))[0]-m.phi(np.array([kold]))[0])
    eps=1e-10;numeric=(value(delta+eps)-value(delta-eps))/(2*eps);dd=float(m.bonds._law(np.full(m.bonds.N,max(kold,delta)))[0]);analytic=(1-dd)*K*delta
    error=abs(numeric-analytic)/T;assert error<3e-5
    scalar.append(dict(case=label,error_over_peak=error))
result=dict(status='PASS',checks=checks,coupled_cases=results,stored_only_active_derivative_relative_error=missing_phi_error,scalar_branch_checks=scalar,
            cells=len(m.base.cells),full_nodes=len(m.base.X),periodic_nodes=len(m.X),radial_subdivisions=m.base.subdivisions.tolist(),
            prepared_slit_phi_J=m.initial_phi,prepared_slit_area_m2=float(m.area[m.starter].sum()),
            accepted_history_readonly=not h0.maximum.flags.writeable and not h0.damage.flags.writeable,
            wall_s=time.monotonic()-start,peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
(ROOT/'receipts/incremental_potential_audit.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
