"""Bounded native-versus-frozen audit; material points and prism contraction only."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
import ctypes, hashlib, json, resource, sys, time
import numpy as np
from native_directional import NativeEnvelope, NativePrisms, DOUBLE, INT64, _ptr
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT.parent
sys.path.insert(0,str(BASE/'bark-constitutive-incremental/source'))
from incremental_envelope import IncrementalEnvelope
from directional_candidate import synthetic_candidate
from quadratic_prism_frozen import QuadraticPrism,flat_reference

FROZEN=[BASE/'bark-constitutive-upgrade/source/directional_candidate.py',BASE/'bark-constitutive-upgrade/source/quadratic_prism_frozen.py',BASE/'bark-constitutive-incremental/source/incremental_envelope.py',BASE/'causal-bark/volume-core/source/batch_volume_incremental.cpp',BASE/'causal-bark/volume-core/source/native_volume.py']
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def frozen(axes=None):
    l=synthetic_candidate(axes)
    return IncrementalEnvelope(l.tension,l.compression,l.G,l.axes,l.minimum_stretch,l.maximum_stretch)
def rotation(rng):
    q,_=np.linalg.qr(rng.normal(size=(3,3)));q[:,0]*=np.linalg.det(q);return q

def main():
    started=time.perf_counter();before={str(p.relative_to(BASE)):sha(p) for p in FROZEN}
    rng=np.random.default_rng(612879);native=NativeEnvelope();law=frozen();checks={};metrics={}
    def check(name,value,limit=None):
        ok=bool(value) if limit is None else bool(value<=limit)
        checks[name]=ok
        if limit is not None:metrics[name]={'observed':float(value),'upper_limit':float(limit)}
        if not ok:raise AssertionError((name,value,limit))
    def reject(name,fn):
        try:fn()
        except (ValueError,FloatingPointError):check(name,True)
        else:check(name,False)
    w0,p0,_=native.evaluate_incremental(np.zeros((3,3)))
    check('identity_exact_zero',float(w0)==0 and np.array_equal(p0,np.zeros((3,3))))
    h=rng.normal(size=(512,3,3))*.025
    ws,ps,ds=native.evaluate_incremental(h);wr,pr,dr=law.evaluate_incremental(h)
    check('batch_energy_numpy_relative',np.linalg.norm(ws-wr)/np.linalg.norm(wr),2e-13)
    check('batch_stress_numpy_relative',np.linalg.norm(ps-pr)/np.linalg.norm(pr),2e-13)
    check('batch_diagnostics_numpy_scaled',max(np.max(abs(ds[k]-dr[k]))/(1e6 if 'energy' in k else 1) for k in ds),2e-13)
    frames=[];finite_cases=0
    for _ in range(8):
        a=rotation(rng);q=rotation(rng);reference=frozen(a)
        fields=np.repeat(np.eye(3)[None],36,axis=0)
        fields[:,0,0]=np.linspace(.2,1.08,36);fields[:,1,1]=np.linspace(1.09,.75,36)
        fields[:,2,2]=np.linspace(.85,1.1,36);fields[:,0,1]=.035;fields[:,2,1]=-.025
        f=q@fields@a.T;h=f-np.eye(3)
        w,p,d=native.evaluate_incremental(h,a);we,pe,de=reference.evaluate_incremental(h)
        frames.append(max(np.max(abs(w-we))/1e6,np.max(abs(p-pe))/1e6));finite_cases+=len(h)
        wbase,pbase,_=native.evaluate(fields)
        frames.append(max(np.max(abs(w-wbase))/1e6,np.max(abs(p-q@pbase@a.T))/1e6))
    check('finite_strain_rotation_and_axis_parity',max(frames),3e-12)
    conditioned=[]
    for epsilon in (1e-2,1e-3,1e-4,1e-5,1e-6,1e-7):
        f=np.array([[.5,.35,.2],[0,epsilon,.1],[0,0,.7]])
        w,p,d=native.evaluate(f);we,pe,de=law.evaluate(f)
        conditioned.append(dict(transverse_component=epsilon,J=float(d['J']),energy_relative_error=abs(float(w-we))/float(we),stress_relative_error=float(np.linalg.norm(p-pe)/np.linalg.norm(pe))))
    check('nearly_coplanar_accepted_stress_numpy_parity',max(r['stress_relative_error'] for r in conditioned),2e-8)
    check('nearly_coplanar_accepted_energy_numpy_parity',max(r['energy_relative_error'] for r in conditioned),2e-13)
    metrics['nearly_coplanar_points']=conditioned
    bounds_match=True
    for boundary in (.15,1.15):
        for offset in (-1e-12,0.,1e-12):
            f=np.diag([boundary+offset,1.,1.]);accepted=[]
            for evaluator in (native,law):
                try:evaluator.evaluate(f);accepted.append(True)
                except ValueError:accepted.append(False)
            bounds_match &= accepted[0]==accepted[1]
    check('exact_and_neighboring_domain_boundary_parity',bounds_match)

    metrics['finite_rotated_points']=finite_cases
    rotations=np.array([rotation(rng) for _ in range(40)])
    w,p,_=native.evaluate(rotations)
    check('proper_rotation_zero_energy_stress_scaled',max(np.max(abs(w))/1e6,np.max(abs(p))/1e6),2e-14)
    symmetry=0;fd_error=0
    for _ in range(18):
        h=rng.normal(size=(3,3))*.024;v=rng.normal(size=(3,3));v/=np.linalg.norm(v)
        w,p,d=native.evaluate_incremental(h);step=1e-6
        fd=(native.evaluate_incremental(h+step*v)[0]-native.evaluate_incremental(h-step*v)[0])/(2*step)
        fd_error=max(fd_error,abs(float(fd)-np.sum(p*v))/1e6)
        sigma=p@(np.eye(3)+h).T/float(d['J']);symmetry=max(symmetry,np.max(abs(sigma-sigma.T))/1e6)
    check('material_energy_directional_gradient',fd_error,2e-8)
    check('Cauchy_stress_symmetry',symmetry,2e-13)
    small=[];covariance=[]
    for exponent in range(3,13):
        eps=10.**-exponent
        for sign in (-1.,1.):
            h=np.diag([sign*eps,0.,0.]);w,p,d=native.evaluate_incremental(h)
            expected=.5*1e6*eps*eps
            small.append(dict(kind='axial',eps=sign*eps,energy=float(w),relative_error=abs(float(w)-expected)/expected))
            for i in range(3):
                for j in range(3):
                    if i==j:continue
                    h=np.zeros((3,3));h[i,j]=sign*eps;w,p,d=native.evaluate_incremental(h);we,pe,de=law.evaluate_incremental(h)
                    expected=.5*law.G*np.log1p(eps*eps)
                    small.append(dict(kind=f'shear_{i}{j}',eps=sign*eps,energy=float(w),angular_relative_error=abs(float(d['angular_energy_J_m3'])-expected)/expected,stress_relative_error=float(np.linalg.norm(p-pe)/np.linalg.norm(pe)),energy_relative_error=abs(float(w-we))/abs(float(we))))
        a=rotation(rng);h=np.diag([eps,-.5*eps,.2*eps]);h[0,1]=.3*eps
        w,p,_=native.evaluate_incremental(h);wa,pa,_=native.evaluate_incremental(a@h@a.T,a)
        covariance.append(max(abs(float(w-wa))/float(w),np.linalg.norm(pa-a@p@a.T)/np.linalg.norm(p)))
    check('axial_energy_through_1e_minus_12',max(r['relative_error'] for r in small if r['kind']=='axial'),1e-12)
    check('all_signed_shear_energies_through_1e_minus_12',max(r['angular_relative_error'] for r in small if r['kind']!='axial'),1e-12)
    check('tiny_shear_numpy_stress_parity',max(r['stress_relative_error'] for r in small if r['kind']!='axial'),1e-12)
    check('tiny_shear_numpy_energy_parity',max(r['energy_relative_error'] for r in small if r['kind']!='axial'),1e-12)
    check('tiny_strain_material_spatial_covariance',max(covariance),1e-12)
    knots=[]
    for axis,knot in [(0,-.08),(0,-.55),(1,.03)]:
        for delta in (-1e-8,-1e-12,0,1e-12,1e-8):
            h=np.zeros((3,3));h[axis,axis]=knot+delta
            w,p,_=native.evaluate_incremental(h);we,pe,_=law.evaluate_incremental(h)
            knots.append(max(abs(float(w-we))/1e6,np.max(abs(p-pe))/1e6))
    # Include the series/log1p densification arithmetic junction u=0.001.
    for delta in (-1e-12,0,1e-12):
        h=np.diag([-(.55+.45*.001+delta),0,0]);w,p,_=native.evaluate_incremental(h);we,pe,_=law.evaluate_incremental(h)
        knots.append(max(abs(float(w-we))/1e6,np.max(abs(p-pe))/1e6))
    check('branch_and_series_junction_numpy_parity',max(knots),1e-13)
    for name,f in [('reflection',np.diag([-1,1,1])),('singular',np.diag([0,1,1])),('below_domain',np.diag([.149,1,1])),('above_domain',np.diag([1,1.151,1])),('nonfinite',np.full((3,3),np.nan)),('overflow',np.full((3,3),1e308)),('near_coplanar',np.array([[.5,.5,0],[0,1e-12,0],[0,0,1.]]))]:
        reject('reject_'+name,lambda f=f:native.evaluate(f))
    for name,a in [('left_handed',np.diag([-1,1,1])),('nonorthonormal',np.diag([1.001,1,1])),('nonfinite',np.full((3,3),np.inf))]:
        reject('reject_axes_'+name,lambda a=a:native.evaluate_incremental(np.zeros((3,3)),a))
    for name,h in [('empty',np.empty((0,3,3))),('shape',np.zeros((3,2))),('scalar',1.)]:reject('reject_input_'+name,lambda h=h:native.evaluate_incremental(h))
    params=native.parameters.copy();params[30]=-1.;bad=np.array([-1],np.int64);out=np.empty(1);p=np.empty((1,3,3));d=np.empty((1,6));h=np.zeros((1,3,3));a=np.eye(3)[None].copy()
    rc=native.point_function(1,_ptr(h),_ptr(a),_ptr(params),_ptr(out),_ptr(p),_ptr(d),bad.ctypes.data_as(INT64))
    check('raw_native_invalid_parameter_rejection',rc==-3)

    X=flat_reference(np.array([[0.,0.],[.001,0.],[0.,.001]]),.001)
    cells=[QuadraticPrism(X,1e6,0,240.,nq=4),QuadraticPrism(1.2*X,1e6,0,310.,nq=4)]
    gradients=np.array([c.grad for c in cells]);weights=np.array([c.weights for c in cells]);rho=np.array([c.rho for c in cells])
    axes=np.array([np.eye(3),rotation(rng)])[:,None]
    f0=np.array([[.74,.03,.01],[0,1.024,.01],[.01,0,.97]])
    h0=np.broadcast_to(f0-np.eye(3),(2,len(cells[0].weights),3,3)).copy()
    # The second material frame is independently rotated; all axes are explicit.
    prism=NativePrisms(gradients,weights,h0,axes,rho)
    du=rng.normal(size=(2,18,3))*1e-6
    energy,gradient,diag=prism.evaluate(du)
    def expected(u,origin=h0):
        h=origin+np.einsum('cia,cqib->cqab',u-u.mean(axis=1,keepdims=True),gradients)
        ws=[];ps=[];js=[]
        for c in range(2):
            wc,pc,dc=frozen(axes[c,0]).evaluate_incremental(h[c]);ws.append(wc);ps.append(pc);js.append(dc['J'])
        eg=np.einsum('cqib,cqab,cq->cia',gradients,np.array(ps),weights);eg-=eg.mean(axis=1,keepdims=True)
        return float(np.sum(weights*np.array(ws))),eg,np.array(js)
    en,gr,j=expected(du)
    check('nonaffine_prism_numpy_energy_relative',abs(energy-en)/abs(en),3e-13)
    check('nonaffine_prism_numpy_gradient_relative',np.linalg.norm(gradient-gr)/np.linalg.norm(gr),3e-13)
    check('nonaffine_prism_numpy_J',np.max(abs(diag['J']-j)),3e-14)
    fd=np.empty_like(du);step=2e-10
    for c in range(2):
        for i in range(18):
            for a in range(3):
                perturb=np.zeros_like(du);perturb[c,i,a]=step
                fd[c,i,a]=(prism.evaluate(du+perturb)[0]-prism.evaluate(du-perturb)[0])/(2*step)
    check('all_108_nodal_energy_gradients_scaled',np.max(abs(fd-gradient)),2e-7)
    check('balanced_nodal_force_scaled',np.max(abs(gradient.sum(axis=1))),2e-13)
    # Test the exact derivative of centering even if supplied gradients have a
    # small nonzero mean; no hidden assumption about partition-of-unity cleanup.
    disturbed=gradients+.125
    centered=NativePrisms(disturbed,weights,h0,axes,rho);v=rng.normal(size=du.shape);v/=np.linalg.norm(v)
    ec,gc,_=centered.evaluate(du);ep=centered.evaluate(du+step*v)[0];em=centered.evaluate(du-step*v)[0]
    check('centering_derivative_for_nonzero_gradient_sum',abs((ep-em)/(2*step)-np.sum(gc*v)),2e-7)
    dz=np.zeros_like(du);e0,g0,d0=prism.evaluate(dz);ew,gw,_=expected(dz)
    check('origin_prestrain_preserved',abs(e0-ew)/abs(ew),2e-13)
    check('origin_stress_preserved',np.linalg.norm(g0-gw)/np.linalg.norm(gw),2e-13)
    reset=NativePrisms(gradients,weights,np.zeros_like(h0),axes,rho)
    er,_,_=reset.evaluate(dz)
    check('prestrain_is_not_reset_to_stress_free',e0>0 and er==0.)
    du1=du*.37;du2=du-du1
    next_h0=h0+np.einsum('cia,cqib->cqab',du1-du1.mean(axis=1,keepdims=True),gradients)
    transferred=NativePrisms(gradients,weights,next_h0,axes,rho)
    et,gt,dt=transferred.evaluate(du2)
    check('origin_transfer_energy_relative',abs(et-energy)/energy,3e-13)
    check('origin_transfer_gradient_relative',np.linalg.norm(gt-gradient)/np.linalg.norm(gradient),3e-13)
    check('origin_transfer_does_not_reset_mass',np.array_equal(dt['cell_material_mass_kg'],diag['cell_material_mass_kg']))
    trans=du+np.array([.013,-.027,.019])
    et,gt,_=prism.evaluate(trans)
    check('local_increment_translation_energy_relative',abs(et-energy)/energy,1e-11)
    check('local_increment_translation_gradient_relative',np.linalg.norm(gt-gradient)/np.linalg.norm(gradient),1e-11)
    check('fixed_reference_mass',np.max(abs(diag['cell_material_mass_kg']-rho*weights.sum(axis=1))),0.)
    integrated=np.sum(diag['quadrature_current_density_kg_m3']*diag['J']*weights,axis=1)
    check('current_density_integrates_fixed_mass_relative',np.max(abs(integrated/prism.mass-1)),3e-14)
    check('mean_density_times_current_volume_relative',np.max(abs(diag['cell_mean_density_kg_m3']*diag['cell_volume_m3']/prism.mass-1)),3e-14)
    twice=NativePrisms(gradients,weights,h0,axes,2*rho);e2,g2,d2=twice.evaluate(du)
    check('density_only_changes_mass_and_density',energy==e2 and np.array_equal(gradient,g2) and np.array_equal(d2['cell_material_mass_kg'],2*diag['cell_material_mass_kg']) and np.array_equal(d2['quadrature_current_density_kg_m3'],2*diag['quadrature_current_density_kg_m3']))
    small_prism=[]
    single=NativePrisms(gradients[:1],weights[:1],np.zeros_like(h0[:1]),np.eye(3),rho[:1])
    for eps in (1e-6,1e-9,1e-12):
        h=np.zeros((3,3));h[0,1]=eps;u=(X@h.T)[None]
        es,gs,ds=single.evaluate(u);wm,pm,_=law.evaluate_incremental(h)
        eg=np.einsum('qib,ab,q->ia',cells[0].grad,pm,cells[0].weights)
        small_prism.append(dict(eps=eps,energy_J=es,energy_relative_error=abs(es-float(wm)*cells[0].reference_volume)/(float(wm)*cells[0].reference_volume),gradient_relative_error=float(np.linalg.norm(gs[0]-eg)/np.linalg.norm(eg))))
    check('prism_tiny_shear_energy_through_1e_minus_12',max(r['energy_relative_error'] for r in small_prism),2e-12)
    check('prism_tiny_shear_gradient_through_1e_minus_12',max(r['gradient_relative_error'] for r in small_prism),2e-12)
    check('reversible_unloading_returns_exact_reference',single.evaluate(np.zeros((1,18,3)))[0]==0 and diag['dissipated_energy_J']==0)
    reject('prism_reject_negative_reference_weight',lambda:NativePrisms(gradients,-weights,h0,axes,rho))
    reject('prism_reject_zero_density',lambda:NativePrisms(gradients,weights,h0,axes,0.))
    reject('prism_reject_inverted_origin',lambda:NativePrisms(gradients,weights,np.diag([-2.,0,0]),axes,rho))
    reject('prism_reject_large_increment',lambda:prism.evaluate(np.broadcast_to(X,(2,18,3))))
    check('frozen_sources_unchanged',{str(p.relative_to(BASE)):sha(p) for p in FROZEN}==before)
    result=dict(status='EXPERIMENTAL_UNCALIBRATED_NATIVE_PORT',passed=all(checks.values()),checks=checks,metrics=metrics,tiny_strain=small,tiny_prism=small_prism,frozen_dependency_sha256=before,seconds=time.perf_counter()-started,maximum_resident_memory_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,scope='Fixed-law batch material points and isolated 18-node prism contractions only; no fit, growth solve, fracture, scene or render.')
    (ROOT/'receipts/native_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('checks','metrics','tiny_strain','frozen_dependency_sha256')},indent=2))
    print(f'{sum(checks.values())}/{len(checks)} checks passed')
if __name__=='__main__':main()
