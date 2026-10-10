"""Bounded reproducible material-point and one-element checks. No morphology solve."""
from pathlib import Path
import hashlib
import importlib.util
import json
import resource
import time
import numpy as np
from scipy.optimize import brentq
from directional_candidate import Branch, DirectionalEnvelope, ExperimentalPrism, synthetic_candidate
from quadratic_prism_frozen import QuadraticPrism, flat_reference, shapes

ROOT = Path(__file__).resolve().parents[1]
RECEIPTS = ROOT/'receipts'
SCALE = 1e6


def rotation(rng):
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    q[:, 0] *= np.linalg.det(q)
    return q


def main():
    start = time.perf_counter()
    rng = np.random.default_rng(713)
    law = synthetic_candidate()
    X = flat_reference(np.array([[0., 0.], [0.001, 0.], [0., 0.001]]), 0.001)
    cell = ExperimentalPrism(X, 240., law)
    checks = {}
    metrics = {}

    def check(name, value, limit=None):
        ok = bool(value) if limit is None else bool(value <= limit)
        checks[name] = ok
        if limit is not None:
            metrics[name] = {'observed': float(value), 'upper_limit': float(limit)}
        if not ok:
            raise AssertionError((name, value, limit))

    W0, P0, _ = law.evaluate(np.eye(3))
    check('reference_energy_and_stress', max(abs(float(W0)), np.max(abs(P0)))/SCALE, 1e-14)
    E0, g0, d0 = cell.evaluate(X)
    check('reference_element_energy_and_force', max(abs(E0)/(SCALE*cell.reference_volume), np.max(abs(g0))/(SCALE*1e-6)), 1e-13)

    cases = [np.diag([s, 1., 1.]) for s in (0.2, 0.35, 0.6, 0.9, 1.02)]
    cases += [np.diag([1., s, 1.]) for s in (0.8, 1.02, 1.04, 1.10)]
    cases += [np.eye(3)+rng.normal(size=(3, 3))*0.035 for _ in range(12)]
    # Evaluate material gradients away from nonsmooth tangent knots.
    fd_error = 0.
    objectivity = 0.
    material_frame = 0.
    stress_symmetry = 0.
    homogeneous = 0.
    component_energy = 0.
    minimum_energy = 0.
    for F in cases:
        W, P, d = law.evaluate(F)
        minimum_energy = min(minimum_energy, float(W))
        h = 1e-6
        fd = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                D = np.zeros((3, 3)); D[i, j] = h
                fd[i, j] = (law.evaluate(F+D)[0]-law.evaluate(F-D)[0])/(2*h)
        fd_error = max(fd_error, float(np.max(abs(fd-P)))/SCALE)
        Q = rotation(rng)
        Wq, Pq, _ = law.evaluate(Q@F)
        objectivity = max(objectivity, abs(float(Wq-W))/SCALE, float(np.max(abs(Pq-Q@P)))/SCALE)
        R = rotation(rng)
        Wr, Pr, _ = synthetic_candidate(R).evaluate(F@R.T)
        material_frame = max(material_frame, abs(float(Wr-W))/SCALE, float(np.max(abs(Pr-P@R.T)))/SCALE)
        sigma = P@F.T/np.linalg.det(F)
        stress_symmetry = max(stress_symmetry, float(np.max(abs(sigma-sigma.T)))/SCALE)
        en, gr, diag = cell.evaluate(X@F.T)
        Ph = gr.T@X/cell.reference_volume
        homogeneous = max(homogeneous, abs(en/cell.reference_volume-float(W))/SCALE, float(np.max(abs(Ph-P)))/SCALE)
        component_energy = max(component_energy, abs(en-diag['axial_energy_J']-diag['angular_energy_J'])/(SCALE*cell.reference_volume))
    check('material_gradient_relative_error', fd_error, 1e-8)
    check('superposed_rotation_energy_and_P', objectivity, 2e-13)
    check('reference_frame_covariance', material_frame, 2e-13)
    check('Cauchy_stress_symmetry', stress_symmetry, 2e-13)
    check('homogeneous_energy_and_stress', homogeneous, 2e-13)
    check('energy_component_accounting', component_energy, 2e-13)
    check('nonnegative_sampled_energy', minimum_energy >= -1e-9)

    # Non-affine one-element forces, conservation, density and interior orientation.
    F = np.array([[0.75, 0.04, 0.01], [0.0, 1.025, 0.02], [0.01, 0.0, 0.99]])
    x = X@F.T+rng.normal(size=X.shape)*1e-6
    energy, grad, diag = cell.evaluate(x)
    h = 2e-10
    fd = np.zeros_like(x)
    for i in range(18):
        for j in range(3):
            D = np.zeros_like(x); D[i, j] = h
            fd[i, j] = (cell.evaluate(x+D)[0]-cell.evaluate(x-D)[0])/(2*h)
    check('nonaffine_nodal_energy_gradient', float(np.max(abs(fd-grad)))/(SCALE*1e-6), 1e-8)
    check('balanced_internal_force', float(np.linalg.norm(grad.sum(axis=0)))/(SCALE*1e-6), 1e-12)
    check('balanced_internal_moment', float(np.linalg.norm(np.cross(x, grad).sum(axis=0)))/(SCALE*1e-9), 1e-12)
    Q = rotation(rng)
    eq, gq, _ = cell.evaluate(x@Q.T+np.array([0.12, -0.04, 0.07]))
    check('element_rigid_motion_energy_and_force', max(abs(eq-energy)/(SCALE*cell.reference_volume), float(np.max(abs(gq-grad@Q.T)))/(SCALE*1e-6)), 1e-10)
    Fq = np.einsum('ia,qib->qab', x-x.mean(axis=0), cell.grad)
    Jq = np.linalg.det(Fq)
    integrated_mass = np.sum(diag['quadrature_current_density_kg_m3']*Jq*cell.weights)
    check('reference_mass_equals_current_density_integral', abs(integrated_mass-cell.mass)/cell.mass, 1e-13)
    check('mass_matrix_constant_density_integral', abs(cell.mass_scalar.sum()-cell.mass)/cell.mass, 1e-13)
    check('mean_density_is_mass_over_current_volume', abs(diag['mapped_mean_density_kg_m3']*diag['mapped_volume_m3']/cell.mass-1), 1e-13)
    determinants = []
    for ir in range(7):
        for js in range(7-ir):
            for z in np.linspace(-1, 1, 9):
                _, dn = shapes(ir/6, js/6, z)
                Jref = (X-X.mean(axis=0)).T@dn
                Fsample = (x-x.mean(axis=0)).T@dn@np.linalg.inv(Jref)
                determinants.append(np.linalg.det(Fsample))
    check('positive_orientation_at_252_extra_points', min(determinants) > 0)
    metrics['extra_orientation_sample_minimum_det_F'] = float(min(determinants))

    scaled = ExperimentalPrism(X*2, 240., law)
    es, gs, ds = scaled.evaluate(x*2)
    check('SI_length_scaling_energy_volume_mass', max(abs(es/energy-8)/8, abs(ds['mapped_volume_m3']/diag['mapped_volume_m3']-8)/8, abs(scaled.mass/cell.mass-8)/8), 1e-13)
    check('SI_length_scaling_gradient', float(np.max(abs(gs-4*grad)))/(SCALE*1e-6), 1e-12)
    density_only = ExperimentalPrism(X, 480., law)
    ed, gd, dd = density_only.evaluate(x)
    check('reference_density_changes_mass_only', ed == energy and np.array_equal(gd, grad) and density_only.mass == 2*cell.mass and dd['mapped_mean_density_kg_m3'] == 2*diag['mapped_mean_density_kg_m3'])

    bad_cases = [np.diag([-0.8, 1, 1]), np.diag([0, 1, 1]), np.diag([0.1, 1, 1]), np.diag([1, 1.2, 1]), np.full((3, 3), np.nan)]
    rejects = 0
    for bad in bad_cases:
        try:
            law.evaluate(bad)
        except ValueError:
            rejects += 1
    check('inversion_nonfinite_and_domain_rejection', rejects == len(bad_cases))

    # C1 energy/stress across the explicitly piecewise tangent junctions.
    knot_errors = []
    for branch, knot in [(law.compression[0], 0.08), (law.compression[0], 0.55), (law.tension[1], 0.03)]:
        wm, pm, _ = branch.evaluate(knot-1e-9)
        wp, pp, _ = branch.evaluate(knot+1e-9)
        knot_errors.append(max(abs(float(wp-wm)), abs(float(pp-pm)))/SCALE)
    check('C1_branch_junction_continuity', max(knot_errors), 3e-9)
    _, _, low = law.compression[0].evaluate(0.04)
    _, _, plateau = law.compression[0].evaluate(0.3)
    _, _, dense = law.compression[0].evaluate(0.75)
    check('radial_initial_plateau_densification_tangent_order', float(plateau) < float(low) < float(dense))
    _, _, t0 = law.tension[1].evaluate(0.015)
    _, _, t1 = law.tension[1].evaluate(0.05)
    check('tangential_decreasing_but_positive_tangent', 0 < float(t1) < float(t0))
    check('radial_tension_compression_asymmetry', abs(float(law.evaluate(np.diag([1.1, 1, 1]))[1][0, 0])) > abs(float(law.evaluate(np.diag([0.9, 1, 1]))[1][0, 0])))
    check('radial_vs_tangential_compression_distinct', abs(float(law.evaluate(np.diag([0.6, 1, 1]))[1][0, 0])) < abs(float(law.evaluate(np.diag([1, 0.6, 1]))[1][1, 1])))

    # Conservative closed-cycle audit and loading work refinement. No dissipation.
    work_rows = []
    final_energy = float(law.evaluate(np.diag([0.2, 1, 1]))[0])
    for n in (200, 400, 800):
        stretches = np.linspace(1, 0.2, n+1)
        fields = np.repeat(np.eye(3)[None], n+1, axis=0)
        fields[:, 0, 0] = stretches
        _, stresses, _ = law.evaluate(fields)
        loading = np.trapezoid(stresses[:, 0, 0], stretches)
        unloading = np.trapezoid(stresses[::-1, 0, 0], stretches[::-1])
        work_rows.append(dict(increments=n, loading_work_J_m3=float(loading), stored_energy_J_m3=final_energy, relative_gap=float(abs(loading-final_energy)/final_energy), cycle_work_J_m3=float(loading+unloading)))
    check('loading_work_refines_to_stored_energy', work_rows[-1]['relative_gap'] < work_rows[-2]['relative_gap'] < work_rows[0]['relative_gap'])
    check('loading_work_relative_gap_800', work_rows[-1]['relative_gap'], 2e-5)
    check('closed_cycle_zero_work', max(abs(row['cycle_work_J_m3']) for row in work_rows)/SCALE, 1e-13)

    # A small sampled rank-one tangent audit, not a global ellipticity certificate.
    minimum_rank_one = np.inf
    for F in cases:
        for _ in range(60):
            a = rng.normal(size=3); a /= np.linalg.norm(a)
            b = rng.normal(size=3); b /= np.linalg.norm(b)
            D = np.outer(a, b); h = 1e-5
            deltaP = (law.evaluate(F+h*D)[1]-law.evaluate(F-h*D)[1])/(2*h)
            minimum_rank_one = min(minimum_rank_one, float(np.sum(deltaP*D)))
    metrics['sampled_rank_one_tangent_minimum_Pa'] = float(minimum_rank_one)
    checks['sampled_rank_one_tangents_positive'] = bool(minimum_rank_one > 0)

    # Old benchmark preserved byte-for-byte and produces identical results.
    original_path = ROOT.parent/'causal-bark/volume-core/source/quadratic_prism.py'
    frozen_path = ROOT/'source/quadratic_prism_frozen.py'
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    check('old_source_snapshot_hash_exact', sha(original_path) == sha(frozen_path))
    spec = importlib.util.spec_from_file_location('original_quadratic_prism', original_path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    old = QuadraticPrism(X, 2.8e6, 0.28, 240.)
    matched = QuadraticPrism(X, 1e6, 0., 240.)
    orig = module.QuadraticPrism(X, 2.8e6, 0.28, 240.)
    oe, og, od = old.evaluate(x)
    re, rg, rd = orig.evaluate(x)
    check('old_law_exact_output_preserved', oe == re and np.array_equal(og, rg) and od == rd)

    # Same infinitesimal E=1 MPa, nu=0, G=0.5 MPa control avoids a modulus confound.
    linear_limit_errors = []
    for D in (np.diag([1., 0., 0.]), np.diag([0., 1., 0.]), np.array([[0., 1., 0.], [0., 0., 0.], [0., 0., 0.]])):
        h = 1e-6
        _, pp, _ = law.evaluate(np.eye(3)+h*D)
        _, pm, _ = law.evaluate(np.eye(3)-h*D)
        candidate_tangent = (pp-pm)/(2*h)
        reference_tangent = matched.mu*(D+D.T)+matched.lam*np.trace(D)*np.eye(3)
        linear_limit_errors.append(float(np.max(abs(candidate_tangent-reference_tangent)))/SCALE)
    check('matched_neo_hookean_infinitesimal_normal_and_shear_tangent', max(linear_limit_errors), 1e-8)

    rows = []
    for axis, label in enumerate(('radial', 'tangential', 'axial')):
        for strain in (-0.8, -0.7, -0.55, -0.4, -0.08, -0.04, 0., 0.015, 0.03, 0.05, 0.1):
            F = np.eye(3); F[axis, axis] = 1+strain
            w, p, d = law.evaluate(F)
            sigma = p@F.T/float(d['J'])
            free = [i for i in range(3) if i != axis]
            check(f'candidate_free_lateral_{axis}_{strain}', float(np.max(abs(p[free])))/SCALE, 1e-13)
            # Existing benchmark both fixed and traction-free lateral conditions.
            old_results = {}
            for lateral in ('fixed', 'traction_free'):
                a = 1. if lateral == 'fixed' else brentq(lambda a: old.mu*(a*a-1)+old.lam*np.log((1+strain)*a*a), 0.1, 4., xtol=1e-14)
                Fold = np.eye(3)*a; Fold[axis, axis] = 1+strain
                ew, gw, dw = old.evaluate(X@Fold.T)
                Pold = gw.T@X/old.reference_volume
                old_results[lateral] = dict(lateral_stretch=float(a), nominal_stress_Pa=float(Pold[axis, axis]), energy_J_m3=ew/old.reference_volume, relative_density=1/np.linalg.det(Fold))
            me, mg, _ = matched.evaluate(X@F.T)
            mp = mg.T@X/matched.reference_volume
            rows.append(dict(axis=label, engineering_strain=strain, candidate_nominal_stress_Pa=float(p[axis, axis]), candidate_Cauchy_stress_Pa=float(sigma[axis, axis]), candidate_energy_J_m3=float(w), candidate_transverse_stretches=[1., 1.], candidate_relative_density=1/float(d['J']), original_neo_hookean=old_results, matched_small_strain_neo_hookean=dict(E_Pa=1e6, nu=0., nominal_stress_Pa=float(mp[axis, axis]), energy_J_m3=me/matched.reference_volume)))

    RECEIPTS.mkdir(exist_ok=True)
    out = dict(status='EXPERIMENTAL_UNCALIBRATED', passed=bool(all(checks.values())), checks=checks, metrics=metrics, loading_work=work_rows, original_source_sha256=sha(original_path), seconds=time.perf_counter()-start, maximum_resident_memory_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, scope='Material points and one quadratic prism only. No full solve or render. No unloading calibration, fracture or growth simulation.')
    (RECEIPTS/'candidate_audit.json').write_text(json.dumps(out, indent=2)+'\n')
    (RECEIPTS/'uniaxial_comparison.json').write_text(json.dumps(dict(status='synthetic demonstration; numbers are not a cork fit', old_law_parameters=dict(E_Pa=2.8e6, nu=0.28), rows=rows), indent=2)+'\n')
    print(json.dumps({k: v for k, v in out.items() if k not in ('checks', 'metrics')}, indent=2))
    print(f'{sum(checks.values())}/{len(checks)} checks passed; rank-one sample minimum = {minimum_rank_one:.6g} Pa')
    if not out['passed']:
        raise AssertionError('Audit includes a failed gate')


if __name__ == '__main__':
    main()
