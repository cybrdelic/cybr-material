"""Independent identities before comparing the notched FE quantities."""
from pathlib import Path
import json
import resource
import sys
import time
import unittest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'source'))
from notch_reference import Model, element_rule, relative, shape


class ReferenceTests(unittest.TestCase):
    def test_partition_and_orientation(self):
        N, dN = shape(0.17, -0.39)
        self.assertAlmostEqual(N.sum(), 1)
        np.testing.assert_allclose(dN.sum(axis=0), 0, atol=1e-15)
        corners = np.array([[0, 0], [.002, 0], [.002, .003], [0, .003]])
        self.assertGreater(np.linalg.det(corners.T @ dN), 0)

    def test_element_energy_force_rigid_modes(self):
        hx, hy, t = .002, .003, .001
        D = np.diag([5.6e6, 2.8e6, .7e6])
        ke, rule = element_rule(hx, hy, D, t)
        xy = np.array([[0, 0], [hx, 0], [hx, hy], [0, hy]])
        L = np.array([[.003, .002], [-.001, -.0004]])
        u = (xy @ L.T).ravel()
        eps = np.array([L[0, 0], L[1, 1], L[0, 1]+L[1, 0]])
        exact_energy = .5*eps @ D @ eps*hx*hy*t
        self.assertLess(relative(.5*u @ ke @ u, exact_energy), 1e-12)
        # Analytical nodal energy gradient vs independently perturbed energies.
        direction = np.array([1., -2., 3., -4., 4., 2., -1., 3.])
        perturb = 1e-8
        fd = (.5*(u+perturb*direction) @ ke @ (u+perturb*direction)
              -.5*(u-perturb*direction) @ ke @ (u-perturb*direction))/(2*perturb)
        self.assertLess(relative(fd, direction @ (ke @ u)), 1e-10)
        for rigid in [np.tile([1., 0], 4), np.tile([0., 1.], 4),
                      np.column_stack([-xy[:, 1], xy[:, 0]]).ravel()]:
            self.assertLess(np.linalg.norm(ke @ rigid)/(np.linalg.norm(ke)*np.linalg.norm(rigid)), 1e-14)
        np.testing.assert_allclose(ke, ke.T, rtol=1e-14, atol=1e-12)

    def test_intact_coupon_analytic_solution(self):
        for ex, ey, shear in [(2.8e6, 2.8e6, 1.4e6), (5.6e6, 2.8e6, .7e6)]:
            m = Model(8, .04, .04, .001, np.diag([ex, ey, shear]))
            s = m.solve(0, 4e-6)
            expected_P = ey*4e-6/.04*.04*.001
            self.assertLess(relative(s['P_N'], expected_P), 1e-12)
            np.testing.assert_allclose(s['u'][1::2], m.xy[:, 1]*4e-6/.04, atol=1e-17)
            self.assertLess(s['relative_work_error'], 1e-12)

    def test_thickness_and_quadratic_load_scaling(self):
        values = []
        for thickness, delta in [(.001, 4e-6), (.003, 4e-6), (.001, 8e-6)]:
            m = Model(16, .04, .04, thickness, np.diag([2.8e6, 2.8e6, 1.4e6]))
            s = m.solve(8, delta)
            values.append((s, m.domain_J(s, [.01])[0]))
        (s, j), (s3, j3), (s2, j2) = values
        self.assertLess(relative(s3['P_N'], 3*s['P_N']), 1e-12)
        self.assertLess(relative(s3['energy_full_J'], 3*s['energy_full_J']), 1e-12)
        self.assertLess(relative(j3, j), 1e-12)
        self.assertLess(relative(s2['P_N'], 2*s['P_N']), 1e-12)
        self.assertLess(relative(j2, 4*j), 1e-12)

    def test_analytic_mode_I_domain_sign_and_factor(self):
        # Exact isotropic nu=0 crack displacement field, polar quadrature.
        # Independent from the FE solve and nodal q implementation.
        E, KI, kappa = 2.8e6, 100.0, 3.0
        D = np.diag([E, E, E/2])
        gauss, weights = np.polynomial.legendre.leggauss(80)
        theta = (gauss+1)*np.pi/2  # Upper half [0, pi].
        wt = weights*np.pi/2
        for rin, rout in [(.003, .006), (.005, .01), (.007, .014)]:
            radii = (gauss+1)*(rout-rin)/2+rin
            wr = weights*(rout-rin)/2
            r, th = np.meshgrid(radii, theta, indexing='ij')
            c, s, ch, sh = np.cos(th), np.sin(th), np.cos(th/2), np.sin(th/2)
            A = KI/E/np.sqrt(2*np.pi)
            common = kappa-c
            ur = A/(2*np.sqrt(r))[..., None]*np.stack([ch*common, sh*common], axis=-1)
            uth = A*np.sqrt(r)[..., None]*np.stack([-.5*sh*common+ch*s, .5*ch*common+sh*s], axis=-1)
            ux = ur*c[..., None]-uth*s[..., None]/r[..., None]
            uy = ur*s[..., None]+uth*c[..., None]/r[..., None]
            eps = np.stack([ux[..., 0], uy[..., 1], uy[..., 0]+ux[..., 1]], axis=-1)
            sig = eps @ D
            W = .5*np.sum(eps*sig, axis=-1)
            A1 = sig[..., 0]*ux[..., 0]+sig[..., 2]*ux[..., 1]-W
            A2 = sig[..., 2]*ux[..., 0]+sig[..., 1]*ux[..., 1]
            integrand = -(A1*c+A2*s)/(rout-rin)*r
            J = 2*np.einsum('i,ij,j->', wr, integrand, wt)
            self.assertLess(relative(J, KI**2/E), 1e-12)

    def test_FE_domain_implementation_against_exact_crack_field(self):
        # This feeds the ACTUAL domain_J routine an independently known field.
        # It includes Q4 interpolation error, so use the predeclared 3% J gate,
        # not the algebraic-identity tolerance.
        cfg = json.loads((Path(__file__).resolve().parents[1]/'PROTOCOL.json').read_text())
        E, KI, a = 2.8e6, 100., .02
        m = Model(64, .04, .04, .001, np.diag([E, E, E/2]))
        dx, dy = m.xy[:, 0]-a, m.xy[:, 1]
        r, theta = np.hypot(dx, dy), np.arctan2(dy, dx)
        common = KI/E*np.sqrt(r/(2*np.pi))*(3-np.cos(theta))
        u = np.column_stack([common*np.cos(theta/2), common*np.sin(theta/2)]).ravel()
        Js = m.domain_J({'u': u, 'a_m': a}, cfg['domain_outer_radii_m'])
        for J in Js:
            self.assertLess(relative(J, KI**2/E), cfg['gates']['finest_J_vs_fixed_displacement_energy_relative'])


if __name__ == '__main__':
    start = time.monotonic()
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(ReferenceTests))
    receipt = {'tests_run': result.testsRun, 'failures': len(result.failures),
               'errors': len(result.errors), 'passed': result.wasSuccessful(),
               'wall_seconds': time.monotonic()-start,
               'maximum_resident_memory_KiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
               'identities': ['Q4 partition and positive orientation', 'affine energy and analytic force gradient',
                              'rigid modes and stiffness symmetry', 'intact-coupon closed-form force and displacement',
                              'thickness and quadratic load scaling', 'exact Mode I polar domain integral',
                              'actual FE domain routine evaluated on exact crack-field nodal data']}
    (Path(__file__).resolve().parents[1]/'receipts'/'test_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    raise SystemExit(0 if result.wasSuccessful() else 1)
