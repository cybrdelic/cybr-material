"""EXPERIMENTAL, UNCALIBRATED monotonic-envelope material-point candidate.

Reference energy W(F,A) = sum_i phi_i(|F a_i|-1)
                         - G log(det(F) / product_i |F a_i|).
A has the orthonormal reference radial, tangential, axial directions as columns.
No cell simulation, growth, fracture, damage, plasticity or hysteresis is implied.
The standalone prism adapter reuses an unchanged numerical benchmark snapshot.
"""
from dataclasses import dataclass
import numpy as np
from quadratic_prism_frozen import QuadraticPrism


@dataclass(frozen=True)
class Branch:
    """Positive stress magnitude versus positive engineering strain magnitude.

    Initial modulus E; tangent H after y. Optional compression-only densification
    D*(s-d)^2/(1-s) after d. All coefficients are Pa; y,d,s are dimensionless.
    """
    E: float
    H: float
    y: float
    D: float = 0.0
    d: float = 0.8

    def __post_init__(self):
        if not np.isfinite([self.E, self.H, self.y, self.D, self.d]).all():
            raise ValueError("Nonfinite branch parameter")
        if self.E <= 0 or self.H < 0 or self.y <= 0 or self.D < 0 or not 0 < self.d < 1:
            raise ValueError("Invalid branch parameter")
        if self.D and self.y >= self.d:
            raise ValueError("Densification must follow the initial branch")

    def evaluate(self, s):
        s = np.asarray(s, dtype=float)
        if np.any(s < 0) or not np.isfinite(s).all():
            raise ValueError("Branch needs nonnegative finite strain magnitude")
        first = np.minimum(s, self.y)
        post = np.maximum(s - self.y, 0)
        energy = 0.5*self.E*first**2 + self.E*self.y*post + 0.5*self.H*post**2
        stress = self.E*first + self.H*post
        tangent = np.where(s < self.y, self.E, self.H)
        if self.D:
            if np.any(s >= 1):
                raise ValueError("Densification requires compression magnitude < 1")
            t = np.maximum(s-self.d, 0)
            a = 1-self.d
            u = t/a
            # -log(1-u)-u-u^2/2, avoiding cancellation at onset.
            remainder = -np.log1p(-u)-u-0.5*u*u
            small = u < 0.001
            series = sum(u**k/k for k in range(3, 11))
            remainder = np.where(small, series, remainder)
            energy = energy + self.D*a*a*remainder
            stress = stress + self.D*t*t/(a-t)
            tangent = tangent + self.D*t*(2*a-t)/(a-t)**2
        return energy, stress, tangent


class DirectionalEnvelope:
    def __init__(self, tension, compression, shear_Pa, axes=None,
                 minimum_stretch=0.15, maximum_stretch=1.15):
        self.tension = tuple(tension)
        self.compression = tuple(compression)
        self.G = float(shear_Pa)
        self.axes = np.eye(3) if axes is None else np.asarray(axes, dtype=float).copy()
        if len(self.tension) != 3 or len(self.compression) != 3:
            raise ValueError("Need radial, tangential and axial branches")
        if not np.isfinite(self.G) or self.G <= 0:
            raise ValueError("Shear modulus must be positive")
        if self.axes.shape != (3, 3) or not np.isfinite(self.axes).all():
            raise ValueError("Invalid material frame")
        if not np.allclose(self.axes.T@self.axes, np.eye(3), atol=1e-12, rtol=0):
            raise ValueError("Material frame must be orthonormal")
        if np.linalg.det(self.axes) < 0:
            raise ValueError("Material frame must be right handed")
        if not np.isfinite([minimum_stretch, maximum_stretch]).all() or not 0 < minimum_stretch < 1 < maximum_stretch:
            raise ValueError("Invalid operational stretch bounds")
        if any(b.D for b in self.tension):
            raise ValueError("Compression densification is not a tensile law")
        self.minimum_stretch = float(minimum_stretch)
        self.maximum_stretch = float(maximum_stretch)
        self.axes.setflags(write=False)

    def evaluate(self, F):
        F = np.asarray(F, dtype=float)
        if F.shape[-2:] != (3, 3) or not np.isfinite(F).all():
            raise ValueError("Expected finite 3x3 deformation gradients")
        J = np.linalg.det(F)
        if np.any(J <= 0):
            raise ValueError("Inverted or singular deformation")
        B = F@self.axes
        stretch = np.linalg.norm(B, axis=-2)
        if np.any(stretch < self.minimum_stretch) or np.any(stretch > self.maximum_stretch):
            raise ValueError("Outside declared experimental stretch domain")
        e = stretch-1
        axial_energy = np.zeros_like(J)
        stress = np.zeros_like(e)
        for i in range(3):
            ut, pt, _ = self.tension[i].evaluate(np.maximum(e[..., i], 0))
            uc, pc, _ = self.compression[i].evaluate(np.maximum(-e[..., i], 0))
            axial_energy += ut+uc
            stress[..., i] = pt-pc
        # Hadamard's inequality makes the angular energy nonnegative for J>0.
        log_angle_volume = np.linalg.slogdet(F)[1]-np.log(stretch).sum(axis=-1)
        angular_energy = -self.G*log_angle_volume
        invT = np.linalg.inv(F).swapaxes(-2, -1)
        coefficient = stress/stretch + self.G/stretch**2
        P = (B*coefficient[..., None, :])@self.axes.T-self.G*invT
        W = axial_energy+angular_energy
        return W, P, dict(J=J, stretches=stretch,
                         axial_energy_J_m3=axial_energy,
                         angular_energy_J_m3=angular_energy)


def synthetic_candidate(axes=None):
    """Numerical fixture, deliberately NOT a fit to published summary values."""
    linear = Branch(1e6, 1e6, 0.2)
    radial_compression = Branch(1e6, 2e4, 0.08, D=2e6, d=0.55)
    tangential_tension = Branch(1e6, 2.5e5, 0.03)
    return DirectionalEnvelope(
        tension=(linear, tangential_tension, linear),
        compression=(radial_compression, linear, linear),
        shear_Pa=0.5e6, axes=axes)


class ExperimentalPrism(QuadraticPrism):
    """Separate adapter; no existing production caller is changed."""
    def __init__(self, reference_nodes, density, law, nq=5):
        # Inherited E,nu only construct unchanged geometry/mass machinery.
        # Their values do not enter the overridden candidate energy/gradient.
        super().__init__(reference_nodes, E=1.0, nu=0.0, density=density, nq=nq)
        self.law = law

    def evaluate(self, nodes):
        x = np.asarray(nodes, dtype=float).reshape(18, 3)
        F = np.einsum('ia,qib->qab', x-x.mean(axis=0), self.grad)
        W, P, d = self.law.evaluate(F)
        J = d['J']
        volume = float(self.weights@J)
        gradient = np.einsum('qib,qab,q->ia', self.grad, P, self.weights)
        return float(self.weights@W), gradient, {
            'mapped_volume_m3': volume,
            'minimum_det_F': float(J.min()),
            'maximum_det_F': float(J.max()),
            'material_mass_kg': self.mass,
            'mapped_mean_density_kg_m3': self.mass/volume,
            'quadrature_current_density_kg_m3': self.rho/J,
            'axial_energy_J': float(self.weights@d['axial_energy_J_m3']),
            'angular_energy_J': float(self.weights@d['angular_energy_J_m3']),
            'dissipated_energy_J': 0.0,
        }
