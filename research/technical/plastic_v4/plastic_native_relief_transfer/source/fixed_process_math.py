"""Small, independently testable corrections to the existing production math.

These functions do not add a material law or calibrate a parameter.  They repair
the transfer/accounting of a law already declared in the frozen candidates.
"""
import numpy as np


def archard_depth(pressure_sliding_Pa_m, coefficient, hardness_Pa):
    """d = k/H * integral(p ds), in metres; no second independent grit budget."""
    work = np.asarray(pressure_sliding_Pa_m, dtype=np.float64)
    if hardness_Pa <= 0 or coefficient < 0 or np.any(work < 0):
        raise ValueError('Nonnegative work/coefficient and positive hardness required')
    return work * (coefficient / hardness_Pa)


def finite_height_gradient(height_m, pixel_pitch_m):
    """Second-order derivatives on a finite coupon, never an opposite-edge wrap.

    Image rows are y-down; callers own the documented image-to-UV axis mapping.
    Periodic stock textures may still use periodic derivatives independently.
    """
    h = np.asarray(height_m, dtype=np.float64)
    if h.ndim != 2 or min(h.shape) < 3 or pixel_pitch_m <= 0:
        raise ValueError('Need a 2D field of at least 3x3 and positive SI pitch')
    gy, gx = np.gradient(h, pixel_pitch_m, edge_order=2)
    return gx, gy


def petg_bead_relief_m(phase, bead_width_m=.00045):
    """Exact selected r5 wall profile, excluding its separate start/stop seam."""
    phase = np.asarray(phase, dtype=np.float64)
    bead = np.sqrt(np.maximum(0., 1. - ((phase - .5) / .52) ** 2))
    return (bead - .65) * bead_width_m * .155


def encode_petg_relief(phase, height_scale_m=.00016, bead_width_m=.00045):
    """Alternative flat-surface relief, exactly sharing the resolved bead state.

    No old independent .008*fine height is added: it is absent from the selected
    mesh. A future unresolved component needs its own representation contract.
    """
    return .5 + petg_bead_relief_m(phase, bead_width_m) / height_scale_m


def mean_film_thickness(nominal_m, deposit_mean_m):
    """Retain the zero mode when top and substrate deviations are centered."""
    if nominal_m <= 0 or deposit_mean_m < 0:
        raise ValueError('Expected positive film and nonnegative deposited volume')
    return nominal_m + deposit_mean_m


def elastic_pressure_transition(s_m, tool_depth_m, pressure0_Pa, pressure1_Pa,
                                layer_m, modulus_Pa):
    """Quasistatic linear pressure ramp at fixed viscous state (no dashpot time).

    The elastic reservoir changes; pressure work is its exact difference.  A
    truly instantaneous dynamic load would require inertia/damping instead.
    This routine explicitly chooses the same quasistatic model as the core.
    """
    s = np.asarray(s_m, dtype=np.float64)
    H = np.asarray(tool_depth_m, dtype=np.float64)
    if min(pressure0_Pa, pressure1_Pa) < 0 or layer_m <= 0 or modulus_Pa <= 0:
        raise ValueError('Invalid pressure or elastic parameters')
    u0 = np.minimum(H, s + pressure0_Pa * layer_m / modulus_Pa)
    u1 = np.minimum(H, s + pressure1_Pa * layer_m / modulus_Pa)
    e0 = .5 * modulus_Pa / layer_m * (u0 - s) ** 2
    e1 = .5 * modulus_Pa / layer_m * (u1 - s) ** 2
    return dict(displacement_m=u1, pressure_work_J_m2=e1-e0,
                elastic_energy_before_J_m2=e0, elastic_energy_after_J_m2=e1,
                dissipated_energy_J_m2=np.zeros_like(s))
