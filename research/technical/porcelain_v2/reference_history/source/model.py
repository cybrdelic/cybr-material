"""Fixed, explicitly imposed Q92 reference history, in SI units. No look fitting."""
from pathlib import Path
import hashlib
import json
import math
import numpy as np
from scipy.integrate import quad, simpson

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
ARRIVALS = ROOT / 'experiments/porcelain_deposition_arrivals'
GAMMAS = (0.15, 0.3, 0.6)
START_C, END_C, RAMP_K_PER_S = 1020.0, 1040.0, 5.0 / 60.0
R0 = 0.17984575

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def eta(T_C):
    """Published log10(eta_Poise)=-.299+2512/(T_K-789.2); 1 P=.1 Pa s."""
    return 10.0 ** (-1.299 + 2512.0 / (np.asarray(T_C) + 273.15 - 789.2))

def finite_depth(x):
    """Stable no-slip finite-depth Stokes factor; x=angular wave number*depth."""
    x = np.asarray(x, dtype=np.float64)
    result = np.empty_like(x)
    small = x < 1e-3
    z = x[small]
    result[small] = (2*z**3/3)*(1 - 9*z*z/5)
    z = x[~small]
    q = np.exp(-2*z)
    result[~small] = (1-q*q-4*z*q)/(1+q*q+(4*z*z+2)*q)
    return result

def bvp_rate(x):
    sh, ch = math.sinh(x), math.cosh(x)
    a, b = np.linalg.solve([[-2*x*ch, 2*ch+2*x*sh],
                           [2*ch-2*x*sh, 2*x*ch]], [0., -1.])
    return -(a*(sh-x*ch)+b*x*sh)

def setup():
    original = json.loads((ARRIVALS / 'receipts/generation.json').read_text())
    h = original['retained_mean_film_m']
    p = original['parameters']
    duration = (END_C - START_C) / RAMP_K_PER_S
    J, err = quad(lambda t: 1/float(eta(START_C+RAMP_K_PER_S*t)), 0, duration,
                  epsabs=1e-14, epsrel=1e-12)
    refinement = []
    for intervals in (16, 32, 64, 128, 256):
        t = np.linspace(0, duration, intervals+1)
        v = float(simpson(1/eta(START_C+RAMP_K_PER_S*t), x=t))
        refinement.append({'time_intervals': intervals, 'J_per_Pa': v,
                           'relative_error_vs_adaptive': abs(v/J-1)})
    assert refinement[-1]['relative_error_vs_adaptive'] < 1e-10
    xs = [0.1, 0.3, 1., 2., 3.3316292189473677, 5.]
    bvp_errors = [abs(bvp_rate(x)/(float(finite_depth(x))/2)-1) for x in xs]
    assert max(bvp_errors) < 1e-10
    small_error = abs(float(finite_depth(1e-4))/(2e-12/3)-1)
    assert small_error < 2e-8 and float(finite_depth(30)) == 1.0
    assert float(finite_depth(0)) == 0
    rows = []
    for wavelength in (p['tile_m'], p['deposition_cell_m'], p['pressed_powder_cell_m']):
        k = 2*np.pi/wavelength
        for gamma in GAMMAS:
            exponent = gamma*J*k*float(finite_depth(k*h))/2
            rows.append({'wavelength_m': wavelength, 'gamma_N_per_m': gamma,
                         'kh': k*h, 'decay_exponent': exponent,
                         'amplitude_survival': math.exp(-exponent)})
    sample_slope = np.array([0., .1, .5, 1., 2., 3.])
    rough_audit = []
    for amplitude in (1., 1e-3, 1e-6, 1e-9, 1e-12, 1e-15, 0.):
        slope = amplitude*sample_slope
        r = np.clip(R0+.008*(slope/(slope.mean()+1e-12)-1), R0-.008, R0+.016)
        rough_audit.append({'height_amplitude_multiplier': amplitude,
                            'roughness_min': float(r.min()), 'roughness_max': float(r.max()),
                            'roughness_std': float(r.std()), 'roughness_mean': float(r.mean())})
    inputs = [ARRIVALS/'source/generate_control.py', ARRIVALS/'state/arrival_packets.npz',
              ARRIVALS/'receipts/generation.json', ARRIVALS/'maps/Substrate_Height.png',
              ROOT/'experiments/porcelain_finite_depth_audit/audit.py',
              ROOT/'experiments/porcelain_finite_depth_audit/results.json',
              ROOT/'experiments/porcelain_finite_depth_audit/commercial_glaze_2015_addendum.md',
              ROOT/'experiments/porcelain_glaze_layer_assignment/frozen_handoff.json']
    return {
        'experiment': 'porcelain_reference_history', 'native_resolution': p['resolution'],
        'classification': 'Representative imposed commercial-glaze process; not calibrated porcelain chemistry or firing',
        'rheology': {'material': 'Q92 commercial alkaline glaze, CRECER POLES S.A.',
            'source': 'https://doi.org/10.1016/j.mspro.2015.05.031',
            'source_type': 'HSM/VFT fit, not independent rheometry',
            'equation_SI': 'eta_Pa_s=10**(-1.299+2512/(T_C+273.15-789.2))',
            'A_in_Poise': -.299, 'B_K': 2512., 'T0_K': 789.2,
            'temperature_C': [START_C, END_C], 'eta_range_Pa_s': [float(eta(END_C)), float(eta(START_C))]},
        'history': {'start_C': START_C, 'end_C': END_C, 'ramp_K_per_s': RAMP_K_PER_S,
                    'duration_s': duration, 'isothermal_dwell_s': 0,
                    'scope': 'Only the documented 5 C/min ramp through the manufacturer-stated 1020-1040 C window; earlier heating and cooling excluded',
                    'J_per_Pa': J, 'adaptive_absolute_error_estimate_per_Pa': err,
                    'time_quadrature_refinement': refinement},
        'surface_tension': {'central_N_per_m': .3, 'sensitivity_N_per_m': list(GAMMAS),
            'status': 'Glass-literature PROXY, not measured Q92',
            'reference': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC4373153/',
            'envelope_status': 'Declared half/double parameter sensitivity, not calibrated uncertainty interval'},
        'exposure_m': {str(gamma): gamma*J for gamma in GAMMAS},
        'retained_mean_film_m': h, 'scalar_survival': rows,
        'checks': {'independent_BVP_max_relative_error': max(bvp_errors),
                   'small_kh_limit_relative_error': small_error, 'large_kh_limit': True,
                   'zero_mode_transfer_exactly_one': True, 'time_quadrature_refinement_pass': True},
        'optical_accounting_change': {'status': 'HELD separate variant; first history binding retains exact old roughness image', 'intrinsic_source_roughness': R0,
            'new_formula': 'roughness=source_roughness, with no extra resolved-slope variance',
            'old_formula': 'clip(R0+.008*(slope/(mean(slope)+1e-12)-1), R0-.008, R0+.016)',
            'practical_defect': 'Normalized-slope contrast is amplitude-invariant while mean(slope)>>1e-12; decreasing resolved height does not suppress it in that regime.',
            'strict_limit': 'Because of +1e-12, variation tends to zero as amplitude tends exactly to zero, but at R0-.008 rather than R0.',
            'audit': rough_audit, 'separate_from_process_change': True,
            'interpretation': 'R0 retained as unresolved intrinsic finish; resolved normal variance is represented once by the native normal field. No estimated subpixel variance is added.'},
        'assumptions': ['Constant mean film depth and volume; no evaporation, reaction, density or thickness change.',
            'Homogeneous Newtonian connected liquid throughout this imposed interval; the paper does not establish this in the retained porcelain.',
            'Linear small-slope finite-depth Stokes mobility for a flat no-slip support, applied to the retained topography over a corrugated substrate.',
            'Negligible inertia, gravity and viscoelasticity; constant proxy gamma.',
            'No curved-film edge flow, glaze-body chemistry, calibrated spectral optics, or inference that Q92 is a porcelain glaze.',
            'Template geometry and retained optical constants are inherited approximations.',
            'No look-selected coefficients or thermal hold; sensitivity outcomes never select the central gamma.'],
        'input_sha256': {str(q.relative_to(ROOT)): sha(q) for q in inputs},
        'selected': False, 'visual_acceptance': False, 'rendered': False}

if __name__ == '__main__':
    report = setup()
    (HERE/'receipts/setup.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'J_per_Pa': report['history']['J_per_Pa'], 'exposure_m': report['exposure_m'],
                      'checks': report['checks'], 'scalar_survival': report['scalar_survival']}, indent=2))
