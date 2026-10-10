"""Read-only source/parameter audit plus explicitly conditional dimensional arithmetic.
No material fitting, finite-element solve, parameter modification or image output.
"""
from pathlib import Path
import hashlib
import json
import math
import resource
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
WORK = ROOT.parent


def main():
    started = time.perf_counter()
    frozen = json.loads((ROOT/'receipts/final_manifest.json').read_text())
    verified = {f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() == v['sha256']
                for f, v in frozen['files'].items()}
    assert all(verified.values()), 'The original 66-test candidate must remain frozen'
    old_path = WORK/'bark-panel-resume/receipts/cohesive_resolution_audit.json'
    old = json.loads(old_path.read_text())
    old_rows = []
    for p in old['layers']:
        E, G, T = p['E_Pa'], p['Gc_J_per_m2'], p['cohesive_strength_Pa']
        ell = E*G/T**2
        assert math.isclose(ell, p['energy_strength_characteristic_length_m'], rel_tol=1e-14)
        old_rows.append(dict(layer=p['layer'], E_Pa=E, Gc_J_m2=G, peak_traction_Pa=T,
            isotropic_plane_stress_K_Pa_sqrt_m=math.sqrt(E*G),
            isotropic_plane_strain_K_Pa_sqrt_m=math.sqrt(E/(1-0.28**2)*G),
            nominal_characteristic_length_m=ell,
            bilinear_final_opening_m=2*G/T,
            old_14mm_pitch_over_length=0.014/ell,
            current_0_525mm_pitch_over_length=0.000525/ell))

    # Primary 2010 Table 2 group means only. These are not raw paired specimens.
    tensile_means = [
        ('class1_inner',23.72,0.75), ('class1_mid',23.88,0.62),
        ('class1_outer',18.62,0.50), ('class4_inner',24.31,0.70),
        ('class4_mid',19.84,0.56), ('class4_outer',20.86,0.57)]
    K_nonradial = 94e3
    K_radial = 125e3
    K_nr_display_endpoints = (78e3, 110e3)
    K_r_display_endpoints = (111e3, 139e3)
    conditional_rows = []
    for group, E_MPa, f_MPa in tensile_means:
        Eproxy, Tproxy = E_MPa*1e6, f_MPa*1e6
        Gproxy = K_nonradial**2/Eproxy
        ellproxy = Eproxy*Gproxy/Tproxy**2
        assert math.isclose(ellproxy, (K_nonradial/Tproxy)**2, rel_tol=1e-14)
        conditional_rows.append(dict(group=group, measured_tangential_E_mean_Pa=Eproxy,
            measured_nominal_fracture_stress_mean_Pa=Tproxy,
            assumed_effective_crack_modulus_Pa=Eproxy,
            assumed_cohesive_peak_traction_Pa=Tproxy,
            mixed_study_G_proxy_J_m2=Gproxy,
            mixed_study_length_proxy_m=ellproxy,
            hypothetical_bilinear_final_opening_m=2*Gproxy/Tproxy))

    E_range = [min(v[1] for v in tensile_means)*1e6,max(v[1] for v in tensile_means)*1e6]
    T_range = [min(v[2] for v in tensile_means)*1e6,max(v[2] for v in tensile_means)*1e6]
    conditional_ranges = dict(
        status='Arithmetic envelope of incompatible study/group means and displayed K +/- endpoints; NOT a statistical or constitutive bound',
        assumptions=['1991 nonradial K applies to 2010 tangential samples',
                     'effective anisotropic crack modulus is replaced by measured uniaxial E',
                     'cohesive peak traction is replaced by unnotched nominal fracture stress',
                     'length prefactor is one and small-scale/linear-elastic assumptions apply'],
        G_proxy_J_m2=[K_nr_display_endpoints[0]**2/E_range[1],K_nr_display_endpoints[1]**2/E_range[0]],
        length_proxy_m=[(K_nr_display_endpoints[0]/T_range[1])**2,(K_nr_display_endpoints[1]/T_range[0])**2],
        mean_E_over_old_E_extreme_ratios=[E_range[0]/2.8e6,E_range[1]/1.8e6],
        nominal_strength_mean_over_old_peak_extreme_ratios=[T_range[0]/180e3,T_range[1]/90e3],
        central_nonradial_K_with_old_strength_length_m=[(K_nonradial/180e3)**2,(K_nonradial/90e3)**2])
    sensitivities = []
    for direction, K in [('nonradial',K_nonradial),('radial',K_radial)]:
        for Eeffective in (10e6,20e6,40e6):
            sensitivities.append(dict(direction=direction,assumed_effective_crack_modulus_Pa=Eeffective,
                                     conditional_G_J_m2=K*K/Eeffective))
    out = dict(
        status='EVIDENCE_AND_CONDITIONAL_SCALE_ONLY_NO_PARAMETER_CHANGE',
        frozen_candidate_hash_checks=verified,
        sources={
            'toughness_1991': dict(url='https://link.springer.com/article/10.1007/BF00576525',
                access='Official publisher abstract verified 2026-10-06; full specimen methods not accessed',
                K_nonradial_Pa_sqrt_m=K_nonradial,K_nonradial_reported_plus_minus_Pa_sqrt_m=16e3,
                K_radial_Pa_sqrt_m=K_radial,K_radial_reported_plus_minus_Pa_sqrt_m=14e3,
                statistical_meaning='Abstract does not define the +/- statistic; endpoints are not confidence bounds',
                orientation='Three independent mode-I systems grouped into two mechanism/value classes; full crack-plane/front mapping needed'),
            'tension_2010': dict(url='https://www.researchgate.net/publication/233765103_Tensile_properties_of_cork_in_the_tangential_direction_Variation_with_quality_porosity_density_and_radial_position_in_the_cork_plank',
                locator='Table 2; group means for twelve samples, with SD reported separately',
                caveat='Body class-4 overall E=18.3 MPa does not equal the simple mean of displayed Table-2 subgroup means; do not reconcile without original data',
                specimen='Boiled/air-dried processed planks, about 7% mean moisture; not age-indexed bark cohorts'),
            'cohesive_scale_method': dict(url='https://ntrs.nasa.gov/api/citations/20080020385/downloads/20080020385.pdf',
                locator='Song, Davila and Rose, Guidelines and Parameter Selection for the Simulation of Progressive Delamination, equations 3-4; numerical methodology, not cork data'),
            'old_parameter_receipt': dict(path=str(old_path.relative_to(WORK)),sha256=hashlib.sha256(old_path.read_bytes()).hexdigest()),
            'old_parameter_source': dict(path='bark-panel-resume/source/small_coupon_r1.py',
                sha256=hashlib.sha256((WORK/'bark-panel-resume/source/small_coupon_r1.py').read_bytes()).hexdigest())},
        formulas=dict(anisotropic_mode_I='G_I = H_I(C, crack normal, propagation direction, front, plane state) K_I^2; define E_I_eff=1/H_I; H_I is unknown here',
                      characteristic_length='ell_proxy = E_I_eff G_Ic / T0^2; actual l_cz = shape/material/geometry factor times ell_proxy',
                      mixed_study_proxy='ell_proxy=(K_Ic/T0)^2 only under the declared same-effective-modulus and strength substitutions',
                      bilinear_law='delta_f=2 G_Ic/T0, delta_0=T0/K_pen; K_pen is Pa/m, distinct from toughness K_Ic'),
        old_isotropic_benchmark=old_rows,
        nonradial_cross_study_conditional_cases=conditional_rows,
        conditional_arithmetic_envelope=conditional_ranges,
        effective_crack_modulus_scenarios=dict(status='Hypothetical sensitivities only; chosen 10/20/40 MPa values are NOT measured anisotropic crack-modulus bounds',rows=sensitivities),
        recommendation='Retain old numerical refinement conclusions for the old law. Use a single homogeneous, single-notched direction-specific coupon to identify compliance/energy and the actual traction-zone length before assigning a physical mesh target. Do not swap Gc or reduce peak strength for speed.',
        seconds=time.perf_counter()-started,
        maximum_resident_memory_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (HERE/'fracture_scale_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(frozen_files_verified=len(verified),old_length_mm=[r['nominal_characteristic_length_m']*1000 for r in old_rows],conditional_nonradial_length_mm=[v*1000 for v in conditional_ranges['length_proxy_m']],conditional_G_J_m2=conditional_ranges['G_proxy_J_m2'],seconds=out['seconds']),indent=2))


if __name__ == '__main__':
    main()
