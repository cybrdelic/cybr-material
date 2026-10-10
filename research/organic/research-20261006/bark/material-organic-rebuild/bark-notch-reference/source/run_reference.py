"""Bounded experiment; protocol/gates are read from a frozen pre-run file."""
from pathlib import Path
import hashlib
import json
import platform
import resource
import time
import numpy as np
import scipy
from notch_reference import Model, relative, serial_state

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    started = time.monotonic()
    cfg = json.loads((ROOT/'PROTOCOL.json').read_text())
    geo, gates = cfg['geometry_SI'], cfg['gates']
    material_results = {}
    for name, mat in cfg['materials_Pa'].items():
        records = []
        for n in cfg['meshes_nx_ny_half']:
            model = Model(n, geo['width_m'], geo['full_height_m']/2,
                          geo['thickness_m'], np.diag([mat['Ex'], mat['Ey'], mat['Gxy']]))
            tip = round(geo['crack_length_m']/model.hx)
            offsets = [-1, 0, 1]
            if n == cfg['meshes_nx_ny_half'][-1]:
                offsets = [-4, -2, -1, 0, 1, 2, 4]
            states = {i: model.solve(tip+i, geo['half_top_displacement_m']) for i in offsets}
            base = states[0]
            J = model.domain_J(base, cfg['domain_outer_radii_m'], cfg['domain_inner_radius_fraction'])
            diffs = []
            for span in [1, 2, 4] if len(offsets) == 7 else [1]:
                lo, hi = states[-span], states[span]
                da = 2*span*model.hx
                Gdisp = -(hi['energy_full_J']-lo['energy_full_J'])/(da*model.thickness)
                dC = (hi['compliance_m_N']-lo['compliance_m_N'])/da
                Gforce = base['P_N']**2*dC/(2*model.thickness)
                diffs.append({'half_span_elements': span, 'half_span_m': span*model.hx,
                              'G_fixed_displacement_J_m2': Gdisp,
                              'G_matched_force_J_m2': Gforce,
                              'fixed_force_P_N': base['P_N'],
                              'force_vs_displacement_relative': relative(Gforce, Gdisp)})
            row = {'n_half': n, 'elements_half': len(model.conn), 'dofs_half': len(base['u']),
                   'hx_m': model.hx, 'hy_m': model.hy,
                   'positive_element_Jacobian_m2': model.hx*model.hy/4,
                   'baseline': serial_state(base),
                   'crack_length_states': {str(k): serial_state(s) for k, s in states.items()},
                   'domain_outer_radii_m': cfg['domain_outer_radii_m'],
                   'domain_J_J_m2': J.tolist(), 'domain_J_mean_J_m2': float(J.mean()),
                   'domain_relative_spread': float((J.max()-J.min())/J.mean()),
                   'energy_derivatives': diffs,
                   'domain_J_vs_Gdisp_max_relative': float(np.max(np.abs(J-diffs[0]['G_fixed_displacement_J_m2']))/diffs[0]['G_fixed_displacement_J_m2'])}
            records.append(row)
            print(name, n, 'J=', J.tolist(), 'G=', diffs[0]['G_fixed_displacement_J_m2'], flush=True)
            del model, states, base
        coarse, fine = records[-2:]
        final_diffs = fine['energy_derivatives']
        checks = {
            'positive_J_and_G': all(r['domain_J_mean_J_m2'] > 0 and r['energy_derivatives'][0]['G_fixed_displacement_J_m2'] > 0 for r in records),
            'J_vs_energy': fine['domain_J_vs_Gdisp_max_relative'] <= gates['finest_J_vs_fixed_displacement_energy_relative'],
            'domain_dependence': fine['domain_relative_spread'] <= gates['finest_J_domain_spread_relative'],
            'J_refinement': relative(fine['domain_J_mean_J_m2'], coarse['domain_J_mean_J_m2']) <= gates['final_two_meshes_J_relative'],
            'G_refinement': relative(final_diffs[0]['G_fixed_displacement_J_m2'], coarse['energy_derivatives'][0]['G_fixed_displacement_J_m2']) <= gates['final_two_meshes_G_relative'],
            'energy_increment': relative(final_diffs[0]['G_fixed_displacement_J_m2'], final_diffs[1]['G_fixed_displacement_J_m2']) <= gates['finest_h_vs_2h_energy_derivative_relative'],
            'matched_controls': final_diffs[0]['force_vs_displacement_relative'] <= gates['finest_fixed_force_vs_fixed_displacement_G_relative'],
            'all_states_equilibrium_and_work': all(max(s['relative_free_residual'], s['relative_force_balance'], s['relative_work_error'], abs(s['horizontal_anchor_reaction_N']/s['P_N'])) <= gates['equilibrium_work_and_reaction_relative'] for r in records for s in r['crack_length_states'].values()),
        }
        material_results[name] = {'material': mat, 'meshes': records, 'checks': checks,
                                  'J_final_refinement_relative': relative(fine['domain_J_mean_J_m2'], coarse['domain_J_mean_J_m2']),
                                  'G_final_refinement_relative': relative(final_diffs[0]['G_fixed_displacement_J_m2'], coarse['energy_derivatives'][0]['G_fixed_displacement_J_m2']),
                                  'G_h_vs_2h_relative': relative(final_diffs[0]['G_fixed_displacement_J_m2'], final_diffs[1]['G_fixed_displacement_J_m2']),
                                  'passed': all(checks.values())}
    evidence = ROOT.parent/'bark-constitutive-upgrade'/'fracture-scale-evidence'
    files = [ROOT/'PROTOCOL.json', ROOT/'SOURCES.md', ROOT/'README.md', *sorted((ROOT/'source').glob('*.py')),
             *sorted((ROOT/'tests').glob('*.py')), evidence/'DECISION.md', evidence/'fracture_scale_receipt.json']
    elapsed = time.monotonic()-started
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    receipt = {'status': 'PASS_NUMERICAL_REFERENCE_ONLY' if all(r['passed'] for r in material_results.values()) else 'FAIL_DECLARED_NUMERICAL_GATES',
               'limitations': 'Independent 2D elastic Q4 reference only. Does not validate the production 3D prism/cohesive solver, identify bark elastic/fracture data, resolve a cohesive process zone, predict morphology, or determine a physical mesh pitch.',
               'protocol': cfg, 'materials': material_results,
               'primary_derivation_source': {'url': 'https://esag.harvard.edu/rice/015_Rice_PathIndepInt_JAM68.pdf',
                   'verified': 'Indexed primary first-page text, equation (2), orientation and assumptions, 2026-10-06',
                   'direct_PDF_access': 'HTTP 502; no PDF bytes or external PDF hash available'},
               'source_sha256': {str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p.relative_to(ROOT.parent)): sha(p) for p in files},
               'environment': {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__},
               'wall_seconds': elapsed, 'maximum_resident_memory_KiB': rss,
               'resource_gates_passed': elapsed < gates['process_wall_seconds'] and rss < gates['process_maximum_resident_KiB']}
    (ROOT/'receipts'/'reference_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(receipt['status'], 'seconds=', elapsed, 'peak_KiB=', rss, flush=True)
    if not all(r['passed'] for r in material_results.values()) or not receipt['resource_gates_passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
