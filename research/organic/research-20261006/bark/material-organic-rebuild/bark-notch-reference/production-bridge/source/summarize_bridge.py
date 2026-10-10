"""Audit both limits, retain the failed first ladder and bound the conclusion."""
import json
from pathlib import Path
from run_mesh import ROOT, check_frozen, rel, sha


def main():
    check_frozen()
    cfg=json.loads((ROOT/'PROTOCOL.json').read_text());gates=cfg['gates']
    meshes=[json.loads((ROOT/'receipts'/f'mesh_{n}.json').read_text()) for n in cfg['meshes_n_square_half']]
    extra40=json.loads((ROOT/'receipts/mesh_40_tied.json').read_text())
    fine=json.loads((ROOT/'receipts/mesh_64_tied.json').read_text())
    for record in meshes+[extra40,fine]:
        hashes=record.get('bridge_source_sha256',record.get('source_sha256'))
        for name,h in hashes.items():assert sha(ROOT/'source'/name)==h, name
    quantities=['P_N','compliance_m_N','energy_full_J']
    tied16,tied32=[r['cases'][0]['production_equilibrium'] for r in meshes[-2:]]
    original_change={k:rel(tied32[k],tied16[k]) for k in quantities}
    lastpen,prevpen=meshes[-1]['cases'][-1]['production_equilibrium'],meshes[-1]['cases'][-2]['production_equilibrium']
    penalty_error={k:rel(lastpen[k],tied32[k]) for k in quantities}
    penalty_change={k:rel(lastpen[k],prevpen[k]) for k in quantities}
    ref=json.loads((ROOT.parent/'receipts/reference_receipt.json').read_text())['materials']['isotropic_nu0']['meshes']
    ref_change={k:rel(ref[-1]['baseline'][k],ref[-2]['baseline'][k]) for k in quantities}
    checks={
        'all_local_native_operator_checks':all(c['passed'] for r in meshes for c in r['cases']) and all(fine['checks'].values()) and all(extra40['checks'].values()),
        'finest_tied_matches_frozen_reference':max(fine['relative_difference_from_reference'].values())<=gates['finest_tied_reaction_compliance_energy_vs_frozen_reference_relative'],
        '32_to_64_tied_mesh_limit':max(fine['relative_change_from_32'].values())<=gates['final_two_tied_meshes_relative'],
        '32_smallest_penalty_vs_exact_tie':max(penalty_error.values())<=gates['finest_smallest_penalty_vs_exact_tie_relative'],
        '32_final_two_penalties':max(penalty_change.values())<=gates['finest_final_two_penalties_relative'],
        'all_completed_process_resource_gates':all(r['resource_gates_passed'] for r in meshes+[extra40,fine]),
        'frozen_reference_and_original_production_sources_unchanged':True,
    }
    summary={'status':'PASS_CONSTRAINED_SMALL_LOAD_PRODUCTION_BRIDGE' if all(checks.values()) else 'FAIL_DECLARED_BRIDGE_GATES',
             'checks':checks,'original_16_to_32_mesh_change':original_change,
             'original_16_to_32_gate_passed':max(original_change.values())<=gates['final_two_tied_meshes_relative'],
             'final_32_to_64_tied_mesh_change':fine['relative_change_from_32'],
             'final_tied_vs_frozen_reference':fine['relative_difference_from_reference'],
             '32_smallest_penalty_vs_exact_tie':penalty_error,
             '32_final_two_penalty_change':penalty_change,
             'frozen_reference_64_to_128_mesh_change':ref_change,
             'interpretation':'Cross-solver agreement compares two finite discretizations, not a rigorous continuum error. The frozen Q4 reaction/compliance/energy changes about 0.35% at its last refinement. Finite penalty softening can cancel bulk discretization error and is kept separate.',
             'first_64_attempt':'Stopped during a redundant 30.1MiB operator-array copy by the unchanged 1GiB virtual-memory limit. Retained log: receipts/mesh_64_attempt1_memory.log. Removed duplicate storage only; final run completed within the same cap.',
             'maximum_completed_process_wall_seconds':max(r['wall_seconds'] for r in meshes+[extra40,fine]),
             'maximum_completed_process_resident_KiB':max(r['maximum_resident_memory_KiB'] for r in meshes+[extra40,fine]),
             'limitations':['Only isotropic nu=0, constrained in-plane, small-load, undamaged response is tested.',
                            'No production J/domain integral or crack-advance energy-release calculation is tested in this bridge.',
                            'No nonlinear damage, dissipation, mixed-mode failure, arbitrary finite deformation, general 3D behavior, growth/birth mechanics or morphology validation.',
                            'No bark compliance, toughness, cohesive strength, crack-system calibration or physical mesh recommendation.',
                            'The 64 mesh is exact-tie only; finite-penalty refinement reaches 32 and is compared against its own tied control.'],
             'reproduction':'Run bash run_mesh.sh for n=4,8,16,32, then bash run_tied_refinement.sh and bash run_tied_64.sh, then python source/summarize_bridge.py under the same one-thread/resource environment.',
             'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in sorted((ROOT/'source').glob('*.py'))}}
    (ROOT/'receipts/bridge_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    delivered=[p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='manifest.json']
    manifest={'status':summary['status'],'originals_unchanged':True,
              'sha256':{str(p.relative_to(ROOT)):sha(p) for p in sorted(delivered)}}
    (ROOT/'receipts/manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    if not all(checks.values()):raise SystemExit(1)


if __name__=='__main__':main()
