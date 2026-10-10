from pathlib import Path
import argparse
import hashlib
import json
import resource
import time
import numpy as np
from numpy.polynomial.legendre import leggauss
from bridge import Bridge, ROOT, serial


def rel(a,b):return float(abs(a-b)/abs(b))
def normrel(a,b):return float(np.linalg.norm(a-b)/np.linalg.norm(b))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def check_frozen():
    ref=ROOT.parent
    snapshot=json.loads((ROOT/'receipts/frozen_reference_snapshot.json').read_text())
    assert sha(ref/'receipts/manifest.json')==snapshot['manifest_sha256']
    for path,h in snapshot['sha256'].items():assert sha(ref/path)==h,path
    source=json.loads((ROOT/'receipts/production_source_snapshot.json').read_text())
    for path,h in source['sha256'].items():
        assert sha(Path(source['origin'])/path)==h,path
        assert sha(ROOT/'frozen-production'/path)==h,path


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--n',type=int,required=True)
    args=parser.parse_args();start=time.monotonic();check_frozen()
    cfg=json.loads((ROOT/'PROTOCOL.json').read_text());gate=cfg['gates']
    ref=json.loads((ROOT.parent/'receipts/reference_receipt.json').read_text())['materials']['isotropic_nu0']['meshes'][-1]['baseline']
    scale=cfg['frozen_reference_load_scale']
    targets={'P_N':ref['P_N']*scale,'compliance_m_N':ref['compliance_m_N'],
             'energy_full_J':ref['energy_full_J']*scale**2}
    m=Bridge(args.n,cfg);rows=[]
    for eta in [None]+cfg['finite_penalty_eta']:
        case_start=time.monotonic();m.prepare(eta);ev=m.solve(cfg['half_grip_displacement_m'])
        linear=m.evaluate(ev['linear_u']);assembled=m.K@ev['linear_u']
        nonlinear_vs_linear=normrel(linear['g'],assembled)
        # Central tangent check uses actual native forces, separately for bulk
        # and interfaces, so stiff penalties cannot hide a wrong bulk tangent.
        rng=np.random.default_rng(482);v=rng.normal(size=m.size)*m.h*1e-7
        plus=m.evaluate(v);minus=m.evaluate(-v)
        bulk_tangent_error=normrel((plus['gb']-minus['gb'])/2,m.Kb@v)
        if eta is None:
            interface_error=float(np.linalg.norm(plus['gc'])/np.linalg.norm(m.Kb@v))
        else:interface_error=normrel(plus['gc'],m.Kc@v)
        native_energy_tangent=rel((plus['energy_half_J']+minus['energy_half_J'])/2,.5*v@(m.K@v))
        eps=1e-4
        Ep=m.evaluate(ev['u']*(1+eps))['energy_half_J'];Em=m.evaluate(ev['u']*(1-eps))['energy_half_J']
        derivative=(Ep-Em)/(2*eps);analytic=float(ev['g']@ev['u'])
        energy_gradient_error=rel(derivative,analytic)
        # Exact conservative work over the equilibrated loading branch. The
        # nonlinear production law must not be judged solely by 1/2 P Delta.
        z,w=leggauss(3);work=0.;ramp=[]
        for lam,weight in zip((z+1)/2,w/2):
            state=m.solve(cfg['half_grip_displacement_m']*lam)
            work+=2*cfg['half_grip_displacement_m']*weight*state['P_N']
            ramp.append({'load_fraction':float(lam),'P_N':state['P_N'],
                         'residual':state['free_residual_relative']})
        checks={
            'actual_bulk_tangent':bulk_tangent_error<=gate['native_force_tangent_central_difference_relative'],
            'actual_interface_force':interface_error<=gate['native_interface_force_vs_assembled_force_relative'],
            'actual_energy_tangent':native_energy_tangent<=gate['native_force_tangent_central_difference_relative'],
            'small_load_tangent_limit':nonlinear_vs_linear<=gate['native_nonlinear_force_vs_assembled_tangent_at_small_load_relative'],
            'actual_energy_gradient':energy_gradient_error<=gate['native_total_energy_directional_derivative_relative'],
            'actual_equilibrium':ev['free_residual_relative']<=gate['nonlinear_free_residual_relative'],
            'reaction_balance':ev['reaction_balance_relative']<=gate['reaction_balance_relative'],
            'integrated_work':rel(work,ev['energy_full_J'])<=gate['integrated_grip_work_vs_actual_energy_relative'],
            'undamaged':ev['max_opening_over_initiation']<gate['maximum_opening_over_initiation'],
            'positive_orientation':ev['min_det_F']>0,
        }
        row={'eta':eta,'case':'exact_tie' if eta is None else 'finite_penalty',
             'dofs':m.size,'production_equilibrium':serial(ev),
             'targets_scaled_frozen_reference':targets,
             'relative_difference_from_reference':{k:rel(ev[k],v) for k,v in targets.items()},
             'force_checks':{'bulk_central_tangent_relative':bulk_tangent_error,
                            'interface_assembled_relative':interface_error,
                            'total_energy_central_tangent_relative':native_energy_tangent,
                            'nonlinear_vs_tangent_at_linear_solution_relative':nonlinear_vs_linear,
                            'actual_energy_directional_derivative_relative':energy_gradient_error},
             'external_work_full_J':work,'work_error_relative':rel(work,ev['energy_full_J']),
             'ramp_samples':ramp,'checks':checks,'passed':all(checks.values()),
             'case_wall_seconds':time.monotonic()-case_start}
        rows.append(row)
        print(json.dumps({'n':args.n,'eta':eta,'P_N':ev['P_N'],'reference_errors':row['relative_difference_from_reference'],
                          'max_opening_over_initiation':ev['max_opening_over_initiation'],
                          'checks':checks,'seconds':row['case_wall_seconds']}),flush=True)
    check_frozen()
    result={'status':'PASS_LOCAL_OPERATOR_CHECKS' if all(x['passed'] for x in rows) else 'FAIL_LOCAL_OPERATOR_CHECKS',
            'n':args.n,'prisms':m.N,'quadrature_points_per_prism':len(m.cells[0].weights),
            'cohesive_quadrature_points':len(m.ops_unit['area']),
            'protocol_sha256':sha(ROOT/'PROTOCOL.json'),
            'bridge_source_sha256':{p.name:sha(p) for p in sorted((ROOT/'source').glob('*.py'))},
            'frozen_reference_and_original_sources_unchanged':True,'cases':rows,
            'wall_seconds':time.monotonic()-start,
            'maximum_resident_memory_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    result['resource_gates_passed']=(result['wall_seconds']<gate['wall_seconds_each_process'] and
                                   result['maximum_resident_memory_KiB']<gate['maximum_resident_memory_KiB'])
    (ROOT/'receipts'/f'mesh_{args.n}.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['wall_seconds'],result['maximum_resident_memory_KiB'],flush=True)
    if result['status']!='PASS_LOCAL_OPERATOR_CHECKS' or not result['resource_gates_passed']:raise SystemExit(1)


if __name__=='__main__':main()
