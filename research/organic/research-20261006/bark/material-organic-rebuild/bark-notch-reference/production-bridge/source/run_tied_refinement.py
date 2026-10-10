"""Additional 40-mesh sensitivity; preserves the original doubling-gate failure."""
import argparse
import json
import resource
import time
import numpy as np
from scipy.sparse import csr_matrix
from numpy.polynomial.legendre import leggauss
from bridge import Bridge, ROOT, serial
from run_mesh import check_frozen, normrel, rel, sha


class ExactTieBridge(Bridge):
    def __init__(self,n,cfg):
        super().__init__(n,cfg)
        # The exact production quadrature has already been packed into bulk
        # arrays and its tangent assembled. Drop redundant per-cell storage;
        # retain one cell only for the quadrature-count receipt.
        self.cells=self.cells[:1]
        self.bulk.cells=()

    def interface_tangent(self):
        # All trace DOFs will be identified. The exact constrained interface
        # Hessian is zero. Native interface evaluation remains active.
        return csr_matrix((self.N*12,self.N*12))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--protocol',default='ADDITIONAL_REFINEMENT.json')
    args=parser.parse_args()
    start=time.monotonic();check_frozen()
    cfg=json.loads((ROOT/'PROTOCOL.json').read_text())
    extra=json.loads((ROOT/args.protocol).read_text())
    ref=json.loads((ROOT.parent/'receipts/reference_receipt.json').read_text())['materials']['isotropic_nu0']['meshes'][-1]['baseline']
    targets={'P_N':ref['P_N']*.01,'compliance_m_N':ref['compliance_m_N'],'energy_full_J':ref['energy_full_J']*.0001}
    m=ExactTieBridge(extra['mesh_n'],cfg);m.prepare(None)
    ev=m.solve(cfg['half_grip_displacement_m'])
    rng=np.random.default_rng(482);v=rng.normal(size=m.size)*m.h*1e-7
    plus,minus=m.evaluate(v),m.evaluate(-v)
    tangent=normrel((plus['gb']-minus['gb'])/2,m.Kb@v)
    interface=float(np.linalg.norm(plus['gc'])/np.linalg.norm(m.Kb@v))
    z,w=leggauss(3);work=0.
    for lam,weight in zip((z+1)/2,w/2):
        state=m.solve(cfg['half_grip_displacement_m']*lam)
        work+=2*cfg['half_grip_displacement_m']*weight*state['P_N']
    previous=json.loads((ROOT/'receipts/mesh_32.json').read_text())['cases'][0]['production_equilibrium']
    changes={k:rel(ev[k],previous[k]) for k in targets}
    checks={'actual_bulk_tangent':tangent<cfg['gates']['native_force_tangent_central_difference_relative'],
            'actual_interface_force':interface<cfg['gates']['native_interface_force_vs_assembled_force_relative'],
            'actual_equilibrium':ev['free_residual_relative']<cfg['gates']['nonlinear_free_residual_relative'],
            'actual_work':rel(work,ev['energy_full_J'])<cfg['gates']['integrated_grip_work_vs_actual_energy_relative'],
            'reference_match':max(rel(ev[k],v) for k,v in targets.items())<cfg['gates']['finest_tied_reaction_compliance_energy_vs_frozen_reference_relative']}
    check_frozen()
    result={'status':'ADDITIONAL_SENSITIVITY_LOCAL_CHECKS_PASS' if all(checks.values()) else 'ADDITIONAL_SENSITIVITY_FAIL',
            'n':extra['mesh_n'],'prisms':m.N,'dofs':m.size,'production_equilibrium':serial(ev),
            'relative_difference_from_reference':{k:rel(ev[k],v) for k,v in targets.items()},
            'relative_change_from_32':changes,'original_16_to_32_doubling_gate_remains_failed':True,
            'bulk_tangent_relative':tangent,'interface_force_relative':interface,
            'external_work_full_J':work,'work_error_relative':rel(work,ev['energy_full_J']),
            'checks':checks,'frozen_reference_and_production_sources_unchanged':True,
            'wall_seconds':time.monotonic()-start,
            'maximum_resident_memory_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'source_sha256':{p.name:sha(p) for p in (ROOT/'source').glob('*.py')},
            'extra_protocol_sha256':sha(ROOT/args.protocol)}
    result['resource_gates_passed']=result['wall_seconds']<115 and result['maximum_resident_memory_KiB']<1048576
    (ROOT/'receipts'/f"mesh_{extra['mesh_n']}_tied.json").write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('production_equilibrium','source_sha256')},indent=2),flush=True)
    if not all(checks.values()) or not result['resource_gates_passed']:raise SystemExit(1)


if __name__=='__main__':main()
