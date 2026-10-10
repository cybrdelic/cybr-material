"""Report independent numerical changes without treating damage as crack length."""
import argparse,json
import numpy as np
from opening_model import ROOT
from checkpoint_io import read_checkpoint

def compare(left,right,kind):
    aa,a=read_checkpoint(ROOT/'state'/left/'restart.npz');bb,b=read_checkpoint(ROOT/'state'/right/'restart.npz')
    assert a['status']==b['status']=='completed'
    ta=a['trace'];tb=b['trace'];fa=ta[-1];fb=tb[-1]
    la=np.array([r['load'] for r in ta]);lb=np.array([r['load'] for r in tb])
    grid=np.unique(np.r_[la,lb]);Qa=np.interp(grid,la,[r['eigenstrain_conjugate_J'] for r in ta]);Qb=np.interp(grid,lb,[r['eigenstrain_conjugate_J'] for r in tb])
    def outer_gap(data):
        X=data['reference_nodes'];u=data['full_displacement'];ids=np.flatnonzero((abs(X[:,0])<1e-12)&(abs(X[:,2]-.105)<1e-12));values=[]
        for y in np.unique(X[ids,1]):
            pair=ids[abs(X[ids,1]-y)<1e-12];assert len(pair)==2
            values.append(float(np.linalg.norm(u[pair[1]]-u[pair[0]])))
        return max(values)
    ga=outer_gap(aa);gb=outer_gap(bb);gape=abs(ga-gb)/gb
    we=abs(fa['work_J']-fb['work_J'])/fb['work_J']
    pe=float(abs(max(Qa)-max(Qb))/max(Qb))
    result=dict(left=left,right=right,kind=kind,load=fa['load'],same_physical_scenario=kind!='physical_scenario',
        maximum_gap_m=[ga,gb],maximum_gap_relative_change=gape,
        gap_metric='Maximum solved outer-edge gap at common material nodes; avoids comparing different quadrature sampling positions',
        maximum_quadrature_sample_gap_m=[fa['maximum_opening_m'],fb['maximum_opening_m']],
        total_work_J=[fa['work_J'],fb['work_J']],total_work_relative_change=we,maximum_work_conjugate_relative_change=pe,
        work_conjugate_curve_relative_RMS=float(np.linalg.norm(Qa-Qb)/np.linalg.norm(Qb)),
        new_dissipation_J=[fa['new_dissipation_from_initial_J'],fb['new_dissipation_from_initial_J']],
        stored_energy_J=[fa['stored_energy_J'],fb['stored_energy_J']],
        work_ledger_relative_errors=[abs(fa['energy_work_error_J'])/fa['work_J'],abs(fb['energy_work_error_J'])/fb['work_J']],
        minimum_physical_J=[fa['minimum_physical_J'],fb['minimum_physical_J']],
        maximum_combined_free_force_N=[fa['maximum_free_force_N'],fb['maximum_free_force_N']],
        interpretation='Numerical changes in the same conditional opening model. First nonzero integration-point damage is not a converged crack-initiation load; fully separated extension and active zone are inspected separately.')
    if kind!='physical_scenario':
        result['engineering_comparison_gates']=dict(gap_under_two_percent=gape<.02,work_under_one_percent=we<.01,peak_conjugate_under_two_percent=pe<.02,
            both_work_ledgers_under_one_per_mille=max(result['work_ledger_relative_errors'])<.001,
            both_combined_forces_under_five_e_minus5_N=max(result['maximum_combined_free_force_N'])<5e-5,
            both_physical_orientations_positive=min(result['minimum_physical_J'])>0)
        result['gate_scope']='Declared numerical engineering tolerances for a conditional mechanism control, not biological validation or appearance acceptance.'
    (ROOT/'receipts'/f'compare_{kind}.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('left');p.add_argument('right');p.add_argument('kind');a=p.parse_args();compare(a.left,a.right,a.kind)
