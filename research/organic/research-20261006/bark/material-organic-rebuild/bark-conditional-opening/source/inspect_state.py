"""Actual gaps, tractions, dissipated fraction and resolved-zone diagnostics."""
import argparse,csv,json
import numpy as np
from opening_model import OpeningModel,History,ROOT
from checkpoint_io import read_checkpoint

def inspect(label):
    a,meta=read_checkpoint(ROOT/'state'/label/'restart.npz');p=meta['parameters']
    m=OpeningModel(p['radial_pitch'],p['quadrature'],p['penalty'],p['Gc']);u=a['displacement'];h=History.make(a['maximum'],a['damage']);m.set_load(meta['trace'][-1]['load'])
    _,g,_,r,opening=m.evaluate(u,h,work=True)
    o=m.seam.geometry;full=m.seam.expand(u);local=full[m.seam.cell_map]
    pair=np.concatenate((local[o['pairs'][:,0]],local[o['pairs'][:,1]]),axis=1);pair-=pair[:,:1]
    jump=np.einsum('qi,qia->qa',o['J'],pair)
    t1=o['tref'][:,:3]+np.einsum('qi,qia->qa',o['T1'],pair);t2=o['tref'][:,3:]+np.einsum('qi,qia->qa',o['T2'],pair)
    cross=np.cross(t1,t2);normal=o['normal_sign'][:,None]*cross/np.linalg.norm(cross,axis=1)[:,None]
    dn=np.einsum('qa,qa->q',jump,normal);jt=jump-dn[:,None]*normal
    K=m.contract.penalty_Pa_per_m;D=h.damage
    normal_traction=K*((1-D)*dn+D*np.minimum(dn,0.));shear_traction=(1-D)*K*np.linalg.norm(jt,axis=1)
    equivalent=np.sqrt(np.sum(jt*jt,axis=1)+np.maximum(dn,0.)**2)
    assert np.max(abs(equivalent-opening))<1e-12
    fraction=m.phi(h.maximum)/m.contract.Gc_J_m2
    new=~m.starter
    separated=new&(h.maximum>=m.contract.delta_f*(1-1e-10))&(D>=1-1e-10)&(dn>1e-8)&(np.maximum(normal_traction,0.)<m.contract.peak_Pa*1e-4)&(shear_traction<m.contract.peak_Pa*1e-4)
    active=new&(fraction>1e-8)&(fraction<1-1e-8)&(normal_traction>0.)
    radial=m.interface_reference[:,2];inner=m.base.config['outer_radius_m']-sum(m.base.config['cohort_thickness_m']);pitch=p['radial_pitch']
    slab=np.floor((radial-inner)/pitch+1e-7).astype(int)
    if np.any(active):
        span=float(radial[active].max()-radial[active].min());count=len(np.unique(slab[active]));bounds=[float(radial[active].min()),float(radial[active].max())]
    else:span=0.;count=0;bounds=[]
    work=meta['trace'][-1]['work_J'];balance=r['stored_energy_J']+r['new_dissipation_from_initial_J']-work
    rows=[]
    for k in np.unique(slab):
        mask=slab==k;w=m.area[mask];mean=lambda v:float(w@v[mask]/w.sum())
        rows.append(dict(radial_slab=int(k),inner_m=float(inner+k*pitch),outer_m=float(inner+(k+1)*pitch),
            mean_normal_opening_m=mean(dn),maximum_normal_opening_m=float(dn[mask].max()),mean_normal_traction_Pa=mean(normal_traction),
            minimum_normal_traction_Pa=float(normal_traction[mask].min()),maximum_normal_traction_Pa=float(normal_traction[mask].max()),
            mean_stiffness_damage=mean(D),mean_spent_Gc_fraction=mean(fraction),prepared=bool(np.all(m.starter[mask])),
            active_fraction=float(w@active[mask]/w.sum()),open_fully_separated_fraction=float(w@separated[mask]/w.sum())))
    result=dict(label=label,status=meta['status'],load=m.load,reference_mass_kg=m.seam.mass,maximum_equivalent_opening_m=float(opening.max()),
        maximum_normal_opening_m=float(dn.max()),maximum_stiffness_damage_new_ligament=float(D[new].max()),
        maximum_spent_Gc_fraction_new_ligament=float(fraction[new].max()),
        fully_separated_new_area_m2=float(m.area[separated].sum()),open_fully_separated_extension_equivalent_m=float(m.area[separated].sum()/m.base.config['height_m']),
        prepared_slit_m=.001,active_zone_radial_bounds_m=bounds,active_zone_point_span_m=span,active_radial_intervals=count,
        active_zone_point_span_over_pitch=span/pitch,active_zone_minimum_six_intervals_met=span/pitch>=6,
        stored_energy_J=r['stored_energy_J'],new_dissipation_J=r['new_dissipation_from_initial_J'],work_J=work,
        energy_work_error_J=balance,relative_work_error=abs(balance)/max(abs(work),1e-20),maximum_combined_free_force_N=r['maximum_free_force_N'],
        minimum_physical_J=r['minimum_physical_J'],maximum_tangential_strain_by_cohort=r['maximum_tangential_strain_by_cohort'],
        profile=rows,interpretation='D is stiffness loss, not crack advance. Fully separated new extension additionally requires saturation, a positive solved normal gap, and negligible cohesive traction. The active cohesive zone is reported separately. All geometry must use actual nodal displacement, never a D threshold.')
    with (ROOT/'receipts'/f'{label}_radial_profile.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    (ROOT/'receipts'/f'{label}_inspection.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k!='profile'},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('label');inspect(p.parse_args().label)
