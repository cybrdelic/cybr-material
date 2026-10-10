"""Common-load convergence and direct state/strain boundary diagnosis."""
import hashlib,json,resource,time
from pathlib import Path
import numpy as np
from model import Coupon,ROOT,DomainError,shapes


def surface_sample(model,u,n=49):
    xs=np.linspace(model.x[0],model.x[-1],n);ys=np.linspace(model.y[0],model.y[-1],n)
    reference=[];displacement=[];coordinates=[]
    nx=len(model.x)-1;ny=len(model.y)-1
    for y in ys:
        for x in xs:
            i=min(nx-1,int((x-model.x[0])/(model.x[-1]-model.x[0])*nx))
            j=min(ny-1,int((y-model.y[0])/(model.y[-1]-model.y[0])*ny))
            candidates=(2*(j*nx+i),2*(j*nx+i)+1)
            for triangle in candidates:
                p=model.xy[triangle]
                r,s=np.linalg.solve(np.column_stack((p[1]-p[0],p[2]-p[0])),np.array([x,y])-p[0])
                if r>=-1e-10 and s>=-1e-10 and r+s<=1+1e-10:break
            else:raise RuntimeError('Surface mapping failed')
            cell=2*model.T+triangle;N,_=shapes(r,s,1.)
            reference.append(N@model.Xcell[cell]);displacement.append(N@u[model.map[cell]]);coordinates.append([x,y])
    X=np.array(reference);U=np.array(displacement);xy=np.array(coordinates)
    radial=X.copy();radial[:,1]=0;radial/=np.linalg.norm(radial,axis=1)[:,None]
    dr=np.einsum('ia,ia->i',U,radial)
    return X,U,xy,dr


def strain_location(model,u):
    local=u[model.map];H=np.einsum('cia,cqib->cqab',local-local.mean(axis=1,keepdims=True),model.grad)
    F=H+np.eye(3);A=np.eye(3)+model.load*model.Aprime
    B=np.einsum('cqab,cqbd,cqde->cqae',F,A,model.axes)
    e=np.linalg.norm(B,axis=-2)-1
    qX=np.einsum('cqi,cia->cqa',np.array([c.N for c in model.cells]),model.Xcell)
    reports=[]
    for axis,name in enumerate(('radial','tangential','axial')):
        cell,point=np.unravel_index(np.argmax(e[...,axis]),e.shape[:2]);X=qX[cell,point]
        xy=np.array([np.arctan2(X[0],X[2])*model.config['outer_radius_m'],X[1]])
        reports.append(dict(direction=name,strain=float(e[cell,point,axis]),cohort=int(cell//model.T),cell=int(cell),quadrature_point=int(point),
                            reference_position_m=X.tolist(),surface_coordinates_m=xy.tolist(),
                            distance_to_fixed_substrate_m=float(np.linalg.norm(X[[0,2]])-(model.config['outer_radius_m']-sum(model.config['cohort_thickness_m']))),
                            distance_to_nearest_lateral_cut_m=float(min(model.config['width_m']/2-abs(xy[0]),model.config['height_m']/2-abs(xy[1])))))
    return reports


def main():
    started=time.monotonic(); rows=[]; samples={}
    labels=['coarse16_half','medium16_half','fine16_half','medium32_half']
    for label in labels:
        receipt=json.loads((ROOT/'state'/label/'restart.json').read_text());assert receipt['status']=='completed'
        src=ROOT/'state'/label/'restart.npz';assert hashlib.sha256(src.read_bytes()).hexdigest()==receipt['state_sha256']
        a=np.load(src);m=Coupon(receipt['parameters']['pitch_m'],precondition=False);m.set_load(.5);u=a['displacement'];X,U,xy,dr=surface_sample(m,u)
        interior=(abs(xy[:,0])<.012)&(abs(xy[:,1])<.012)
        row=dict(label=label,pitch_m=m.pitch,load_steps=receipt['parameters']['steps'],state_sha256=receipt['state_sha256'],
                 final=receipt['trace'][-1],maximum_surface_displacement_m=float(np.linalg.norm(U,axis=1).max()),
                 surface_radial_displacement_range_m=[float(dr.min()),float(dr.max())],
                 central_surface_radial_range_m=[float(dr[interior].min()),float(dr[interior].max())],
                 maximum_strain_locations=strain_location(m,u),peak_rss_MiB=receipt['peak_rss_MiB'])
        rows.append(row);samples[label]=U
        np.savez_compressed(ROOT/'state'/label/'surface_samples.npz',reference=X,displacement=U,material_xy=xy,radial_displacement=dr)
    comparisons=[]
    for left,right in [('coarse16_half','medium16_half'),('medium16_half','fine16_half'),('medium16_half','medium32_half')]:
        a=next(r for r in rows if r['label']==left);b=next(r for r in rows if r['label']==right)
        ua=samples[left];ub=samples[right];difference=ua-ub
        comparisons.append(dict(left=left,right=right,
            energy_relative_change=abs(a['final']['stored_energy_J']-b['final']['stored_energy_J'])/b['final']['stored_energy_J'],
            eigenstrain_work_relative_change=abs(a['final']['eigenstrain_work_J']-b['final']['eigenstrain_work_J'])/b['final']['eigenstrain_work_J'],
            surface_displacement_rms_difference_m=float(np.sqrt(np.mean(np.sum(difference**2,axis=1)))),
            surface_displacement_max_difference_m=float(np.linalg.norm(difference,axis=1).max()),
            surface_displacement_relative_rms_difference=float(np.linalg.norm(difference)/np.linalg.norm(ub))))
    # Reconstruct the preserved guard failure, capturing only a rejected trial.
    m=Coupon(.006);u=np.load(ROOT/'state/coarse16/restart.npz')['displacement'];m.set_load(.875)
    base=m.evaluate;invalid=[]
    def diagnostic_evaluate(candidate,*args,**kwargs):
        try:return base(candidate,*args,**kwargs)
        except DomainError as exc:
            invalid[:]=[(candidate.copy(),str(exc))];raise
    m.evaluate=diagnostic_evaluate
    try:m.equilibrate(u)
    except RuntimeError as exc:failure=str(exc)
    else:raise AssertionError('Expected preserved guard failure')
    failure_info=dict(load=.875,message=failure,rejected_trial_only=True,maximum_strain_locations=strain_location(m,invalid[-1][0]),
        interpretation='A declared 6% extrapolation scope in an unmeasured direction, not a measured material failure. Hard fixation meeting a free cut face causes a concentrated boundary response. It is a deliberate rigid-substrate coupon idealization, not a resolved physiological attachment law.')
    gates=dict(common_load_all_equilibria_converged=True,
               final_mesh_energy_below_one_percent=comparisons[1]['energy_relative_change']<.01,
               final_mesh_surface_rms_below_one_percent=comparisons[1]['surface_displacement_relative_rms_difference']<.01,
               final_increment_work_below_one_per_mille=comparisons[2]['eigenstrain_work_relative_change']<.001)
    result=dict(rows=rows,comparisons=comparisons,gates=gates,guard_failure=failure_info,
                morphology_claim='No fissures or biological bark formation are predicted: the qualified deformation is a continuous prefracture layer response. No material scene or color promotion is authorized by these numerical results.',
                rate_scope='The law is rate independent and these are quasistatic load-increment comparisons, not a measured strain-rate comparison.',
                wall_s=time.monotonic()-started,peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
    (ROOT/'receipts/comparison_and_limit.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(dict(comparisons=comparisons,gates=gates,guard_failure=failure_info),indent=2))

if __name__=='__main__':main()
