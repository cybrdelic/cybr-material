"""An exact missed P2 mode and an independent nonlinear quadrature comparison."""
import json,time
import numpy as np
from periodic_sector import PeriodicSector,ROOT
from single_seam import SingleSeam,ConditionalFracture

start=time.monotonic();control=PeriodicSector();fixture=ConditionalFracture(1.,1.,1000.,'Quadrature-only unit fixture, not a physical cork material.',evidence_status='synthetic unit test only')

def reduced(seam,full):
    counts=np.asarray(seam.T.power(2).sum(axis=0)).ravel()
    return (np.asarray(seam.T.T@full.ravel())/counts).reshape(seam.X.shape)

def nodal_field(seam,nonlinear=False):
    c=seam.config;t=np.array(c['cohort_thickness_m']);inner=c['outer_radius_m']-t.sum();bounds=inner+np.r_[0.,t.cumsum()]
    full=np.zeros_like(seam.base.X);ids=seam.duplicated_nodes;X=seam.base.X[ids]
    layer=np.minimum(np.searchsorted(bounds[1:],X[:,2],side='right'),2)
    zeta=2*(X[:,2]-.5*(bounds[layer]+bounds[layer+1]))/t[layer]
    if not nonlinear:
        jump=np.zeros_like(X);jump[:,0]=1e-5*(zeta*zeta-1/3);common=np.zeros_like(X)
    else:
        wave=np.sin(2*np.pi*X[:,1]/c['height_m'])
        jump=np.column_stack((-1e-4*(1+.4*zeta*zeta),np.zeros(len(X)),2e-4*wave*(zeta*zeta+.3)))
        common=np.column_stack((1e-3*(wave+.2*zeta*zeta),np.zeros(len(X)),np.zeros(len(X))))
    full[ids]=common-.5*jump;full[len(control.base.X):]=common+.5*jump
    return reduced(seam,full)

null=[];nonlinear=[];gradients=[]
for order in (2,3,4,5,7,9):
    seam=SingleSeam(control,face_order=order,allow_underintegrated_control=(order==2));bonds=seam.bind(fixture)
    u=nodal_field(seam);energy,gradient,opening=seam.interface(u,bonds)
    exact=.5*fixture.penalty_Pa_per_m*float(seam.geometry['area'].sum())*1e-10*(4/45)
    null.append(dict(order=order,energy_J=energy,analytic_energy_J=exact,relative_error=abs(energy-exact)/exact,
                     maximum_sampled_opening_m=float(opening.max())))
    u=nodal_field(seam,nonlinear=True);energy,g,opening=seam.interface(u,bonds,np.full(bonds.N,.55))
    nonlinear.append(dict(order=order,energy_J=energy));gradients.append(g)
for row,g in zip(nonlinear,gradients):
    row['energy_relative_difference_from_order9']=abs(row['energy_J']-nonlinear[-1]['energy_J'])/nonlinear[-1]['energy_J']
    row['gradient_relative_difference_from_order9']=float(np.linalg.norm(g-gradients[-1])/np.linalg.norm(gradients[-1]))
assert null[0]['relative_error']>.999999
assert all(row['relative_error']<1e-11 for row in null[1:])
result=dict(status='PASS: two-point rule rejected; quadratic undamaged face energy resolved at order 3+',
            exact_negative_control=null,nonlinear_damaged_comparison=nonlinear,
            analytical_mean_gap_squared='mean[(zeta^2-1/3)^2] = 4/45; the 2-point nodes zeta=+-1/sqrt(3) see zero',
            field='Same P2 nodal geometry at every order. Nonlinear control has finite common face bending, compressive/mixed jump and fixed damage 0.55. No damage history or material solve.',
            limitation='Order 3 exactly integrates the planar undamaged quadratic-gap spring energy. Nonlinear normals, closure boundaries and evolving damage require separate quadrature convergence; same-rule gradient checks do not establish integration accuracy.',
            old_attempt='attempts/two_point retains the source, topology and initially passing gradient receipts now known to miss this mode.',
            selected_default_order=4,selected_nonlinear_gradient_relative_error=nonlinear[2]['gradient_relative_difference_from_order9'],
            selected_sample_gradient_gate=1e-4,wall_s=time.monotonic()-start)
assert result['selected_nonlinear_gradient_relative_error']<result['selected_sample_gradient_gate']
(ROOT/'receipts/interface_quadrature.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
