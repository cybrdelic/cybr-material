from pathlib import Path
import numpy as np,json,hashlib
from scipy.spatial.transform import Rotation
from quadratic_prism import QuadraticPrism,flat_reference,shapes
R=Path(__file__).resolve().parents[1];xy=np.array([[0.,0.],[.002,0.],[.0003,.0015]]);t=.0002;X=flat_reference(xy,t);p=QuadraticPrism(X,2.8e6,.28,240);checks=[];rng=np.random.default_rng(721)
errors=[]
for _ in range(30):
 r,s=rng.dirichlet([1,1,1])[:2];z=rng.uniform(-1,1);N,dN=shapes(r,s,z);errors.append(max(abs(N.sum()-1),abs(dN.sum(0)).max()))
checks.append({'case':'shape_partition_and_derivative','maximum_error':max(errors),'pass':max(errors)<1e-14})
volume=abs(np.linalg.det(np.column_stack((xy[1]-xy[0],xy[2]-xy[0]))))/2*t;checks.append({'case':'reference_volume_and_mass','volume_relative_error':abs(p.reference_volume-volume)/volume,'pass':abs(p.reference_volume-volume)/volume<1e-13 and abs(p.mass_scalar.sum()-p.mass)/p.mass<1e-13})
Q=Rotation.from_rotvec([1.5,-.7,.6]).as_matrix();x=X@Q.T+[.003,-.004,.002];E,g,v=p.evaluate(x);checks.append({'case':'arbitrary_finite_rigid_motion','energy_J':E,'gradient_norm_N':float(np.linalg.norm(g)),'pass':E<1e-24 and np.linalg.norm(g)<1e-12})
F=np.array([[1.12,.03,0],[0,.96,.02],[0,0,.91]]);x=X@F.T;E,g,v=p.evaluate(x);J=np.linalg.det(F);psi=p.mu/2*(np.sum(F*F)-3)-p.mu*np.log(J)+p.lam/2*np.log(J)**2;expected=volume*psi;checks.append({'case':'homogeneous_finite_strain_and_volume','energy_relative_error':abs(E-expected)/expected,'volume_relative_error':abs(v['mapped_volume_m3']-volume*J)/(volume*J),'pass':abs(E-expected)/expected<1e-12 and abs(v['mapped_volume_m3']-volume*J)/(volume*J)<1e-13})
x+=rng.normal(size=x.shape)*1e-6;E,g,_=p.evaluate(x);d=rng.normal(size=x.shape);d/=np.linalg.norm(d);eps=1e-9;fd=(p.evaluate(x+eps*d)[0]-p.evaluate(x-eps*d)[0])/(2*eps);exact=float((g*d).sum());err=abs(fd-exact)/abs(exact);checks.append({'case':'analytic_force_gradient','relative_error':err,'pass':err<1e-7})
Er,gr,_=p.evaluate(x@Q.T+[.003,-.004,.002]);checks.append({'case':'deformed_energy_and_force_objectivity','energy_relative_error':abs(Er-E)/E,'force_covariance_error':float(np.linalg.norm(gr-g@Q.T)/np.linalg.norm(g)),'pass':abs(Er-E)/E<1e-12 and np.linalg.norm(gr-g@Q.T)/np.linalg.norm(g)<1e-12})
# Exact P2 pure bending displacement has zero linear shear; verify thin-limit
# bending energy at nu=0 rather than accepting a shear-locked element.
p0=QuadraticPrism(X,2.8e6,0,240);kappa=.01;x=X.copy();x[:,0]-=kappa*X[:,0]*X[:,2];x[:,2]+=.5*kappa*X[:,0]**2;Eb,gb,_=p0.evaluate(x);expected=.5*2.8e6*kappa*kappa*volume*t*t/12;err=abs(Eb-expected)/expected;checks.append({'case':'pure_bending_small_strain_limit_no_shear_lock','relative_error':err,'pass':err<1e-5})
# Material may contract through thickness while preserving birth mass.
x=X.copy();x[:,2]*=.8;_,_,v=p.evaluate(x);checks.append({'case':'actual_thickness_contraction_and_mass','density_ratio':v['mapped_mean_density_kg_m3']/240,'pass':abs(v['mapped_mean_density_kg_m3']/240-1.25)<1e-13 and v['material_mass_kg']==p.mass})
x=X.copy();x[:,2]*=-1
try:p.evaluate(x);rejected=False
except ValueError:rejected=True
checks.append({'case':'inversion_rejected','pass':rejected})
out={'all_pass':all(c['pass'] for c in checks),'checks':checks,'source_sha256':hashlib.sha256((R/'source/quadratic_prism.py').read_bytes()).hexdigest(),'scope':'Finite-volume material-cell kernel. The isotropic neo-Hookean benchmark is not a calibrated cork constitutive law; assembly, cohesive/contact, accretion transport and morphology remain unqualified.'};(R/'receipts/prism_tests.json').write_text(json.dumps(out,indent=2,default=lambda v:v.item()));print(json.dumps(out,indent=2,default=lambda v:v.item()));assert out['all_pass']
