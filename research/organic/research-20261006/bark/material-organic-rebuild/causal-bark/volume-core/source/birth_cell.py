"""Retained angular material cohort -> curved finite-volume reference cell.
Birth geometry and material mass persist; subsequent prescribed annular growth
changes current nodes. Equilibrium/fracture/contact assembly remains separate.
"""
import numpy as np
from quadratic_prism import QuadraticPrism

def cylindrical_nodes(theta_y_triangle,mid_radius,thickness):
 a=np.asarray(theta_y_triangle,dtype=float);six=np.vstack((a,(a[0]+a[1])/2,(a[1]+a[2])/2,(a[2]+a[0])/2));out=[]
 # theta increases from +z toward +x, and y is axial. The r,s,zeta
 # orientation for (theta,y) then radial is right-handed.
 for z in (-thickness/2,0,thickness/2):
  radius=mid_radius+z;out.extend(np.column_stack((radius*np.sin(six[:,0]),six[:,1],radius*np.cos(six[:,0]))))
 return np.asarray(out)

def make_cohort_cell(history,cohort_id,theta_y_triangle,E,nu,nq=5):
 if history.radius.shape!=(1,1):raise ValueError('Local uniform material cohort required')
 index=next((i for i,c in enumerate(history.cohorts) if c.identity==cohort_id),None)
 if index is None:raise ValueError('Unknown material cohort')
 c=history.cohorts[index];state=history.state()[index];a=np.asarray(theta_y_triangle,dtype=float);theta_y_area=.5*np.linalg.det(np.column_stack((a[1]-a[0],a[2]-a[0])))
 if theta_y_area<=0:raise ValueError('Positive material-domain orientation required')
 t=float(c.thickness_m[0,0]);birth_radius=float(c.reference_mid_radius_m[0,0]);current_radius=float((state['inner_radius_m'][0,0]+state['outer_radius_m'][0,0])/2)
 reference=cylindrical_nodes(a,birth_radius,t);current=cylindrical_nodes(a,current_radius,t);mass=float(c.material_mass_kg[0,0])*theta_y_area/(history.da*history.dz)
 provisional=QuadraticPrism(reference,E,nu,1.,nq);rho=mass/provisional.reference_volume;cell=QuadraticPrism(reference,E,nu,rho,nq)
 return cell,current,{'cohort_id':c.identity,'birth_time_s':c.birth_time_s,'current_time_s':history.time_s,'birth_mid_radius_m':birth_radius,'current_mid_radius_m':current_radius,'material_theta_y_area_rad_m':theta_y_area,'retained_material_mass_kg':mass,'analytic_birth_volume_m3':birth_radius*t*theta_y_area,'FE_birth_volume_m3':cell.reference_volume,'geometry_quadrature_density_kg_m3':rho,'reference_geometry_frozen':True}
