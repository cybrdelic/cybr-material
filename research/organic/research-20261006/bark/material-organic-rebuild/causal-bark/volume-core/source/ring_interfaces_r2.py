"""Numerical-continuity penalty scaled in the retained material face-normal metric.

The ring's cylindrical material surface is developed to (R_birth*theta, y).
This intrinsic metric avoids using a curved volume's chord centroid as a normal
thickness. It is specific to these uniform cylindrical cohorts, not a general
curved-polyhedron trace estimate. Finite eta still adds O(eta) compliance and
must be independently refined. No anatomical interface stiffness is inferred.
"""
import numpy as np
from ring_interfaces import build_interfaces as original_build

def _material_triangle(meta,sectors,height):
 a=2*np.pi*meta['sector']/sectors;b=2*np.pi*(meta['sector']+1)/sectors
 return np.array([[a,-height/2],[b,-height/2],[b,height/2]] if meta['triangle']==0 else [[a,-height/2],[b,height/2],[a,height/2]])

def build_interfaces(cells,q0,metadata,sectors,order=3,compliance=.01,strength=1e5,Gc=10.):
 if not np.isfinite(compliance) or compliance<=0:raise ValueError('Positive dimensionless continuity compliance required')
 o=original_build(cells,q0,metadata,sectors,order,compliance,strength,Gc);cache={};length=np.empty(len(o['area']));side_lengths=np.empty((len(o['area']),2));modulus=np.empty(len(o['area']))
 for p,(i,j) in enumerate(o['pairs']):
  key=(int(i),int(j))
  if key not in cache:
   a,b=metadata[i],metadata[j];kind=o['kind'][p]
   if kind=='birth_interface':
    da=.5*a['analytic_birth_volume_m3']/(a['material_theta_y_area_rad_m']*a['birth_mid_radius_m']);db=.5*b['analytic_birth_volume_m3']/(b['material_theta_y_area_rad_m']*b['birth_mid_radius_m'])
   elif kind=='fracture':
    if a['cohort_id']!=b['cohort_id']:raise ValueError('Intralayer face spans different birth metrics')
    ha=np.ptp(cells[i].X[:,1]);hb=np.ptp(cells[j].X[:,1]);A=_material_triangle(a,sectors,ha);B=_material_triangle(b,sectors,hb);B[:,0]+=2*np.pi*round((A[:,0].mean()-B[:,0].mean())/(2*np.pi));radius=a['birth_mid_radius_m'];A[:,0]*=radius;B[:,0]*=radius
    shared=[]
    for vertex in A:
     if np.min(np.linalg.norm(B-vertex,axis=1))<1e-12:shared.append(vertex)
    if len(shared)!=2:raise ValueError('Expected exactly one shared developed material edge')
    edge=shared[1]-shared[0];normal=np.array([-edge[1],edge[0]])/np.linalg.norm(edge);mid=.5*(shared[0]+shared[1]);da=abs(float((A.mean(0)-mid)@normal));db=abs(float((B.mean(0)-mid)@normal))
   else:raise ValueError('Unknown interface role')
   if not np.isfinite([da,db]).all() or min(da,db)<=0:raise ValueError('Invalid material dual thickness')
   # Isotropic equal traction penalties retain the original bilinear strength
   # convention. Direction-dependent physical K/strength is a separate model.
   inv_stiffness=compliance*(da/cells[i].E+db/cells[j].E);cache[key]=(1/inv_stiffness,da,db)
  K,da,db=cache[key];o['kt'][p]=o['kn'][p]=K;length[p]=da+db;side_lengths[p]=(da,db);modulus[p]=K*compliance*(da+db)
 o.update(dual_normal_length_m=length,dual_side_lengths_m=side_lengths,effective_penalty_modulus_Pa=modulus,numerical_continuity_compliance=np.full(len(length),compliance),interface_role=np.full(len(length),'numerical_continuity'))
 return o
