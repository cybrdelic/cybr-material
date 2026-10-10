"""Closed-ring adjacency and birth-reference cohesive quadrature.
No tangential foundation surrogate. Interfaces are tissue-continuity test bonds,
not inferred weak glue planes at each annual birth.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from quadratic_prism import shapes

def build_interfaces(cells,q0,metadata,sectors,order=3,compliance=.01,strength=1e5,Gc=10.):
 pairs=[];J=[];T1=[];T2=[];area=[];kt=[];sgn=[];kinds=[]
 u,w=leggauss(order);u01=(u+1)/2;w01=w/2
 def add(c1,c2,parameters,kind,reference):
  # parameter rows: r1,s1,z1, r2,s2,z2, tangent1 directions, tangent2 directions,weight.
  xref=np.vstack((reference[0],reference[1]));h=.5*(np.ptp(cells[c1].X[:,1])+np.ptp(cells[c2].X[:,1])) if kind=='fracture' else .5*(metadata[c1]['analytic_birth_volume_m3']/metadata[c1]['material_theta_y_area_rad_m']/metadata[c1]['birth_mid_radius_m']+metadata[c2]['analytic_birth_volume_m3']/metadata[c2]['material_theta_y_area_rad_m']/metadata[c2]['birth_mid_radius_m'])
  K=2.8e6/(compliance*h)
  for r1,s1,z1,r2,s2,z2,d11,d12,d21,d22,weight in parameters:
   n1,d1=shapes(r1,s1,z1);n2,d2=shapes(r2,s2,z2);j=np.r_[-n1,n2];a=np.r_[.5*(d1@d11),.5*(d2@d12)];b=np.r_[.5*(d1@d21),.5*(d2@d22)];j-=j.mean();a-=a.mean();b-=b.mean();cross=np.cross(a@xref,b@xref);sign=1.
   if kind=='fracture':sign=1. if cross@(reference[1].mean(0)-reference[0].mean(0))>0 else -1.
   pairs.append((c1,c2));J.append(j);T1.append(a);T2.append(b);area.append(weight*np.linalg.norm(cross));kt.append(K);sgn.append(sign);kinds.append(kind)
 # Same angular triangulation in both cohorts; born interlayer faces coincide at q0.
 per_layer=sectors*2
 for k in range(per_layer):
  parameters=[]
  for r,wr in zip(u01,w01):
   for v,wv in zip(u01,w01):
    s=(1-r)*v;parameters.append((r,s,1.,r,s,-1.,np.array([1.,0.,0.]),np.array([1.,0.,0.]),np.array([0.,1.,0.]),np.array([0.,1.,0.]),wr*wv*(1-r)))
  add(k,k+per_layer,parameters,'birth_interface',(q0[k],q0[k+per_layer]))
 for layer in range(2):
  edge_map={};verts=[]
  for sector in range(sectors):
   a=sector;b=(sector+1)%sectors;verts.extend((((a,0),(b,0),(b,1)),((a,0),(b,1),(a,1))))
  for k,v in enumerate(verts):
   for ia,ib in ((0,1),(1,2),(2,0)):
    key=tuple(sorted((v[ia],v[ib])))
    if key not in edge_map:edge_map[key]=(k,ia,ib)
    else:
     j,ja,jb=edge_map.pop(key);ga,gb=verts[j][ja],verts[j][jb];ia2=v.index(ga);ib2=v.index(gb);e1=np.zeros(3);e2=np.zeros(3);e1[ja]=-1;e1[jb]=1;e2[ia2]=-1;e2[ib2]=1;d11=np.r_[e1[1:],0.];d12=np.r_[e2[1:],0.];parameters=[]
     for t,wt in zip(u01,w01):
      p1=np.zeros(3);p2=np.zeros(3);p1[ja]=1-t;p1[jb]=t;p2[ia2]=1-t;p2[ib2]=t
      for z,wz in zip(u,w):parameters.append((p1[1],p1[2],z,p2[1],p2[2],z,d11,d12,np.array([0.,0.,1.]),np.array([0.,0.,1.]),wt*wz))
     c1=layer*per_layer+j;c2=layer*per_layer+k;add(c1,c2,parameters,'fracture',(cells[c1].X,cells[c2].X))
  # Only the axial top/bottom faces may remain open; the circumferential seam closes.
  if any(a[1]!=b[1] for a,b in edge_map):raise AssertionError('An unpaired circumferential edge remains')
 return {'pairs':np.asarray(pairs,dtype=np.int32),'J':np.asarray(J),'T1':np.asarray(T1),'T2':np.asarray(T2),'area':np.asarray(area),'kt':np.asarray(kt),'kn':np.asarray(kt),'normal_sign':np.asarray(sgn),'kind':np.asarray(kinds),'strength':np.full(len(pairs),strength),'Gc':np.full(len(pairs),Gc)}
