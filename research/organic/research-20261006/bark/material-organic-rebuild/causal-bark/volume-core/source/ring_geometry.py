"""Closed circumferential material geometry with explicit radial birth history.
This supplies geometry/material identities; it does not yet supply solved ring
fracture/contact. Angular closure avoids treating a free flat edge as a tree.
"""
import sys,numpy as np
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'accretion/source'))
from phellem_history import AccretionHistory
from birth_cell import make_cohort_cell

def build(sectors=16,height=.008,stem_increment_m=.001):
 if int(sectors)!=sectors or sectors<4:raise ValueError('At least4 complete angular sectors required')
 h=AccretionHistory(.105,2*np.pi,height);h.advance(1,new_phellem_thickness_m=.0003,new_density_kg_m3=240,identity='older_outer');h.advance(2,stem_increment_m,new_phellem_thickness_m=.00025,new_density_kg_m3=240,identity='younger_inner');cells=[];q=[];metadata=[];triangles=[]
 for layer,identity in enumerate(('younger_inner','older_outer')):
  for sector in range(sectors):
   a=2*np.pi*sector/sectors;b=2*np.pi*(sector+1)/sectors
   for side,triangle in enumerate(([[a,-height/2],[b,-height/2],[b,height/2]],[[a,-height/2],[b,height/2],[a,height/2]])):
    cell,current,meta=make_cohort_cell(h,identity,np.array(triangle),2.8e6,.28);meta.update(layer=layer,sector=sector,triangle=side);cells.append(cell);q.append(current);metadata.append(meta);triangles.append(np.array(triangle))
 return h,cells,np.array(q),metadata,np.array(triangles)
