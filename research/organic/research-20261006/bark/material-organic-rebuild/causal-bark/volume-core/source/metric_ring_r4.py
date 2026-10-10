"""Retained-birth ring with material-normal numerical continuity scaling."""
import numpy as np
from incremental_ring_r3 import IncrementalRing
from ring_interfaces_r2 import build_interfaces
from native_cohesion_r2 import NativeCohesion
class MetricRing(IncrementalRing):
 def __init__(self,sectors=8,stem_increment_m=0.,continuity_compliance=.01):
  super().__init__(sectors,stem_increment_m)
  old=self.ops;new=build_interfaces(self.cells,self.origin,self.metadata,sectors,compliance=continuity_compliance)
  for key in ('pairs','J','T1','T2','area','normal_sign'):
   if not np.array_equal(old[key],new[key]):raise RuntimeError('Unexpected origin/operator change during penalty replacement')
  self.ops=new;self.bonds=NativeCohesion(self.N,new);self.continuity_compliance=float(continuity_compliance)
 def maximum_reference_frequency(self):
  if not np.array_equal(self.bonds.ops['kn'],self.bonds.ops['kt']):raise ValueError('Reference cohesive tangent requires equal normal/shear stiffness; anisotropic tangent not implemented')
  return super().maximum_reference_frequency()
