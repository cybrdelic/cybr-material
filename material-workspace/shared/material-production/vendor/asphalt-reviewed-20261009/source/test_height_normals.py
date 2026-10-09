import unittest,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'base_recipe/source'))
from height_normals import height_normals,encode16
class TestNormals(unittest.TestCase):
 def test_flat(self):self.assertTrue(np.array_equal(height_normals(np.ones((5,5)),tile_m=.25,height_scale_m=.004),np.broadcast_to([0,0,1],(5,5,3))))
 def test_tilt_with_finite_edges(self):
  y,x=np.indices((9,13));h=.2*x/13-.3*y/9;goal=np.array([-.2,-.3,1]);goal/=np.linalg.norm(goal);self.assertLess(np.max(abs(height_normals(h,tile_m=1,height_scale_m=1,periodic=False)-goal)),1e-14)
 def test_periodic_sine(self):
  n=256;x=(np.arange(n)+.5)/n;h=np.broadcast_to(np.sin(2*np.pi*x),(4,n));a=height_normals(h,tile_m=1,height_scale_m=.01);s=-a[...,0]/a[...,2];self.assertLess(np.max(abs(s-.01*2*np.pi*np.cos(2*np.pi*x))),7e-6)
 def test_scale_invariance(self):
  h=np.random.default_rng(13).normal(size=(5,5));self.assertTrue(np.allclose(height_normals(h,tile_m=.2,height_scale_m=.004),height_normals(h,tile_m=.4,height_scale_m=.008)))
 def test_invalid(self):
  for h,t,z in [(np.ones((2,3)),1,1),(np.full((3,3),np.nan),1,1),(np.ones((3,3)),0,1),(np.ones((3,3)),1,-1)]:
   with self.assertRaises(ValueError):height_normals(h,tile_m=t,height_scale_m=z)
 def test_quantization(self):
  n=height_normals(np.random.default_rng(23).normal(size=(9,9)),tile_m=.25,height_scale_m=.004);q=encode16(n).astype(float)/65535*2-1;q/=np.linalg.norm(q,axis=-1,keepdims=True);self.assertLess(np.max(np.degrees(np.arccos(np.clip(np.sum(q*n,axis=-1),-1,1)))),.002)
if __name__=='__main__':unittest.main()
