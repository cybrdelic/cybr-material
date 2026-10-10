import unittest,sys,json,hashlib,struct,zlib
from pathlib import Path
import numpy as np
from victorville_ground import evaluate,DEFAULT,environment

class GroundContracts(unittest.TestCase):
 def test_deterministic_and_shared_state(self):
  a,ma=evaluate({'tile_m':.12},n=64);b,mb=evaluate({'tile_m':.12},n=64)
  self.assertEqual(ma['packing'],mb['packing'])
  for k in a:np.testing.assert_array_equal(a[k],b[k])
  self.assertEqual(a['BaseColor'].shape,(64,64,3));self.assertTrue(np.isfinite(a['Height']).all())
  self.assertTrue(np.all(a['Metallic']==0));self.assertTrue(np.all(a['Opacity']==1));self.assertTrue(np.all(a['AO']==1))
  self.assertGreater(a['Height'].std(),0);self.assertGreater(a['Roughness'].std(),0)
  self.assertTrue(np.all(a['BaseColor']>.05))
 def test_primitive_identity_independent_of_resolution(self):
  _,a=evaluate({'tile_m':.12},n=64);_,b=evaluate({'tile_m':.12},n=96)
  self.assertEqual(a['packing']['primitive_sha256'],b['packing']['primitive_sha256'])
 def test_seed_changes_primitives(self):
  _,a=evaluate({'tile_m':.12},n=64);_,b=evaluate({'seed':92395,'tile_m':.12},n=64)
  self.assertNotEqual(a['packing']['primitive_sha256'],b['packing']['primitive_sha256'])
 def test_sediment_physically_buries_fixed_gravel(self):
  a,ma=evaluate({'tile_m':.12,'sediment_depth_m':.0005},n=128)
  b,mb=evaluate({'tile_m':.12,'sediment_depth_m':.004},n=128)
  np.testing.assert_array_equal(a['SupportHeight'],b['SupportHeight'])
  self.assertEqual(ma['packing']['primitive_sha256'],mb['packing']['primitive_sha256'])
  self.assertLess(mb['sediment']['exposed_gravel_coverage'],ma['sediment']['exposed_gravel_coverage'])
  self.assertTrue(np.all(b['Height']+1e-7>=b['SoilHeight']))
  self.assertGreater(mb['sediment']['soil_sand_dust_coverage'],.9)
 def test_invalid_physical_extent_rejected(self):
  for parameters in ({'tile_m':0},{'height_scale_m':.001},{'gravel_density':0}):
   with self.assertRaises(ValueError):evaluate(parameters,n=64)
 def test_environment_is_exactly_periodic(self):
  x=np.linspace(0,.6,100);y=x*.34
  a=environment(x,y,.6,92394);b=environment(x+.6,y+.6,.6,92394)
  for x,y in zip(a,b):np.testing.assert_allclose(x,y,atol=1e-12)

if __name__=='__main__':unittest.main()
