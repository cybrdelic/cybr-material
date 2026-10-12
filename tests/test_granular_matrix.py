"""Headless contracts for the opt-in matrix profile; selected recipes untouched."""
import importlib.util,json,tempfile,unittest
from pathlib import Path
R=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('granular_matrix',R/'source/granular_matrix.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
class GranularMatrixContract(unittest.TestCase):
 def test_metric_profile(self):
  p=module.load_profile(R/'source/granular_matrix_profiles.json');bands=p['bands'];self.assertEqual(p['units'],'metres');self.assertEqual(len(bands),5);self.assertEqual([b['wavelength_m']for b in bands],sorted([b['wavelength_m']for b in bands],reverse=True));self.assertTrue(all(b['height_coefficient_m']<b['wavelength_m']for b in bands));self.assertLess(sum(b['height_coefficient_m']for b in bands),.005)
 def test_nonfinite_is_rejected(self):
  d=json.loads((R/'source/granular_matrix_profiles.json').read_text());d['profiles']['victorville_du5_optical_matrix']['bands'][0]['height_coefficient_m']=float('nan')
  with tempfile.TemporaryDirectory()as folder:
   p=Path(folder)/'bad.json';p.write_text(json.dumps(d))
   with self.assertRaises(AssertionError):module.load_profile(p)
 def test_scene_authority_is_not_embedded(self):
  text=(R/'source/granular_matrix.py').read_text();self.assertNotIn('du_shelter',text);self.assertNotIn('du_compact',text);self.assertNotIn('Graded shoulder',text);self.assertNotIn('bpy.data',text)
if __name__=='__main__':unittest.main()
