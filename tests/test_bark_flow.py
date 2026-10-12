import unittest,json,importlib.util,math,tempfile
from pathlib import Path
R=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('bark',R/'source/bark_flow.py');bark=importlib.util.module_from_spec(spec);spec.loader.exec_module(bark)
class BarkContracts(unittest.TestCase):
 def test_metric_profile(self):
  p=bark.load_profile(R/'source/bark_flow_profiles.json');self.assertEqual(len(p['bands']),2);self.assertEqual(p['height_units'],'world metres')
 def test_invalid_profile_fails(self):
  v=json.loads((R/'source/bark_flow_profiles.json').read_text());v['profiles']['victorville_dv3_dry_bark']['bands'][0]['height_m']=float('nan')
  with tempfile.TemporaryDirectory()as d:
   q=Path(d)/'invalid.json';q.write_text(json.dumps(v))
   with self.assertRaises(AssertionError):bark.load_profile(q)
 def test_no_scene_authority(self):
  s=(R/'source/bark_flow.py').read_text();self.assertNotIn('bpy.data',s);self.assertNotIn('CK pocket',s);self.assertNotIn('matrix_world',s)
if __name__=='__main__':unittest.main()
