import json,os,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
from materials.path_config import remap_path
from materials.core import resolve,ContractError
class Paths(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory();self.path=Path(self.tmp.name)/'paths.json';self.path.write_text(json.dumps({'schema_version':1,'prefixes':{'/workspace/shared':'/restored/shared','/workspace/shared/special':'/alternate/special'}}));self.env=patch.dict(os.environ,{'MATERIAL_PATH_MAP_FILE':str(self.path)});self.env.start()
 def tearDown(self):self.env.stop();self.tmp.cleanup()
 def test_most_specific_prefix_wins(self):self.assertEqual(str(remap_path('/workspace/shared/special/map.png')),'/alternate/special/map.png')
 def test_normal_prefix_relocation(self):self.assertEqual(str(remap_path('/workspace/shared/other/map.png')),'/restored/shared/other/map.png')
 def test_unmatched_path_is_unchanged(self):self.assertEqual(str(remap_path('/workspace/sharedish/a')),'/workspace/sharedish/a')
 def test_traversal_is_rejected(self):
  with self.assertRaises(ValueError):remap_path('/workspace/shared/../other')
  with self.assertRaises(ContractError):resolve('../outside')
 def test_registry_reference_uses_same_map(self):self.assertEqual(str(resolve('shared/special/map.png')),'/alternate/special/map.png')
if __name__=='__main__':unittest.main()
