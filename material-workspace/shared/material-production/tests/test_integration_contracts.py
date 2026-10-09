import copy,importlib.util,json,os,subprocess,sys,tempfile,unittest
from pathlib import Path
import numpy as np
from materials.core import *
from materials.generation import evaluate,validate_fields

class Contracts(unittest.TestCase):
 def test_registry_35_unique_and_leather_explicitly_unqualified(self):
  self.assertEqual(len(registry()['materials']),35)
  leather=material('08');validate_selection(leather)
  self.assertIn('original procedural appearance reference',leather['quality'])
  self.assertIn('leather_coupon_transfer',leather['selected']['scene'])
  self.assertFalse(leather['reviewed_correction']['formation_simulation_qualified'])
  self.assertFalse(leather['reviewed_correction']['full_material_realism_accepted'])
 def test_no_revision_discovery_or_auto_promotion(self):
  self.assertNotIn('mineral-intergrowth',material('32')['selected']['scene'])
  self.assertIn('top_only_normal_candidate',material('05')['selected']['scene'])
 def test_fleece_wood_anchors_protected(self):
  for i in ['03','04','17','24']:self.assertTrue(material(i)['protected_appearance'])
  self.assertIn('r5_hero',material('17')['selected']['scene'])
 def test_reference_cannot_escape_workspace(self):
  with self.assertRaises(ContractError):resolve('../etc/passwd')
 def test_hash_mismatch_fails_closed(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'source';p.write_bytes(b'old');h=sha(p);p.write_bytes(b'new')
   with self.assertRaises(ContractError):check_hash(p,h)
 def test_output_cannot_overwrite(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'result.json';write_new(p,{'a':1})
   with self.assertRaises(FileExistsError):write_new(p,{'a':2})
 def test_import_has_no_process_or_file_side_effects(self):
  p=subprocess.run([sys.executable,'-c','import materials, materials.core, materials.generation; print("imported")'],capture_output=True,text=True,check=True)
  self.assertEqual(p.stdout.strip(),'imported')
 def test_brass_actual_source_generation_channels_and_units(self):
  a,r=evaluate('05_champagne_brass',32);b,_=evaluate('05_champagne_brass',32)
  self.assertIn('recovered exact selected brass source',r['mode']);self.assertEqual(r['tile_m'],.2)
  for k in a:np.testing.assert_array_equal(a[k],b[k])
  self.assertGreater(float(np.std(a['Height'])),0)
  self.assertEqual(a['BaseColor'].shape,(32,32,3))
 def test_concrete_actual_selected_generation_deterministic(self):
  a,r=evaluate('32_concrete_polished',32);b,_=evaluate('32_concrete_polished',32)
  for k in a:np.testing.assert_array_equal(a[k],b[k])
  self.assertEqual(r['tile_m'],.32);self.assertEqual(r['height_scale_m'],.002)
  self.assertTrue(np.all(a['Metallic']==0));self.assertTrue(np.all(a['Opacity']==1))
  self.assertGreater(float(np.std(a['BaseColor'])),0)
  self.assertEqual(r['recipe']['clipped_height_pixels'],0)
 def test_channel_nan_range_and_shape_rejected(self):
  d={k:np.ones((16,16),np.float32) for k in ['Height','Roughness','Metallic','Opacity']};d['BaseColor']=np.ones((16,16,3),np.float32)
  validate_fields(d,16)
  for val in [float('nan'),-1,2]:
   bad=copy.deepcopy(d);bad['Height'][0,0]=val
   with self.assertRaises(ContractError):validate_fields(bad,16)
 def test_smoke_resolution_cannot_claim_4k(self):
  with self.assertRaises(ContractError):evaluate('05_champagne_brass',4096)
 def test_unqualified_generator_never_falls_back(self):
  with self.assertRaises(ContractError):evaluate('17_blanket',64)
 def test_actual_brass_and_concrete_inspection_contracts(self):
  for key,cid in [('brass','05'),('concrete','32')]:
   d=read(Path(__file__).parent/'fixtures'/f'{key}_inspection.json')
   validate_inspection(material(cid),d)
 def test_wrong_units_and_uv_cannot_build(self):
  d=read(Path(__file__).parent/'fixtures/brass_inspection.json')
  for key,value in [('unit_scale',.001)]:
   bad=copy.deepcopy(d);bad[key]=value
   with self.assertRaises(ContractError):validate_inspection(material('05'),bad)
  for key,value in [('scale',[2,1,1]),('uv_layers',[]),('dimensions_m',[.4,.2,.003]),('material_face_counts',{'0':66561})]:
   bad=copy.deepcopy(d);bad['objects'][0][key]=value
   with self.assertRaises(ContractError):validate_inspection(material('05'),bad)
 def test_double_relief_normal_mode_is_rejected(self):
  d=read(Path(__file__).parent/'fixtures/brass_inspection.json')
  d['materials'][1]['normal_nodes'][0]['space']='TANGENT'
  with self.assertRaises(ContractError):validate_inspection(material('05'),d)
 def test_source_import_namespace_is_restored(self):
  before_path=sys.path[:];sentinel=object();prior=sys.modules.get('surface_math')
  sys.modules['surface_math']=sentinel
  try:
   evaluate('05_champagne_brass',16)
   self.assertIs(sys.modules['surface_math'],sentinel);self.assertEqual(sys.path,before_path)
  finally:
   if prior is None:sys.modules.pop('surface_math',None)
   else:sys.modules['surface_math']=prior
if __name__=='__main__':unittest.main()
