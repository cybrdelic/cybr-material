import copy,json,sys,unittest
from pathlib import Path
import numpy as np
from materials.core import ROOT,ContractError,material,validate_inspection,read
from materials.generation import evaluate

class AsphaltSelection(unittest.TestCase):
 def test_reviewed_selection_and_prior_retained(self):
  m=material('28')
  self.assertEqual(m['selected']['scene_sha256'],'2396d4c9cb4291548106b04891d2c8a9b6f5379128a57e3654ed36e79ed8984e')
  self.assertEqual(m['selected']['version'],'reviewed_fractured_grain_20261009')
  self.assertEqual(m['selected']['capture']['camera'],'Capture / detail')
  self.assertEqual(m['prior_selections'][0]['scene_sha256'],'5ea7091804d7b5e9acf71b7367580237cb0560a3d51e79598f7aa01ad48a0d0d')
  self.assertFalse(m['reviewed_correction']['subsequent_binder_refinement_selected'])
  self.assertFalse(m['reviewed_correction']['formation_simulation_qualified'])
 def test_deterministic_finite_fields_and_units(self):
  a,r=evaluate('28_asphalt',32);b,_=evaluate('28_asphalt',32)
  self.assertEqual((r['tile_m'],r['height_scale_m']),(.25,.004))
  self.assertEqual(r['recipe']['clipped_height_pixels'],0)
  for k in a:np.testing.assert_array_equal(a[k],b[k])
  self.assertTrue(np.all(a['Metallic']==0))
  self.assertTrue(np.all(a['Opacity']==1))
 def test_exact_actual_scene_inspection(self):
  d=read(ROOT/'tests/fixtures/asphalt_inspection.json');validate_inspection(material('28'),d)
  for update in ['normal','resolution','packing','scale','uv','identity']:
   bad=copy.deepcopy(d)
   if update=='identity':bad['render_contract_sha256']='wrong'
   if update=='normal':bad['materials'][0]['normal_nodes'][0]['space']='TANGENT'
   if update=='resolution':bad['images'][0]['resolution']=[1536,1536]
   if update=='packing':bad['images'][0]['packed']=False
   if update=='scale':bad['objects'][0]['scale']=[2,1,1]
   if update=='uv':bad['objects'][0]['uv_layers']=[]
   with self.assertRaises(ContractError):validate_inspection(material('28'),bad)
 def test_source_namespace_is_restored(self):
  before=sys.path[:];old={k:sys.modules.get(k) for k in ['surface_math','fracture_field']};sentinel=object()
  try:
   for k in old:sys.modules[k]=sentinel
   evaluate('28_asphalt',16)
   self.assertEqual(sys.path,before)
   for k in old:self.assertIs(sys.modules[k],sentinel)
  finally:
   for k,v in old.items():
    if v is None:sys.modules.pop(k,None)
    else:sys.modules[k]=v
if __name__=='__main__':unittest.main()
