import copy,io,json,struct,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
from PIL import Image
from materials.core import ROOT,workspace,sha,ContractError
from materials.png_fingerprint import PNGSink,classify
from materials.experiments import experiment,verify_experiment_build

class PNGContracts(unittest.TestCase):
 def png(self,level=6,value=42):
  out=io.BytesIO();Image.fromarray(np.full((64,64),value,np.uint16)).save(out,format='PNG',compress_level=level);return out.getvalue()
 def fingerprint(self,b):
  sink=PNGSink()
  for i in range(0,len(b),7):sink.write(b[i:i+7])
  return sink.fingerprint()
 def test_fragmented_gray16_stream_matches_hash(self):
  b=self.png();f=self.fingerprint(b);self.assertEqual(f['IHDR'][:3],[64,64,16]);self.assertEqual(f['png_sha256'],__import__('hashlib').sha256(b).hexdigest())
 def test_lossless_deflate_change_is_not_pixel_failure(self):
  a=self.fingerprint(self.png(1));b=self.fingerprint(self.png(9))
  self.assertEqual(a['filtered_stream_sha256'],b['filtered_stream_sha256']);self.assertEqual(classify(a,b,True),'deflate_bitstream_only')
 def test_pixel_change_cannot_be_excused_as_codec(self):
  a=self.fingerprint(self.png(value=41));b=self.fingerprint(self.png(value=42));self.assertEqual(classify(a,b,False),'decoded_pixel_mismatch')
 def test_bad_crc_and_truncation_fail(self):
  b=bytearray(self.png());b[29]^=1
  with self.assertRaises(ContractError):self.fingerprint(b)
  with self.assertRaises(ContractError):self.fingerprint(self.png()[:-6])

class ExperimentContracts(unittest.TestCase):
 def setUp(self):
  self.temp=tempfile.TemporaryDirectory(dir=ROOT);self.root=Path(self.temp.name);self.exp=self.root/'experiments/example';self.exp.mkdir(parents=True);source=self.root/'candidate.blend';source.write_bytes(b'test-only')
  self.source=source;self.row={'schema_version':1,'experiment_id':'example','material_id':'32_concrete_polished','selected':False,'quality_acceptance':False,'source_scene':str(source.relative_to(workspace())),'source_sha256':sha(source),'dependencies':[],'capture':{},'purpose':'test','pass_gates':['explicit test']};self.save()
 def tearDown(self):self.temp.cleanup()
 def save(self):(self.exp/'experiment.json').write_text(json.dumps(self.row))
 def test_explicit_experiment_is_not_selected(self):
  with patch('materials.experiments.ROOT',self.root):
   r,p=experiment('example');self.assertFalse(r['selected'])
 def test_self_promotion_forbidden(self):
  self.row['selected']=True;self.save()
  with patch('materials.experiments.ROOT',self.root):
   with self.assertRaises(ContractError):experiment('example')
 def test_source_drift_and_output_override_forbidden(self):
  self.row['capture']['output']='overwrite.png';self.save()
  with patch('materials.experiments.ROOT',self.root):
   with self.assertRaises(ContractError):experiment('example')
  self.row['capture']={};self.save();self.source.write_bytes(b'changed')
  with patch('materials.experiments.ROOT',self.root):
   with self.assertRaises(ContractError):experiment('example')
 def test_name_cannot_escape_experiment_root(self):
  with self.assertRaises(ContractError):experiment('../registry')
if __name__=='__main__':unittest.main()
