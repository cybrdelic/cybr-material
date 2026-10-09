import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from resource_contract import resource_violations
class ResourceTests(unittest.TestCase):
 def base(self):return {'RSS_cap_MiB':5120,'peak_process_tree_RSS_MiB':5000,'kernel_max_child_RSS_MiB':5100,'exit_code':0}
 def test_success(self):self.assertEqual(resource_violations(self.base(),5120),[])
 def test_exact_boundary(self):r=self.base();r['kernel_max_child_RSS_MiB']=5120;self.assertEqual(resource_violations(r,5120),[])
 def test_between_poll_peak(self):r=self.base();r['kernel_max_child_RSS_MiB']=5310;self.assertTrue(resource_violations(r,5120))
 def test_sampled_peak(self):r=self.base();r['peak_process_tree_RSS_MiB']=5200;self.assertTrue(resource_violations(r,5120))
 def test_missing_peak(self):r=self.base();del r['kernel_max_child_RSS_MiB'];self.assertTrue(resource_violations(r,5120))
 def test_nonfinite(self):r=self.base();r['kernel_max_child_RSS_MiB']=float('nan');self.assertTrue(resource_violations(r,5120))
 def test_cap_mismatch(self):self.assertTrue(resource_violations(self.base(),6144))
 def test_process_failure(self):r=self.base();r['exit_code']=-15;self.assertTrue(resource_violations(r,5120))
if __name__=='__main__':unittest.main()
