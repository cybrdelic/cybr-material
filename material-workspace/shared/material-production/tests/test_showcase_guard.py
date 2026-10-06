import unittest
from materials.showcase import preferred_reference,validate_visual_promotion
from materials.core import ContractError

class ShowcaseGuard(unittest.TestCase):
 def setUp(self):
  self.m={'preferred_visual_reference':{'kind':'complete_panel','status':'restored_visual_baseline','library_file_id':'known-baseline'}}
  self.c={'kind':'complete_panel','library_file_id':'known-candidate'}
  self.r={'baseline_library_file_id':'known-baseline','candidate_library_file_id':'known-candidate','result':'clear_visual_improvement','native_pixels_reviewed':True,'complete_view_compared':True,'review_receipt':'checked.json','review_receipt_sha256':'known-receipt-sha'}
 def test_complete_matched_review(self):self.assertTrue(validate_visual_promotion(self.m,self.c,self.r))
 def test_numerical_pass_does_not_promote(self):
  self.r.update(result='numerical_checks_passed')
  with self.assertRaises(ContractError):validate_visual_promotion(self.m,self.c,self.r)
 def test_detail_does_not_replace_panel(self):
  self.c['kind']='detail_study'
  with self.assertRaises(ContractError):validate_visual_promotion(self.m,self.c,self.r)
 def test_wrong_baseline_rejected(self):
  self.r['baseline_library_file_id']='weaker-experiment'
  with self.assertRaises(ContractError):validate_visual_promotion(self.m,self.c,self.r)
 def test_candidate_mismatch_rejected(self):
  self.r['candidate_library_file_id']='unreviewed'
  with self.assertRaises(ContractError):validate_visual_promotion(self.m,self.c,self.r)
 def test_pinned_reference_wins_over_selected_diagnostic(self):
  self.m['selected']={'appearance_anchor':{'library_file_id':'later-detail'}}
  self.assertEqual(preferred_reference(self.m)['library_file_id'],'known-baseline')

if __name__=='__main__':unittest.main()
