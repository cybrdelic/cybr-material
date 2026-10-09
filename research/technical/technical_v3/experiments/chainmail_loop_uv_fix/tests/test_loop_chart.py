import sys, unittest, numpy as np
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"source"))
from loop_chart import ring_loop_chart

class ClosedChart(unittest.TestCase):
    def test_constant_area_and_no_long_quad(self):
        for nt,ns in [(64,16),(17,7),(8,8)]:
            a=ring_loop_chart(nt,ns).astype(float)
            self.assertTrue(np.allclose(np.ptp(a,axis=1),[1/ns,1/nt],atol=1e-7))
            area=.5*np.sum(a[:,:,0]*np.roll(a[:,:,1],-1,axis=1)-a[:,:,1]*np.roll(a[:,:,0],-1,axis=1),axis=1)
            self.assertTrue(np.allclose(area,-1/ns/nt,atol=1e-7))
            self.assertAlmostEqual(float(area.sum()),-1,places=6)
    def test_seams_share_periodic_samples(self):
        nt,ns=64,16; a=ring_loop_chart(nt,ns).reshape(nt,ns,4,2)
        # Geometrically coincident corners agree modulo an integer texture turn.
        self.assertTrue(np.allclose(a[:,:,1] % 1,np.roll(a[:,:,0],-1,axis=0) % 1))
        self.assertTrue(np.allclose(a[:,:,3] % 1,np.roll(a[:,:,0],-1,axis=1) % 1))
    def test_exact_boundary_corners(self):
        a=ring_loop_chart()
        self.assertTrue(np.array_equal(a[-1],[[15/16,63/64],[15/16,1],[1,1],[1,63/64]]))
    def test_exposes_legacy_defect(self):
        nt,ns=64,16
        old=np.array([[(j/ns,i/(nt-1)),(j/ns,((i+1)%nt)/(nt-1)),(((j+1)%ns)/ns,((i+1)%nt)/(nt-1)),(((j+1)%ns)/ns,i/(nt-1))] for i in range(nt) for j in range(ns)])
        bad=np.any(np.ptp(old,axis=1)>[1/ns+1e-8,1/(nt-1)+1e-8],axis=1)
        self.assertEqual(int(bad.sum()),79)
if __name__=="__main__": unittest.main()

