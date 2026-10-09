import sys,unittest,numpy as np
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"source"))
from loop_chart import open_tube_chart
class TestOpenTube(unittest.TestCase):
 def test_side_seam(self):
  a,L=open_tube_chart([[0,0,0],[1,0,0],[3,0,0]],8);side=a[:-16].reshape(2,8,4,2)
  self.assertTrue(np.allclose(np.ptp(side,axis=2)[...,0],1/8));self.assertAlmostEqual(L,3)
  self.assertTrue(np.allclose(side[0,:,0,1],0));self.assertTrue(np.allclose(side[1,:,0,1],1/3))
  self.assertTrue(np.allclose(side[:,:,1,0]%1,np.roll(side[:,:,0,0],-1,axis=1)%1))
 def test_caps_nonzero_opposite_area(self):
  a,_=open_tube_chart([[0,0,0],[1,0,0]],8)
  caps=a[-16:].reshape(2,8,2).astype(float)
  areas=.5*np.sum(caps[:,:,0]*np.roll(caps[:,:,1],-1,axis=1)-caps[:,:,1]*np.roll(caps[:,:,0],-1,axis=1),axis=1)
  self.assertTrue(np.allclose(areas,[-np.sqrt(2)/2,np.sqrt(2)/2],atol=1e-7))
 def test_geometric_arc_length(self):
  a,L=open_tube_chart([[0,0,0],[1,1,0],[1,3,0]],8)
  self.assertAlmostEqual(L,np.sqrt(2)+2);self.assertAlmostEqual(float(a[8*4,1]),np.sqrt(2)/L,places=6)
 def test_degenerate_path_rejected(self):
  with self.assertRaises(ValueError):open_tube_chart([[0,0,0],[0,0,0]])
if __name__=="__main__":unittest.main()

