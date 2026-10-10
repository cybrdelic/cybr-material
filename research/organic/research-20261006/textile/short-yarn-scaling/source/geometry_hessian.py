from loop_bundle import *
def geometry_hessian(rod,w,C,eps=1e-5):
 # The position map is separable in segment frames, so these blocks are local.
 def gradient(a):
  R=exp(a);body=mv(np.swapaxes(R,-1,-2),C);g=rod.ds[:,None]*np.cross(np.array([0.,0.,1.]),body);return mv(np.swapaxes(right_jacobian(a),-1,-2),g)
 H=np.empty((rod.M,3,3))
 for component in range(3):
  d=np.zeros_like(w);d[:,component]=eps;H[:,:,component]=(gradient(w+d)-gradient(w-d))/(2*eps)
 return (H+np.swapaxes(H,-1,-2))*.5

