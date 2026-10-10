from pathlib import Path
import numpy as np
from carriage import Carriage,RodHessianOperator,direction,solve,P
class TangentSuspended(Carriage):
 def predictor_boundary(self,delta_theta):
  x=self.x.copy();theta=self.theta;_,g=self.evaluate(x);pg,lam,J,S=self.projection(g);op=RodHessianOperator(self.c,lam);eps=1e-6;gg=[];hh=[]
  for sign in [1,-1]:
   self.theta=theta+sign*eps;_,gp=self.evaluate(x);gg.append(gp+self.jacobian().T@lam);hh.append(self.constraints().copy())
  self.theta=theta;self.evaluate(x);control_gradient=(gg[0]-gg[1])/(2*eps);control_constraint=(hh[0]-hh[1])/(2*eps);dx,check=direction(op,self,control_gradient*delta_theta,control_constraint*delta_theta,0.)
  self.predictor_diagnostic=dict(method='Implicit equilibrium tangent at fixed suspended boundary and constant axial force',control_derivative_epsilon_rad=eps,linear_check=check,maximum_rotation_prediction_rad=float(np.linalg.norm(dx[:-1].reshape(self.c.w.shape),axis=2).max()),carriage_prediction_m=float(dx[-1]*self.L))
  if check['relative_linear_residual']>1e-7:raise RuntimeError('Equilibrium tangent solve failed')
  return x+dx
