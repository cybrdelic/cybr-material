"""Explicit candidate correction: two material-side histories, each half weighted.
The double normal functional is unchanged; its two equivalent quadrature
orientations are averaged. This does not replace the failed one-sided control.
"""
import numpy as np
from coupled_spans import Spans,exp,log,right_jacobian,S,UNIT
C=exp([0,0,-np.pi/2]);D=C.T

def swap(q):return np.r_[q[6:9],log(exp(q[9:12])@C),q[:3],log(exp(q[3:6])@D)]
def reorder(w):return np.r_[w[6:],w[:6]]
class SymmetricSpans(Spans):
 def state(self,q,elastic_a=None,elastic_b=None):
  a=Spans.state(self,q,elastic_a);b=Spans.state(self,swap(q),elastic_b);return dict(q=np.array(q).copy(),a=a,b=b,identity='symmetric-half-v1:'+self.identity)
 def exchange_state(self,old):return self.state(swap(old['q']),old['b']['elastic'],old['a']['elastic'])
 def trial(self,q,old):
  if old['identity']!='symmetric-half-v1:'+self.identity:raise ValueError('Wrong two-sided history identity')
  a=Spans.trial(self,q,old['a']);b=Spans.trial(self,swap(q),old['b']);w=.5*(a['wrenches']+reorder(b['wrenches']));wn=.5*(a['normal_wrenches']+reorder(b['normal_wrenches']));wf=.5*(a['friction_wrenches']+reorder(b['friction_wrenches']));g=dict(a['state']['geometry']);g['E']=.5*(a['normal_energy_J']+b['normal_energy_J']);g['normal_wrenches']=wn;A=g['Ra']@right_jacobian(q[3:6]);B=g['Rb']@right_jacobian(q[9:12]);res=np.r_[S*w[:3],A.T@w[3:6],S*w[6:9],B.T@w[9:12]]/UNIT;new=dict(q=np.array(q).copy(),a=a['state'],b=b['state'],identity=old['identity'])
  return dict(residual=res,state=new,geometry=g,wrenches=w,normal_wrenches=wn,friction_wrenches=wf,normal_energy_J=g['E'],tangential_energy_J=.5*(a['tangential_energy_J']+b['tangential_energy_J']),friction_heat_J=.5*(a['friction_heat_J']+b['friction_heat_J']),numerical_loss_J=.5*(a['numerical_loss_J']+b['numerical_loss_J']),tangential_endpoint_work_J=.5*(a['tangential_endpoint_work_J']+b['tangential_endpoint_work_J']),tangential_balance_J=.5*(a['tangential_balance_J']+b['tangential_balance_J']),contact_pa=np.concatenate([a['contact_pa'],b['contact_pb']]),contact_pb=np.concatenate([a['contact_pb'],b['contact_pa']]),tau=.5*np.concatenate([a['tau'],-b['tau']]),du=np.concatenate([a['du'],-b['du']]),active_set=np.r_[a['active_set'],b['active_set']],sliding_set=np.r_[a['sliding_set'],b['sliding_set']])
