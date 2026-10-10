"""Fixed material-side contact quadrature on three actual director rods.
Original double normal energy, symmetric histories, exact rod virtual power.
Unique connected target support required. No full bundle solve claim.
"""
import hashlib,numpy as np
from rod_pullback import RodCoordinates,L,S,exp,log,mv,skew
from segment_contact_sine import evaluate as native
KC=1e8;KT=2.5e9;MU=.3
class FlexibleContact(RodCoordinates):
 def __init__(self,n=3,segments=16,order=9,partner_order=24):
  super().__init__(n,segments);x,w=np.polynomial.legendre.leggauss(order);edge=np.linspace(0,L,129);self.s=(.5*(edge[:-1,None]+edge[1:,None])+.5*L/128*x).ravel();self.w=np.tile(.5*L/128*w,128);self.k=KT*self.w;self.x,self.wg=np.polynomial.legendre.leggauss(partner_order);self.pairs=[(i,j) for i in range(n) for j in range(n) if i!=j];self.identity=hashlib.sha256(self.s.tobytes()+self.w.tobytes()+self.rod.s.tobytes()+('geodesic-surface-quadratic-end-v3:'+str(self.pairs)).encode()).hexdigest()
 def side(self,g,i,j):
  pa,Ra,_,_=self.sample(g,i,self.s);a=g['p'][j,:-1];v=np.diff(g['p'][j],axis=0);v2=np.sum(v*v,axis=1);rel=pa[:,None,:]-a[None,:,:];raw=np.sum(rel*v[None,:,:],axis=2)/v2;lam=np.clip(raw,0,1);delta=rel-lam[:,:,None]*v[None,:,:];d2=np.sum(delta*delta,axis=2);near=np.argmin(d2,axis=1);near_u=lam[np.arange(len(pa)),near];sb=self.rod.s[near]+near_u*self.rod.ds[near];pb,Rb,_,_=self.sample(g,j,sb);dn=pa-pb;dd=np.linalg.norm(dn,axis=1)
  if np.min(dd)<2*self.radius*1e-12:raise ValueError('Centreline singularity')
  normal=dn/dd[:,None];ga=np.zeros_like(pa);gp=np.zeros_like(g['p']);E=0.;last=np.full(len(pa),-np.inf);components=np.zeros(len(pa),int);R=2*self.radius
  for k in range(self.M):
   ids=np.flatnonzero(d2[:,k]<R*R)
   if not len(ids):continue
   centre=raw[ids,k];perp=rel[ids,k]-centre[:,None]*v[k];half=np.sqrt(np.maximum((R*R-np.sum(perp*perp,axis=1))/v2[k],0));lo=np.maximum(0,centre-half);hi=np.minimum(1,centre+half);low_arc=self.rod.s[k]+lo*self.rod.ds[k];high_arc=self.rod.s[k]+hi*self.rod.ds[k];components[ids]+=(low_arc>last[ids]+L*1e-12);last[ids]=high_arc
   theta_lo=np.arcsin(np.clip((lo-centre)/half,-1,1));theta_hi=np.arcsin(np.clip((hi-centre)/half,-1,1));theta=.5*(theta_lo+theta_hi)[:,None]+.5*(theta_hi-theta_lo)[:,None]*self.x;u=centre[:,None]+half[:,None]*np.sin(theta);weights=self.w[ids,None]*self.rod.ds[k]*half[:,None]*np.cos(theta)*(.5*(theta_hi-theta_lo)[:,None]*self.wg);target=a[k]+u[:,:,None]*v[k];diff=pa[ids,None,:]-target;rr2=np.sum(diff*diff,axis=2);root=np.sqrt(np.maximum(R*R-rr2,0));z=root/R;term=np.empty_like(z);small=z<.02;z2=z*z;term[small]=R*z[small]*z2[small]*(1/3+z2[small]*(1/5+z2[small]*(1/7+z2[small]/9)));term[~small]=R*np.arccosh(R/np.sqrt(rr2[~small]))-root[~small];E+=float(np.sum(KC/np.pi*term*weights));force=(-KC/np.pi*root/rr2*weights)[:,:,None]*diff;ga[ids]+=force.sum(axis=1);gp[j,k]-=np.sum((1-u[:,:,None])*force,axis=(0,1));gp[j,k+1]-=np.sum(u[:,:,None]*force,axis=(0,1))
  if components.max(initial=0)>1:raise ValueError('Multiple disconnected partner contact components')
  dummy=np.zeros((self.n,self.M,3));self.point_load(g,i,self.s,ga,np.zeros_like(ga),gp,dummy);N=-np.sum(ga*normal,axis=1)
  if N.min(initial=0)<-1e-14:raise ValueError('Normal-pressure reduction is tensile')
  return dict(i=i,j=j,pa=pa,pb=pb,Ra=Ra,Rb=Rb,sb=sb,n=normal,N=np.maximum(N,0),E=E,gp=gp,segment=near,components=components)
 def geometry(self,q):
  g=self.unpack(q);sides=[self.side(g,i,j) for i,j in self.pairs];g['sides']=sides;g['normal_energy_J']=.5*sum(a['E'] for a in sides);g['normal_gp']=.5*sum((a['gp'] for a in sides),np.zeros_like(g['p']));return g
 def state(self,q,elastic=None):
  return dict(q=np.array(q).copy(),g=self.geometry(q),elastic=np.zeros((len(self.pairs),len(self.s),3)) if elastic is None else np.array(elastic).copy(),identity=self.identity)
 def trial(self,q,old):
  if old['identity']!=self.identity:raise ValueError('Material history mismatch')
  g=self.geometry(q);gp=np.zeros_like(g['p']);gm=np.zeros((self.n,self.M,3));e=np.zeros_like(old['elastic']);heat=loss=work=0.;surface=[];active=[];sliding=[];normal_q=self.pullback(g,g['normal_gp'],np.zeros_like(gm))
  for m,(a,b) in enumerate(zip(old['g']['sides'],g['sides'])):
   ids=np.flatnonzero((b['N']>0)|np.any(old['elastic'][m]!=0,axis=1));active.append(b['N']>0);slide=np.zeros(len(self.s),bool);sliding.append(slide)
   if not len(ids):continue
   i,j=b['i'],b['j'];Qi=b['Ra'][ids]@np.swapaxes(a['Ra'][ids],-1,-2);pb1,Rb1,_,_=self.sample(g,j,a['sb'][ids]);pb0,Rb0,_,_=self.sample(old['g'],j,b['sb'][ids]);Qjf=Rb1@np.swapaxes(a['Rb'][ids],-1,-2);Qjb=b['Rb'][ids]@np.swapaxes(Rb0,-1,-2);Qj=Qjf@exp(.5*log(np.swapaxes(Qjf,-1,-2)@Qjb));mean=Qi@exp(.5*log(np.swapaxes(Qi,-1,-2)@Qj));pred=mv(mean,a['n'][ids]);cross=np.cross(pred,b['n'][ids]);c=np.sum(pred*b['n'][ids],axis=1)
   if c.min(initial=1)<-.99:raise ValueError('Subdivide normal rotation')
   K=skew(cross);Q=(np.eye(3)+K+K@K/(1+c)[:,None,None])@mean;m0=.5*(a['pa'][ids]+a['pb'][ids]);m1=.5*(b['pa'][ids]+b['pb'][ids]);forward=b['pa'][ids]+mv(Qi,m0-a['pa'][ids])-pb1-mv(Qjf,m0-a['pb'][ids]);back=a['pa'][ids]+mv(np.swapaxes(Qi,-1,-2),m1-b['pa'][ids])-pb0-mv(np.swapaxes(Qjb,-1,-2),m1-b['pb'][ids]);forward-=b['n'][ids]*np.sum(forward*b['n'][ids],axis=1)[:,None];back-=a['n'][ids]*np.sum(back*a['n'][ids],axis=1)[:,None];du=.5*(forward-mv(Q,back));e0=mv(Q,old['elastic'][m,ids]);et=e0+du;mag=np.linalg.norm(et,axis=1);cap=MU*b['N'][ids];rat=np.minimum(1.,np.divide(cap,self.k[ids]*mag,out=np.ones_like(mag),where=mag>0));e[m,ids]=rat[:,None]*et;slide[ids]=rat<1;dp=et-e[m,ids];tau=.5*self.k[ids,None]*e[m,ids];heat+=.5*np.sum(cap*np.linalg.norm(dp,axis=1));loss+=.25*np.sum(self.k[ids]*np.sum((e[m,ids]-e0)**2,axis=1));work+=np.sum(tau*du);armA=m1-b['pa'][ids];armB=m1-b['pb'][ids];self.point_load(g,i,self.s[ids],tau,np.cross(armA,tau),gp,gm);self.point_load(g,j,b['sb'][ids],-tau,np.cross(armB,-tau),gp,gm);surface.append(dict(i=i,j=j,s=self.s[ids],sb=b['sb'][ids],tau=tau,armA=armA,armB=armB,du=du))
  friction_q=self.pullback(g,gp,gm);Ut=.25*np.sum(self.k[None,:]*np.sum(e*e,axis=2));Uold=.25*np.sum(self.k[None,:]*np.sum(old['elastic']**2,axis=2));new=dict(q=np.array(q).copy(),g=g,elastic=e,identity=self.identity)
  return dict(state=new,residual=normal_q+friction_q,normal_residual=normal_q,friction_residual=friction_q,normal_energy_J=g['normal_energy_J'],tangential_energy_J=float(Ut),friction_heat_J=float(heat),numerical_loss_J=float(loss),tangential_endpoint_work_J=float(work),tangential_balance_J=float(work-(Ut-Uold)-heat-loss),friction_gp=gp,friction_gm=gm,surface=surface,active_set=np.array(active),sliding_set=np.array(sliding),partner_segments=np.array([s['segment'] for s in g['sides']]))
 def native(self,q,order=24):
  g=self.unpack(q);return native(g['p'],self.radius,np.tile(self.rod.ds,(self.n,1)),KC,order=order,self_contact=False)
