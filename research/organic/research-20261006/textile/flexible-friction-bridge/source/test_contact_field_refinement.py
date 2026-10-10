from pathlib import Path
import argparse,json,time,resource,numpy as np
from flexible_contact import FlexibleContact
from rod_pullback import exp,log,right_jacobian,mv,L,S
P=Path(__file__).resolve().parents[1]
def exact(s):return exp(np.c_[.18*np.sin(4*np.pi*s/L),.12*np.cos(4*np.pi*s/L),.04*np.sin(2*np.pi*s/L)])
def run(M,n=32):
 start=time.monotonic();f=FlexibleContact(segments=M,order=33);arc=(f.rod.s[:-1]+f.rod.s[1:])/2;R0=exact(arc);angle=np.arange(3)*2*np.pi/3;o0=np.c_[1.99*f.radius/np.sqrt(3)*np.cos(angle),1.99*f.radius/np.sqrt(3)*np.sin(angle),np.zeros(3)]
 def controls(t):
  wave=np.sin(2*np.pi*t);phi=np.array([.4,-.2,.1])*wave;z=1e-6*wave;R=R0[None]@exp(np.c_[np.zeros(3),np.zeros(3),phi])[:,None];o=o0.copy();o[0,2]+=z;return f.pack(o,R),phi,z
 q,phio,zo=controls(0);old=f.state(q);g=old['g'];D=f.surface(g,0,f.s)[0];_,_,j,_=f.sample(g,0,f.s);deviation=np.linalg.norm(D[:,:,2]-g['R'][0,j,:,2],axis=1);reference=exact(f.s);frameerr=np.linalg.norm(log(np.swapaxes(D,-1,-2)@reference),axis=1);rows=[];totalW=totalD=totalA=totalG=0.;Eold=0.;status='running'
 for k in range(1,n+1):
  if time.monotonic()-start>110:status='bounded_pause';break
  q,phi,z=controls(k/n);t=f.trial(q,old);g=t['state']['g'];r=t['friction_residual'].reshape(3,-1);body=mv(right_jacobian(g['w'],inverse=True),np.broadcast_to([0.,0.,1.],g['w'].shape));reaction_phi=np.sum(r[:,3:].reshape(3,M,3)*body,axis=(1,2));Fz=r[0,2]/S;W=float(Fz*(z-zo)+reaction_phi@(phi-phio));Dheat=t['friction_heat_J'];A=t['numerical_loss_J'];Et=t['tangential_energy_J'];gap=W-(Et-Eold)-Dheat-A;totalW+=W;totalD+=Dheat;totalA+=A;totalG+=gap;rows.append(dict(step=k,z_m=z,phi_rad=phi.tolist(),reaction_phi_Nm=reaction_phi.tolist(),reaction_z_N=float(Fz),friction_work_J=W,heat_J=Dheat,numerical_loss_J=A,stored_J=Et,geometric_work_gap_J=gap,active=int(t['active_set'].sum()),sliding=int(t['sliding_set'].sum())));old=t['state'];phio=phi;zo=z;Eold=Et
 if status=='running':status='completed'
 result=dict(scope='Prescribed roll/slide history, contact-field discretization only; not equilibrium or manufacturing.',segments=M,steps=n,status=status,max_director_tangent_deviation_rad=float(deviation.max()),RMS_director_tangent_deviation_rad=float(np.sqrt(np.sum(f.w*deviation**2)/L)),max_exact_frame_error_rad=float(frameerr.max()),friction_work_J=totalW,heat_J=totalD,numerical_loss_J=totalA,stored_J=Eold,ledger_gap_J=totalG,ledger_relative=abs(totalG)/max(sum(abs(a['friction_work_J']) for a in rows),1e-18),rows=rows,wall_s=time.monotonic()-start,max_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024);(P/'receipts'/f'contact_field_M{M}.json').write_text(json.dumps(result,indent=2));print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--segments',type=int,required=True);ap.add_argument('--steps',type=int,default=32);a=ap.parse_args();run(a.segments,a.steps)
