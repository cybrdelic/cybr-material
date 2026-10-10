"""Bounded generalized Hessian estimate of the last accepted local equilibrium."""
import argparse,json,time,resource
import numpy as np
from scipy.sparse.linalg import LinearOperator,eigsh,ArpackNoConvergence
from opening_model import OpeningModel,History,ROOT
from checkpoint_io import read_checkpoint

def main(label):
    start=time.monotonic();dest=ROOT/'state'/label;a,meta=read_checkpoint(dest/'restart.npz');p=meta['parameters']
    m=OpeningModel(p['radial_pitch'],p['quadrature'],p['penalty'],p['Gc']);u=a['displacement'];m.set_load(meta['trace'][-1]['load'])
    if (dest/'previous_accepted.npz').exists():
        old,_=read_checkpoint(dest/'previous_accepted.npz');history=History.make(old['maximum'],old['damage']);origin='previous accepted station'
    else:history=m.initial_history;origin='prepared initial history (pilot has no new damage)'
    E,g,h,r,o=m.evaluate(u,history);factor=m.precondition(h.damage);n=len(m.free);calls=0
    def matvec(v):
        nonlocal calls
        if time.monotonic()-start>85:raise TimeoutError('Stability operator budget reached')
        length=np.linalg.norm(v)
        if length==0:return np.zeros_like(v,dtype=np.float64)
        eps=1e-8/length;du=np.zeros(u.size);du[m.free]=v;du=du.reshape(u.shape)
        gp=m.evaluate(u+eps*du,history)[1].ravel()[m.free];gm=m.evaluate(u-eps*du,history)[1].ravel()[m.free]
        calls+=2
        return (gp-gm)/(2*eps)
    A=LinearOperator((n,n),matvec=matvec,dtype=np.float64);Minv=LinearOperator((n,n),matvec=factor.solve,dtype=np.float64)
    rng=np.random.default_rng(9641);v=rng.normal(size=n);w=rng.normal(size=n);Hv=matvec(v);Hw=matvec(w)
    symmetry=abs(v@Hw-w@Hv)/max(abs(v@Hw),abs(w@Hv),1e-12)
    result=dict(label=label,history_origin=origin,load=m.load,finite_difference_global_displacement_norm_m=1e-8,
                symmetry_bilinear_relative_error=symmetry,operator='Hessian of U_bulk + U_interface + new Phi, generalized by the positive secant preconditioner',
                scope='Numerical local discrete stability estimate at this station, not a global minimum or material ellipticity proof')
    try:
        values,vectors=eigsh(A,k=1,M=m.precondition_matrix,Minv=Minv,which='SA',tol=3e-3,maxiter=80,ncv=16)
        v=vectors[:,0];Hv=matvec(v);Mv=m.precondition_matrix@v
        residual=float(np.linalg.norm(Hv-values[0]*Mv)/max(np.linalg.norm(Hv),np.linalg.norm(values[0]*Mv),1e-30))
        result.update(status='PASS' if values[0]>0 and residual<.02 else 'UNRESOLVED',smallest_generalized_eigenvalue=float(values[0]),eigenpair_relative_residual=residual)
    except (ArpackNoConvergence,TimeoutError,ValueError) as error:
        result.update(status='UNRESOLVED',error=str(error))
    result.update(force_evaluations=calls,wall_s=time.monotonic()-start,peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
    (ROOT/'receipts'/f'{label}_stability.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('label');main(p.parse_args().label)
