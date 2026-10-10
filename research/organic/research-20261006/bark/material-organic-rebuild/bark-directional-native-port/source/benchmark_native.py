"""Single-process, single-thread bounded throughput measurements, not solver timings."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
import json,resource,statistics,sys,time
import numpy as np
from native_directional import NativeEnvelope,NativePrisms
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT.parent
sys.path.insert(0,str(BASE/'bark-constitutive-incremental/source'))
from incremental_envelope import IncrementalEnvelope
from directional_candidate import synthetic_candidate
from quadratic_prism_frozen import QuadraticPrism,flat_reference

def measure(fn,count,repeats=5):
    fn();times=[]
    for _ in range(repeats):
        t=time.perf_counter();fn();times.append(time.perf_counter()-t)
    median=statistics.median(times)
    return dict(count_per_call=count,repeats=repeats,seconds=times,median_seconds=median,items_per_second=count/median)
def main():
    started=time.perf_counter();rng=np.random.default_rng(641)
    l=synthetic_candidate();reference=IncrementalEnvelope(l.tension,l.compression,l.G,l.axes,l.minimum_stretch,l.maximum_stretch);native=NativeEnvelope(l)
    n=65536;h=rng.normal(size=(n,3,3))*.02
    point_native=measure(lambda:native.evaluate_incremental(h),n)
    point_numpy=measure(lambda:reference.evaluate_incremental(h),n)
    scalar_count=256
    point_scalar=measure(lambda:[reference.evaluate_incremental(v) for v in h[:scalar_count]],scalar_count,3)
    X=flat_reference(np.array([[0.,0.],[.001,0.],[0.,.001]]),.001);cell=QuadraticPrism(X,1e6,0.,240.,nq=4)
    count=128;q=len(cell.weights);gradients=np.repeat(cell.grad[None],count,axis=0);weights=np.repeat(cell.weights[None],count,axis=0)
    h0=np.broadcast_to(np.diag([-.1,.02,0.]),(count,q,3,3)).copy();u=rng.normal(size=(count,18,3))*1e-6
    prism=NativePrisms(gradients,weights,h0,np.eye(3),240.)
    def numpy_prisms():
        hh=h0+np.einsum('cia,cqib->cqab',u-u.mean(axis=1,keepdims=True),gradients)
        w,p,d=reference.evaluate_incremental(hh)
        g=np.einsum('cqib,cqab,cq->cia',gradients,p,weights);g-=g.mean(axis=1,keepdims=True)
        return np.sum(weights*w),g,np.sum(weights*d['J'],axis=1)
    en,gn,dn=prism.evaluate(u);er,gr,vr=numpy_prisms()
    assert abs(en-er)/er<1e-12 and np.linalg.norm(gn-gr)/np.linalg.norm(gr)<1e-12
    prism_native=measure(lambda:prism.evaluate(u),count)
    prism_numpy=measure(numpy_prisms,count)
    cpu='unknown'
    for line in Path('/proc/cpuinfo').read_text().splitlines():
        if line.startswith('model name'):cpu=line.split(':',1)[1].strip();break
    receipt=dict(scope='Steady-state end-to-end Python wrapper calls, outputs and input validation included. Native CPU single-threaded; no OpenMP and no solve.',CPU=cpu,thread_environment={k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS')},material_points={'native':point_native,'frozen_numpy_vectorized':point_numpy,'frozen_numpy_scalar_loop':point_scalar,'speedup_vs_vectorized':point_numpy['median_seconds']/point_native['median_seconds'],'throughput_ratio_vs_scalar_loop':point_native['items_per_second']/point_scalar['items_per_second']},prisms={'cells':count,'quadrature_points_per_cell':q,'native':prism_native,'frozen_numpy_vectorized':prism_numpy,'native_quadrature_points_per_second':q*prism_native['items_per_second'],'speedup_vs_vectorized':prism_numpy['median_seconds']/prism_native['median_seconds']},seconds=time.perf_counter()-started,maximum_resident_memory_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (ROOT/'receipts/throughput.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
