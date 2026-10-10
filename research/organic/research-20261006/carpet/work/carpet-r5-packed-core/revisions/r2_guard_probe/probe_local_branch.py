"""Isolated deterministic guard correction. Never loads or edits a Blender scene."""
from pathlib import Path
import sys,json,importlib.util,hashlib,time
import numpy as np
from scipy.optimize import brentq

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
source=ROOT/'revisions/r1_frozen/source/packed_construction.py'
spec=importlib.util.spec_from_file_location('frozen_r1_model',source);model=importlib.util.module_from_spec(spec);spec.loader.exec_module(model)
original=model.trim_inner_crown
COUNTS={}

def trim_local_branch(p,axis,shear=0.,fillet_radius=35e-6):
    x=p@axis-shear*p[:,2];z=p[:,2];dx=np.diff(x)
    bad=np.flatnonzero((dx<0)&(np.arange(len(dx))>30)&(np.arange(len(dx))<len(dx)-30))
    if not len(bad):return p
    lo=int(bad[0]);hi=int(bad[-1])+1
    # R1 vetoed on any reversed segment anywhere before/after the crown.
    # The intersection needs only the maximal increasing branch touching it.
    left_bad=np.flatnonzero(np.diff(x[:lo+1])<=0)
    right_bad=np.flatnonzero(np.diff(x[hi:])<=0)
    left_start=int(left_bad[-1])+1 if len(left_bad) else 0
    right_stop=hi+int(right_bad[0])+1 if len(right_bad) else len(p)
    lx=x[left_start:lo+1];lz=z[left_start:lo+1];rx=x[hi:right_stop];rz=z[hi:right_stop]
    if len(lx)<2 or len(rx)<2:return p
    low=max(lx[0],rx[0]);high=min(lx[-1],rx[-1])
    if low>=high:return p
    fun=lambda xx:np.interp(xx,lx,lz)-np.interp(xx,rx,rz)
    if fun(low)*fun(high)>0:return p
    xx=brentq(fun,low,high,xtol=1e-14)
    il=int(np.searchsorted(lx,xx))+left_start;ir=int(np.searchsorted(rx,xx))+hi
    if il<1 or ir<1:return p
    tl=(xx-x[il-1])/(x[il]-x[il-1]);tr=(xx-x[ir-1])/(x[ir]-x[ir-1])
    cross=(p[il-1]*(1-tl)+p[il]*tl+p[ir-1]*(1-tr)+p[ir]*tr)*.5
    q=np.concatenate((p[:il],cross[None],p[ir:]))
    assert np.array_equal(q[0],p[0]) and np.array_equal(q[-1],p[-1])
    COUNTS['corrected']=COUNTS.get('corrected',0)+1
    if left_start or right_stop<len(p):COUNTS['local_bracket_used']=COUNTS.get('local_bracket_used',0)+1
    return q

def outer_envelope(points,radius,center,body_radii):
    p=points.reshape(-1,3).astype('f8');p=p[p[:,2]>.00013]
    a=center[:-1];v=np.diff(center,axis=0);vv=np.sum(v*v,axis=1)
    t=np.clip(np.sum((p[:,None]-a[None])*v[None],axis=2)/vv[None],0,1)
    d=np.linalg.norm(p[:,None]-a[None]-t[:,:,None]*v[None],axis=2).min(1)
    return float((d+float(body_radii.max())-float(np.mean(radius))).max()*1e6)

def run():
    d=np.load(ROOT/'receipts/baseline_loops.npz');reports=[];start=time.time()
    for i in [230,0,266,672,1006,400,900,1400]:
        model.trim_inner_crown=original;p0,r0=model.construct(d['center'][i],d['radius'][i],i)
        model.trim_inner_crown=trim_local_branch;p,r=model.construct(d['center'][i],d['radius'][i],i)
        endpoints_equal=bool(np.array_equal(p[:,[0,-1]],p0[:,[0,-1]]));assert endpoints_equal
        report={'loop':i,'endpoints_bitwise_equal':endpoints_equal,'fibre_count':len(p),'finite':bool(np.isfinite(p).all()),'max_control_point_change_um':float(np.linalg.norm(p-p0,axis=2).max()*1e6),'maximum_body_edge_excess_over_original_circular_envelope_um':outer_envelope(p,d['radius'][i],d['center'][i],r),'r1':model.clearance(p0,r0),'proposed_r2':model.clearance(p,r)}
        reports.append(report);print(json.dumps({'loop':i,'endpoints_equal':endpoints_equal,'r1_gap_um':report['r1']['min_gap_um'],'r2_gap_um':report['proposed_r2']['min_gap_um'],'r1_kr':report['r1']['max_curvature_times_radius'],'r2_kr':report['proposed_r2']['max_curvature_times_radius'],'envelope_excess_um':report['maximum_body_edge_excess_over_original_circular_envelope_um']}),flush=True)
    out={'scope':'Isolated 8-loop proposal only. No new scene or heavy build. r1 .blend and source snapshots unchanged.','r1_model_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'operation':'Find a locally increasing suffix/prefix adjoining the crown lobe; use it only to bracket the existing 2D intersection solve. Keep all unchanged root prefixes, suffixes, body crimp and fairing.','counts':COUNTS,'seconds':time.time()-start,'loops':reports}
    (HERE/'local_branch_probe.json').write_text(json.dumps(out,indent=2));print('Saved',HERE/'local_branch_probe.json',flush=True)
if __name__=='__main__':run()
