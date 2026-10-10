"""Cheap guide-shape feasibility proposal, not a fibre packing or render pass."""
from pathlib import Path
import json
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.spatial import cKDTree

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT.parent/'carpet-r5-packed-core/receipts/baseline_loops.npz'

def guide(c,r):
    midpoint=(c[0]+c[-1])/2
    axis=c[-1]-c[0];axis[2]=0;span=np.linalg.norm(axis);axis/=span
    height=c[24,2]-midpoint[2]
    drift=c[24]-midpoint;drift[2]=0
    crown=1.25*r;leg_height=height-crown
    assert crown<span/2 and leg_height>0
    q=np.linspace(0,1,513);ease=q*q*(3-2*q)
    x=-span/2+(span/2-crown)*ease;z=leg_height*q
    theta=np.linspace(-np.pi/2,np.pi/2,1025)
    ax=crown*np.sin(theta);az=height-crown+crown*np.cos(theta)
    xx=np.r_[x[:-1],ax,-x[-2::-1]];zz=np.r_[z[:-1],az,z[-2::-1]]
    h=zz/height;d=h*h*(3-2*h)
    p=midpoint+xx[:,None]*axis+zz[:,None]*np.array([0,0,1])+d[:,None]*drift
    p[0]=c[0];p[-1]=c[-1];p[1024]=c[24]
    return p,{'crown_radius_m':float(crown),'root_span_m':float(span),'height_m':float(height),'pose_scope':'Exact root centres and apex; same root-axis yaw and apex lean/sway. Intermediate shoulders/drift are changed.'}

def curvature(p):
    s=np.r_[0,np.cumsum(np.linalg.norm(np.diff(p,axis=0),axis=1))]
    dp=np.gradient(p,s,axis=0);ddp=np.gradient(dp,s,axis=0)
    return np.linalg.norm(np.cross(dp,ddp),axis=1)/np.maximum(np.linalg.norm(dp,axis=1)**3,1e-30)

def run():
    d=np.load(DATA);rows=[]
    for index in (968,955,356):
        c=d['center'][index];r=float(d['radius'][index].mean())
        old=CubicSpline(np.linspace(0,1,len(c)),c,axis=0)(np.linspace(0,1,4097))
        proposed,details=guide(c,r)
        a=cKDTree(old).query(proposed)[0];b=cKDTree(proposed).query(old)[0]
        k=curvature(proposed)
        rows.append({'loop':index,**details,'bundle_radius_m':r,'root_centres_exact':bool(np.array_equal(proposed[[0,-1]],c[[0,-1]])),'apex_exact':bool(np.array_equal(proposed[1024],c[24])),'dense_sample_symmetric_guide_distance_m':float(max(a.max(),b.max())),'maximum_sampled_guide_curvature_times_bundle_radius':float(k.max()*r),'axis_aligned_centerline_extent_change_m':((proposed.max(0)-proposed.min(0))-(old.max(0)-old.min(0))).tolist(),'change_is_not_within_previous_2um_envelope_gate':True})
    result={'scope':'Three-loop numerical scope proposal only. No fibres, wrappers, scene, source replacement or render.','assumption':'Rounded crown radius1.25 times nominal yarn radius; geometric clearance margin, not a wool material calibration.','limits':['Dense-sample guide distance is not an exact Hausdorff certificate.','Curvature is a sampled guide estimate, not a full fibre bound.','The previous frozen wrapper coordinates and2um envelope cannot be claimed preserved.','No packing/contact or inter-yarn clearance is established.'],'rows':rows}
    (ROOT/'receipts/proposed_guide_tradeoff.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':run()
