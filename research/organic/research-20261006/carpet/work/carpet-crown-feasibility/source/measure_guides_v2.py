"""Orthogonal rounded cap proposal; no source yarn or fibres are changed."""
import json
import numpy as np
import measure_guides as v1
from scipy.interpolate import CubicSpline
from scipy.spatial import cKDTree

def guide(c,r):
    midpoint=(c[0]+c[-1])/2
    axis=c[-1]-c[0];span=np.linalg.norm(axis);axis/=span
    rise=c[24]-midpoint;lean=float(rise@axis)
    upright=rise-lean*axis;height=np.linalg.norm(upright);upright/=height
    crown=1.25*r;leg_height=height-crown
    assert crown<span/2 and leg_height>0
    q=np.linspace(0,1,513);ease=q*q*(3-2*q)
    left_x=-span/2+(lean-crown+span/2)*ease
    right_x=span/2+(lean+crown-span/2)*ease
    leg_z=leg_height*q
    theta=np.linspace(-np.pi/2,np.pi/2,1025)
    cap_x=lean+crown*np.sin(theta);cap_z=height-crown+crown*np.cos(theta)
    xx=np.r_[left_x[:-1],cap_x,right_x[-2::-1]]
    zz=np.r_[leg_z[:-1],cap_z,leg_z[-2::-1]]
    p=midpoint+xx[:,None]*axis+zz[:,None]*upright
    p[0]=c[0];p[-1]=c[-1];p[1024]=c[24]
    return p,{'crown_radius_m':float(crown),'root_span_m':float(span),'height_in_loop_plane_m':float(height),'pose_scope':'Exact root centres and apex, hence root-axis yaw and total apex lean/sway. Shoulders change. Circular cap is placed in an orthogonal plane; it is not sheared by a height-varying lean.'}

def run():
    d=np.load(v1.DATA);rows=[]
    for index in (968,955,356):
        c=d['center'][index];r=float(d['radius'][index].mean())
        old=CubicSpline(np.linspace(0,1,len(c)),c,axis=0)(np.linspace(0,1,4097))
        proposed,details=guide(c,r)
        a=cKDTree(old).query(proposed)[0];b=cKDTree(proposed).query(old)[0]
        k=v1.curvature(proposed)
        rows.append({'loop':index,**details,'bundle_radius_m':r,'root_centres_exact':bool(np.array_equal(proposed[[0,-1]],c[[0,-1]])),'apex_exact':bool(np.array_equal(proposed[1024],c[24])),'dense_sample_symmetric_guide_distance_m':float(max(a.max(),b.max())),'maximum_sampled_guide_curvature_times_bundle_radius':float(k.max()*r),'axis_aligned_centerline_extent_change_m':((proposed.max(0)-proposed.min(0))-(old.max(0)-old.min(0))).tolist()})
    report={'scope':'Cheap three-loop guide proposal only; no fibres, wrappers, scene or render.','assumption':'Cap radius1.25 times bundle radius is a geometric margin, not wool calibration.','supersedes_geometry_only_v1':'The height-varying shear in v1 raised guide curvature above1/bundle radius. This version places an actual circle in the plane defined by roots and apex.','limits':['Dense-sample guide distance is not an exact Hausdorff certificate.','Sampled guide curvature is not a fibre/contact proof.','Exact old wrapper coordinates and2micrometre envelope are not preserved.','No inter-yarn or fibre-packing acceptance.'],'rows':rows}
    (v1.ROOT/'receipts/proposed_guide_tradeoff_v2.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':run()
