"""Read-only section trace of one already rejected frozen unit; no projection."""
import json,time,resource
import numpy as np
from scipy.interpolate import CubicSpline
from construct_unit import ROOT,R1,migration,unit

def run():
    start=time.time();index=968;data=np.load(R1/'receipts/baseline_loops.npz');saved=np.load(ROOT/f'arrays/rounded_unit_{index}.npz')
    p=migration.evaluate_catmull(saved['body'],8).astype('f8');guide=saved['guide'];r=float(saved['body_radius']);R=float(data['radius'][index].mean())
    arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(guide,axis=0),axis=1))];spline=CubicSpline(arc,guide,axis=0)
    s=np.linspace(0,arc[-1],193);c=spline(s);velocity=spline(s,1);t=unit(velocity);speed=np.linalg.norm(velocity,axis=1);accel=spline(s,2)
    bend=(accel-t*np.sum(accel*t,axis=1)[:,None])/speed[:,None]**2
    old=data['center'][index];axis=unit(old[-1]-old[0]);normal=unit(np.cross(axis,old[24]-(old[0]+old[-1])*.5));in_plane=unit(np.cross(normal[None],t));a0=data['frame'][index,0,0];alpha=np.arctan2(a0@normal,a0@in_plane[0]);a=in_plane*np.cos(alpha)+normal*np.sin(alpha);b=unit(np.cross(t,a))
    rate=float((data['phase'][index,-1]-data['phase'][index,0])/arc[-1]);phase=float(data['phase'][index,0])+rate*s
    e1=np.cos(phase)[:,None]*a+np.sin(phase)[:,None]*b;e2=-np.sin(phase)[:,None]*a+np.cos(phase)[:,None]*b
    curvature=np.stack((np.sum(bend*e1,axis=1),np.sum(bend*e2,axis=1)),axis=1)
    first,second=np.triu_indices(187,1);rows=[];coordinates=[]
    for k in range(len(s)):
        signed=np.sum((p-c[k])*t[k],axis=2);cross=(signed[:,:-1]*signed[:,1:]<=0)&(signed[:,:-1]!=signed[:,1:]);fi,segment=np.nonzero(cross)
        fraction=signed[fi,segment]/(signed[fi,segment]-signed[fi,segment+1]);hit=p[fi,segment]+fraction[:,None]*(p[fi,segment+1]-p[fi,segment]);distance=np.linalg.norm(hit-c[k],axis=1)
        nearby=distance<=R+20e-6;counts=np.bincount(fi[nearby],minlength=187);selected=np.full((187,3),np.nan)
        order=np.argsort(distance)
        for j in order[::-1]:
            if nearby[j]:selected[fi[j]]=hit[j]
        xy=np.stack(((selected-c[k])@e1[k],(selected-c[k])@e2[k]),axis=1);coordinates.append(xy)
        valid=np.isfinite(xy[first]).all(1)&np.isfinite(xy[second]).all(1);i=first[valid];j=second[valid];d=xy[i]-xy[j];mid=(xy[i]+xy[j])*.5;v=rate*np.stack((-mid[:,1],mid[:,0]),axis=1);axial=1-mid@curvature[k]
        metric=d-v*(np.sum(v*d,axis=1)/(axial*axial+np.sum(v*v,axis=1)))[:,None];gap=np.sqrt(np.maximum(np.sum(d*metric,axis=1),0))-2*r
        worst=int(np.argmin(gap));rows.append({'section':k/2,'location':'original_section_station' if k%2==0 else 'halfway_station','guide_arclength_m':float(s[k]),'guide_z_m':float(c[k,2]),'fibres_with_no_nearby_plane_crossing':int((counts==0).sum()),'fibres_with_multiple_nearby_crossings':int((counts>1).sum()),'minimum_declared_metric_physical_gap_um':float(gap.min()*1e6),'pairs_below_physical_diameter':int((gap<0).sum()),'pairs_below_declared_2um_margin':int((gap<2e-6).sum()),'worst_pair':[int(i[worst]),int(j[worst])]})
    summaries={}
    for label in ('original_section_station','halfway_station'):
        region=[x for x in rows if x['location']==label and x['guide_z_m']>.0002]
        summaries[label]={'stations':len(region),'stations_with_metric_physical_overlap':sum(x['pairs_below_physical_diameter']>0 for x in region),'minimum_metric_gap_um':min(x['minimum_declared_metric_physical_gap_um'] for x in region),'maximum_overlapping_pairs':max(x['pairs_below_physical_diameter'] for x in region),'missing_crossings':sum(x['fibres_with_no_nearby_plane_crossing'] for x in region),'multiple_crossings':sum(x['fibres_with_multiple_nearby_crossings'] for x in region)}
    xy=np.asarray(coordinates);step=np.linalg.norm(np.diff(xy[::2],axis=0),axis=2)
    result={'loop':index,'scope':'Read-only intersections of frozen final native curves with the original 97 guide section planes and their 96 midplanes. No seed, projection or correction rerun.','guide_section_spacing_um':float(arc[-1]/96*1e6),'curve_evaluation_subdivisions':8,'summary':summaries,'maximum_final_transverse_motion_between_original_stations_um':float(np.nanmax(step)*1e6),'median_final_transverse_motion_between_original_stations_um':float(np.nanmedian(step)*1e6),'interpretation':'Failures at original section stations in the final curves disprove a claim that this final field is only failing between otherwise clear stations. Original pre-interpolation XY sections were not saved, so this trace cannot determine whether initial section projection, interpolation or later contact correction first introduced each overlap. The projection routine has no post-clamp separation assertion.','limits':['The reported metric is the same curvature/twist approximation used during construction. Independent migration derivatives are not included in that metric.','Plane intersections approximate the native cubic with eight subdivisions. Multiple nearby plane crossings are separately counted.','No fibre-identity permutation is inferred merely from contact; curves preserve their original array identities.'],'rows':rows,'seconds':time.time()-start,'peak_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024}
    (ROOT/'receipts/frozen_section_trace_968.json').write_text(json.dumps(result,indent=2));print(json.dumps({k:v for k,v in result.items() if k!='rows'}))

if __name__=='__main__':run()
