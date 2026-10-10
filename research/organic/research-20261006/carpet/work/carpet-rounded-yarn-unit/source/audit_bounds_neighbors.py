"""Read-only diagnostics of the already rejected unit arrays; no construction."""
import json
import numpy as np
from scipy.spatial import cKDTree
from construct_unit import ROOT,R1,migration
from probe_unit import old_wraps

def resample(p,count=2049):
    arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(p,axis=0),axis=1))]
    s=np.linspace(0,arc[-1],count)
    return np.stack([np.interp(s,arc,p[:,k]) for k in range(3)],axis=1)

def core_envelope(center,R,data,index,neighbors):
    c=data['center'][neighbors];a=c[:,:-1].reshape(-1,3);b=c[:,1:].reshape(-1,3)
    radius=((data['radius'][neighbors,:-1]+data['radius'][neighbors,1:])*.5).ravel()
    p=resample(center);pa=p[:-1];pb=p[1:];mid=(pa+pb)*.5;visible=mid[:,2]>.00013
    ids=cKDTree((a+b)*.5).query(mid,k=48,workers=1)[1]
    gaps=(migration.old.segment_distance(pa[:,None],pb[:,None],a[ids],b[ids])-R-radius[ids]).min(1)
    return {'exposed_centerline_arclength_fraction_with_envelope_overlap':float(np.mean(gaps[visible]<0)),'minimum_exposed_envelope_gap_um':float(gaps[visible].min()*1e6),'uniform_arclength_segments':len(pa)}

def bounds(paths,radii,visible):
    mins=[];maxs=[]
    for p,r in zip(paths,radii):
        p=p.reshape(-1,3);p=p[p[:,2]>.00013] if visible else p
        mins.append((p-r).min(0));maxs.append((p+r).max(0))
    return np.min(mins,0),np.max(maxs,0)

def run():
    data=np.load(R1/'receipts/baseline_loops.npz');report=json.loads((ROOT/'receipts/rounded_unit_probe.json').read_text());rows=[]
    for entry in report['rows']:
        index=entry['loop'];saved=np.load(ROOT/f'arrays/rounded_unit_{index}.npz');r=float(saved['body_radius']);R=float(data['radius'][index].mean());neighbors=entry['neighboring_corrected']['neighbor_indices']
        body=migration.evaluate_catmull(saved['body'],6);wraps=migration.evaluate_catmull(saved['wraps'],6)
        color=int(data['color'][index]);wi=int(np.count_nonzero(data['color'][:index]==color));old=np.load(R1/f'arrays/colour_{color}_positions.npy',mmap_mode='r')[wi];old=migration.evaluate_catmull(old,6);ow=old_wraps(data,index)
        old_visible=bounds([old,ow],[r,16e-6],True);new_visible=bounds([body,wraps],[r,16e-6],True)
        angles=np.arange(6)*np.pi/3;core=data['center'][index,:,None]+data['radius'][index,:,None,None]*(np.cos(angles)[None,:,None]*data['frame'][index,:,0,None]+np.sin(angles)[None,:,None]*data['frame'][index,:,1,None])
        r5_visible=bounds([core,ow],[0.,16e-6],True)
        endpoints=np.concatenate((saved['body'][:,[0,-1]],saved['wraps'][:,[0,-1]]))
        rows.append({'loop':index,'original_r5_core_envelope':core_envelope(data['center'][index],R,data,index,neighbors),'rounded_core_envelope':core_envelope(saved['guide'],R,data,index,neighbors),'r1_visible_bounds_m':{'minimum':old_visible[0].tolist(),'maximum':old_visible[1].tolist()},'r5_visible_bounds_m':{'minimum':r5_visible[0].tolist(),'maximum':r5_visible[1].tolist()},'new_visible_bounds_m':{'minimum':new_visible[0].tolist(),'maximum':new_visible[1].tolist()},'visible_extent_delta_from_r1_um':((new_visible[1]-new_visible[0])-(old_visible[1]-old_visible[0])).tolist() if False else (((new_visible[1]-new_visible[0])-(old_visible[1]-old_visible[0]))*1e6).tolist(),'visible_extent_delta_from_r5_um':(((new_visible[1]-new_visible[0])-(r5_visible[1]-r5_visible[0]))*1e6).tolist(),'buried_endpoint_bounds_m':{'minimum':endpoints.min((0,1)).tolist(),'maximum':endpoints.max((0,1)).tolist()},'full_bounds_include_long_buried_tangent_extensions':True})
    result={'scope':'Read-only interpretation of the existing rejected three-loop outputs. No new construction or scene.','caveats':['The original_r5_core_envelope comparison uses circular capsules enclosing the original six-sided core. Negative gaps establish envelope overlap, not exact mesh intersection.','Rounded and original core-envelope fractions use uniform centerline arclength; they cannot be compared numerically to the body-fibre segment fractions in the main receipt.','R5 visible bounds use recovered six-sided core vertices and 16 micrometre circumscribed wrapper radii. Candidate and r1 bounds use densely sampled Catmull centerlines plus circular radii.','Full bounds include buried tangent extensions; visible bounds filter centers above the 0.13 millimetre backing threshold.','No baseline layout is certified as nonintersecting. Internal body clearance already rejects every unit independently of the neighbor screen.'],'rows':rows}
    (ROOT/'receipts/bounds_neighbor_interpretation.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))

if __name__=='__main__':run()
