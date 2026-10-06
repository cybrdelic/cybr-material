"""Check only changed detail triangles and prove native top geometry stayed frozen."""
from adaptive_mesh import *
import gc
old=json.loads((ROOT/'receipts/native_coverage.json').read_text());detail=json.loads((ROOT/'receipts/molding_details.json').read_text());aa=np.load(old['arrays_path']);bb=np.load(detail['arrays_path']);q0=aa['q_m'];t0=aa['triangles'];c0=aa['chart_ids'];q=bb['q_m'];tri=bb['triangles'];chart=bb['chart_ids'];assert np.array_equal(q[:len(q0)],q0)
def keys(t):
    x=np.sort(t,axis=1).astype('i8');return (x[:,0]*len(q)+x[:,1])*len(q)+x[:,2]
k0=keys(t0);k=keys(tri);unchanged=np.isin(k,k0,assume_unique=True);changed=~unchanged;top_equal=np.array_equal(np.sort(keys(t0[c0==0])),np.sort(keys(tri[chart==0])));assert top_equal
changed_native_height_inactive=bool(np.all(chart[changed]!=0) and np.all(q[tri[changed],2]<=.001+1e-15));assert changed_native_height_inactive
new_tri=tri[changed].copy();changed_count=len(new_tri);del k0,k,q0,t0,c0,aa,bb;gc.collect();e,check=errors(q,new_tri,dense=True);assert check['passes_sampled_camera_budget']
shut=json.loads((ROOT/'receipts/shutoff_contract.json').read_text());lo,hi=shut['retained_qz_bounds_m'];side=chart>=2;z=q[tri,2];cross=side&(((z.min(1)<lo-1e-15)&(z.max(1)>lo+1e-15))|((z.min(1)<hi-1e-15)&(z.max(1)>hi+1e-15)));assert not cross.any()
result={'top_triangles_and_old_q_bitwise_unchanged':top_equal,'changed_native_height_inactive':changed_native_height_inactive,'changed_triangles':changed_count,'changed_dense_camera_audit':check,'all_other_faces_inherit_prior_lod_audit':True,'shutoff_boundary_crossings':int(cross.sum()),'new_mesh_sha256':detail['arrays_sha256'],'old_mesh_sha256':old['arrays_sha256'],'passes':True,'peak_RSS_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(ROOT/'receipts/detail_lod_verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
