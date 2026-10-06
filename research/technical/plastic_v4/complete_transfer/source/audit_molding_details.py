"""Bounded dedicated preservation audit of retained ejectors and parting witness."""
from adaptive_mesh import *
from scipy.spatial import cKDTree
import gc
DETAIL_BUDGET_M=.5e-6

def locate(q,tri,chart,query,which):
    ids=np.flatnonzero(chart==which);t=tri[ids];uvaxes={1:[1,0],2:[1,2],3:[2,1],4:[0,2],5:[2,0]}[which];uv=q[t][:,:,uvaxes];tree=cKDTree(uv.mean(1));hit=np.empty(len(query),'i4');weights=np.empty((len(query),3))
    for start in range(0,len(query),128):
        x=query[start:start+128][:,uvaxes];remaining=np.arange(len(x));picked=np.empty(len(x),'i4');bary=np.empty((len(x),3))
        for k in [64,256,1024,4096]:
            if not len(remaining):break
            _,near=tree.query(x[remaining],k=min(k,len(ids)));corners=uv[near];r=x[remaining,None,:]-corners[:,:,0];ab=corners[:,:,1]-corners[:,:,0];ac=corners[:,:,2]-corners[:,:,0];det=ab[:,:,0]*ac[:,:,1]-ab[:,:,1]*ac[:,:,0];b=(r[:,:,0]*ac[:,:,1]-r[:,:,1]*ac[:,:,0])/det;c=(ab[:,:,0]*r[:,:,1]-ab[:,:,1]*r[:,:,0])/det;w=np.stack([1-b-c,b,c],-1);quality=w.min(-1);best=quality.argmax(1);inside=quality[np.arange(len(remaining)),best]>=-1e-7;solved=remaining[inside];picked[solved]=ids[near[np.flatnonzero(inside),best[inside]]];bary[solved]=w[np.flatnonzero(inside),best[inside]];remaining=remaining[~inside]
        assert not len(remaining), 'Bounded point-location failed; no geometry conclusion'
        hit[start:start+len(x)]=picked;weights[start:start+len(x)]=bary
    return hit,weights

def queries():
    out=[];angles=np.arange(512)/512*2*np.pi;radii=np.r_[np.linspace(0,.0018,17),np.linspace(.0018,.0026,81)];A,D=np.meshgrid(angles,radii)
    for cx in [-.023,.023]:
        for cy in [-.023,.023]:out.append((f'ejector_{cx}_{cy}',1,np.column_stack([cx+D.ravel()*np.cos(A.ravel()),cy+D.ravel()*np.sin(A.ravel()),np.full(A.size,-.004)])))
    u=np.unique(np.r_[np.linspace(-.04,.04,257),np.linspace(-.04,-.037,65),np.linspace(.037,.04,65)]);z=np.linspace(-.00018,.00018,181);U,Z=np.meshgrid(u,z)
    for c,name in [(2,'right'),(3,'left'),(4,'front'),(5,'back')]:
        q=np.empty((U.size,3));q[:,2]=Z.ravel()
        if c in [2,3]:q[:,0]=.04 if c==2 else -.04;q[:,1]=U.ravel()
        else:q[:,1]=-.04 if c==4 else .04;q[:,0]=U.ravel()
        out.append(('parting_'+name,c,q))
    return out

if __name__=='__main__':
    start=time.monotonic();prior=json.loads((ROOT/'receipts/native_coverage.json').read_text());a=np.load(prior['arrays_path']);q=a['q_m'];tri=a['triangles'];chart=a['chart_ids'];history=[]
    for iteration in range(8):
        p=surface(q).astype('f4').astype('f8');mark=np.zeros(len(tri),bool);rows=[]
        for name,ch,points in queries():
            hit,w=locate(q,tri,chart,points,ch);actual=np.einsum('tij,ti->tj',p[tri[hit]],w);expected=surface(points);e=np.linalg.norm(actual-expected,axis=-1);np.logical_or.at(mark,hit,e>DETAIL_BUDGET_M*.7);rows.append({'feature':name,'samples':len(points),'max_error_m':float(e.max()),'rms_error_m':float(np.sqrt(np.mean(e*e)))})
        history.append({'iteration':iteration,'triangles':len(tri),'rows':rows});print(json.dumps(history[-1]),flush=True)
        if max(r['max_error_m'] for r in rows)<=DETAIL_BUDGET_M:break
        q,tri,chart=refine(q,tri,chart,mark);assert len(tri)<MAX_TRIANGLES
    result={'history':history,'detail_budget_m':DETAIL_BUDGET_M,'passes_sampled_details':max(r['max_error_m'] for r in rows)<=DETAIL_BUDGET_M,'triangles':len(tri),'vertices':len(q),'topology':topology(q,tri),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'peak_RSS_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'seconds':time.monotonic()-start,'sampled_not_supremum':True,'no_parameter_change':True}
    if result['passes_sampled_details']:
        path=ROOT/'receipts/camera_lod_mesh_detail_preserved.npz';np.savez_compressed(path,q_m=q,position_m=surface(q).astype('f4'),triangles=tri,chart_ids=chart);result.update(arrays_path=str(path),arrays_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),arrays_bytes=path.stat().st_size)
    (ROOT/'receipts/molding_details.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
