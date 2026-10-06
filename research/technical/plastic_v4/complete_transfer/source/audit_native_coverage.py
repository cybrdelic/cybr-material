"""Independent native-pitch q sweep, stronger than training barycentric probes.
No claim that q lattice points equal physical XY texel knots on curved charts.
"""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/plastic-mpl-cache')
from adaptive_mesh import *
import matplotlib.tri as mtri

def audit(q,tri,chart):
    p=surface(q).astype('f4').astype('f8');err=np.zeros(len(tri));rows=[]
    # 4096-pitch q sweep for the top; original witness/ejector scales independently.
    for ch in range(6):
        ids=np.flatnonzero(chart==ch);t=tri[ids];axes={0:(0,1),1:(1,0),2:(1,2),3:(2,1),4:(0,2),5:(2,0)}[ch];u,v=axes
        mesh=mtri.Triangulation(q[:,u],q[:,v],t);finder=mesh.get_trifinder();maxabs=0.;counters={k:{'max_absolute_m':0.,'max_projected_pixels':0.,'samples':0} for k in CAMERAS}
        if ch==0:au=(np.arange(4096)+.5)/4096*.08-.04;av=au
        elif ch==1:
            # Broad backside plus explicit concentric ejector annuli.
            au=np.linspace(-.04,.04,1025);av=au
        elif u==2:au=np.unique(np.r_[np.linspace(-.004,.004,513),np.linspace(-.0002,.0002,161)]);av=(np.arange(2048)+.5)/2048*.08-.04
        else:au=(np.arange(2048)+.5)/2048*.08-.04;av=np.unique(np.r_[np.linspace(-.004,.004,513),np.linspace(-.0002,.0002,161)])
        for start in range(0,len(av),8):
            X,Y=np.meshgrid(au,av[start:start+8]);x=X.ravel();y=Y.ravel();idx=finder(x,y);assert np.all(idx>=0);tt=t[idx];A=q[tt[:,0]][:,[u,v]];B=q[tt[:,1]][:,[u,v]];C=q[tt[:,2]][:,[u,v]];xy=np.column_stack([x,y]);ab=B-A;ac=C-A;r=xy-A;det=ab[:,0]*ac[:,1]-ab[:,1]*ac[:,0];wb=(r[:,0]*ac[:,1]-r[:,1]*ac[:,0])/det;wc=(ab[:,0]*r[:,1]-ab[:,1]*r[:,0])/det;w=np.column_stack([1-wb-wc,wb,wc]);query=np.einsum('tij,ti->tj',q[tt],w);exact=surface(query);actual=np.einsum('tij,ti->tj',p[tt],w);delta=actual-exact;absolute=np.linalg.norm(delta,axis=-1);maxabs=max(maxabs,float(absolute.max()))
            for name,c in CAMERAS.items():
                m=np.array(c['matrix']).reshape(4,4);field=c['data']['ortho_scale'];xy=(exact-m[:3,3])@m[:3,:2];visible=np.max(np.abs(xy),axis=-1)<field/2+.00075;bound=absolute/field*RESOLUTION;projected=np.linalg.norm(delta@m[:3,:2],axis=-1)/field*RESOLUTION;np.maximum.at(err,ids[idx],np.where(visible,bound,0))
                if visible.any():
                    rr=counters[name];rr['max_absolute_m']=max(rr['max_absolute_m'],float(absolute[visible].max()));rr['max_projected_pixels']=max(rr['max_projected_pixels'],float(projected[visible].max()));rr['samples']+=int(visible.sum())
        row={'chart':ch,'lattice_samples':len(au)*len(av),'all_samples_max_absolute_m':maxabs,'cameras':counters};rows.append(row);print(json.dumps(row),flush=True);del finder,mesh
    return err,rows

if __name__=='__main__':
    start=time.monotonic();src=ROOT/'receipts/camera_lod_mesh.npz';a=np.load(src);q=a['q_m'];tri=a['triangles'];chart=a['chart_ids'];history=[]
    for iteration in range(5):
        e,rows=audit(q,tri,chart);history.append({'iteration':iteration,'vertices':len(q),'triangles':len(tri),'max_conservative_pixels':float(e.max()),'rows':rows});print(json.dumps({'iteration':iteration,'max_pixels':float(e.max()),'triangles':len(tri)}),flush=True)
        if e.max()<=POSITION_PIXELS:break
        q,tri,chart=refine(q,tri,chart,e>REFINE_PIXELS);assert len(tri)<=MAX_TRIANGLES
    dense,check=errors(q,tri,dense=True);topo=topology(q,tri);passes=e.max()<=POSITION_PIXELS and check['passes_sampled_camera_budget']
    receipt={'lattice_history':history,'final_dense':check,'topology':topo,'passes_numerical_lod':bool(passes),'vertices':len(q),'triangles':len(tri),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'adaptive_source_sha256':hashlib.sha256((ROOT/'source/adaptive_mesh.py').read_bytes()).hexdigest(),'surface_contract_sha256':hashlib.sha256((ROOT/'source/surface_contract.py').read_bytes()).hexdigest(),'peak_RSS_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'seconds':time.monotonic()-start,'not_global20nm_qualification':True,'not_supremum_proof':True}
    if passes:
        path=ROOT/'receipts/camera_lod_mesh_final.npz';np.savez_compressed(path,q_m=q,position_m=surface(q).astype('f4'),triangles=tri,chart_ids=chart);receipt.update(arrays_path=str(path),arrays_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),arrays_bytes=path.stat().st_size)
    (ROOT/'receipts/native_coverage.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
