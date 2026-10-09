"""Closed six-chart camera-dependent approximation of the retained surface.
Explicitly not a 20nm full-part representation. Native state is never filtered.
"""
from surface_contract import *
import json, time, resource
POSITION_PIXELS=.25
REFINE_PIXELS=.14
RESOLUTION=960
MAX_TRIANGLES=2000000
PROVENANCE=json.loads((EXPERIMENTS/'plastic_native_relief_transfer/inputs/production_scene.json').read_text())
CAMERAS={k:v for k,v in PROVENANCE['studio_signature']['objects'].items() if v['type']=='CAMERA'}


def seed_mesh(n=128):
    axis=np.unique(np.r_[np.linspace(-.04,.04,n+1),-.037,.037])
    zs=np.unique(np.r_[np.linspace(-.004,-.001,25),np.linspace(-.001,.001,9),np.linspace(.001,.004,25),np.linspace(-.00016,.00016,65)])
    nx=len(axis)-1;nz=len(zs)-1;keys={};q=[];tri=[];charts=[]
    def vert(i,j,k):
        key=(i,j,k)
        if key not in keys:keys[key]=len(q);q.append((axis[i],axis[j],zs[k]))
        return keys[key]
    def quad(v,c):
        tri.extend([(v[0],v[1],v[2]),(v[0],v[2],v[3])]);charts.extend([c,c])
    for j in range(nx):
        for i in range(nx):
            quad([vert(i,j,nz),vert(i+1,j,nz),vert(i+1,j+1,nz),vert(i,j+1,nz)],0)
            quad([vert(i,j,0),vert(i,j+1,0),vert(i+1,j+1,0),vert(i+1,j,0)],1)
    for k in range(nz):
        for i in range(nx):
            quad([vert(nx,i,k),vert(nx,i+1,k),vert(nx,i+1,k+1),vert(nx,i,k+1)],2)
            quad([vert(0,i,k),vert(0,i,k+1),vert(0,i+1,k+1),vert(0,i+1,k)],3)
            quad([vert(i,0,k),vert(i+1,0,k),vert(i+1,0,k+1),vert(i,0,k+1)],4)
            quad([vert(i,nx,k),vert(i,nx,k+1),vert(i+1,nx,k+1),vert(i+1,nx,k)],5)
    return np.array(q,'f8'),np.array(tri,'i4'),np.array(charts,'i1')


def weights(dense=False):
    if not dense:return np.array([[1,1,1],[1.5,1.5,0],[1.5,0,1.5],[0,1.5,1.5],[2, .5,.5],[.5,2,.5],[.5,.5,2]])/3
    grid=np.array([(i,j,6-i-j) for i in range(7) for j in range(7-i)],dtype='f8')/6
    rng=np.random.default_rng(1179);r=rng.exponential(size=(12,3));r/=r.sum(axis=1,keepdims=True)
    return np.r_[grid,r]


def errors(q,tri,dense=False):
    p=surface(q).astype('f4').astype('f8');maxpx=np.zeros(len(tri));report={k:{'max_absolute_m':0.,'max_projected_pixels':0.,'samples':0} for k in CAMERAS};absolute=0.
    for start in range(0,len(tri),4096):
        t=tri[start:start+4096];tq=q[t];tp=p[t]
        for w in weights(dense):
            x=np.einsum('tvk,v->tk',tq,w);exact=surface(x);actual=np.einsum('tvk,v->tk',tp,w);delta=actual-exact;absolute=max(absolute,float(np.linalg.norm(delta,axis=-1).max()))
            for name,c in CAMERAS.items():
                m=np.array(c['matrix']).reshape(4,4);field=c['data']['ortho_scale'];rel=exact-m[:3,3];xy=rel@m[:3,:2]
                # Include 0.75mm guard beyond view so silhouettes remain sampled.
                visible=np.max(np.abs(xy),axis=1)<field/2+.00075
                projected=np.linalg.norm(delta@m[:3,:2],axis=-1)/field*RESOLUTION
                # Refinement uses conservative 3D distance; receipt reports both.
                bound=np.linalg.norm(delta,axis=-1)/field*RESOLUTION
                maxpx[start:start+len(t)]=np.maximum(maxpx[start:start+len(t)],np.where(visible,bound,0))
                if visible.any():
                    row=report[name];row['max_absolute_m']=max(row['max_absolute_m'],float(np.linalg.norm(delta[visible],axis=-1).max()));row['max_projected_pixels']=max(row['max_projected_pixels'],float(projected[visible].max()));row['samples']+=int(visible.sum())
    return maxpx,{'cameras':report,'all_chart_sampled_max_absolute_m':absolute,'max_camera_conservative_pixels':float(maxpx.max()),'samples_per_triangle':len(weights(dense)),'passes_sampled_camera_budget':bool(maxpx.max()<=POSITION_PIXELS)}


def refine(q,tri,chart,mark):
    # Unique undirected edges support shared mids across all six chart boundaries.
    e=np.sort(np.stack([tri[:,[0,1]],tri[:,[1,2]],tri[:,[2,0]]],axis=1),axis=-1)
    codes=e[...,0].astype('i8')*len(q)+e[...,1];unique,inv=np.unique(codes,return_inverse=True);inv=inv.reshape(-1,3)
    split=np.zeros(len(unique),bool);split[inv[mark].ravel()]=True
    mids=np.full(len(unique),-1,'i4');mids[split]=np.arange(len(q),len(q)+split.sum());selected=unique[split];pair=np.column_stack([selected//len(q),selected%len(q)])
    q=np.r_[q,q[pair].mean(axis=1)];mi=mids[inv];bits=(mi>=0)@np.array([1,2,4]);new=[];cs=[]
    for pat in range(8):
        ids=np.flatnonzero(bits==pat);a,b,c=tri[ids].T;ab,bc,ca=mi[ids].T
        if pat==0:ts=[(a,b,c)]
        elif pat==1:ts=[(a,ab,c),(ab,b,c)]
        elif pat==2:ts=[(b,bc,a),(bc,c,a)]
        elif pat==4:ts=[(c,ca,b),(ca,a,b)]
        elif pat==3:ts=[(b,bc,ab),(a,ab,c),(ab,bc,c)]
        elif pat==6:ts=[(c,ca,bc),(b,bc,a),(bc,ca,a)]
        elif pat==5:ts=[(a,ab,ca),(c,ca,b),(ca,ab,b)]
        else:ts=[(a,ab,ca),(ab,b,bc),(ca,bc,c),(ab,bc,ca)]
        for t in ts:new.append(np.column_stack(t));cs.append(chart[ids])
    return q,np.concatenate(new),np.concatenate(cs)


def topology(q,tri):
    e=np.concatenate([tri[:,[0,1]],tri[:,[1,2]],tri[:,[2,0]]]);codes=np.minimum(e[:,0],e[:,1]).astype('i8')*len(q)+np.maximum(e[:,0],e[:,1]);order=np.argsort(codes);ordered=codes[order]
    closed=bool(np.all(ordered[::2]==ordered[1::2]) and np.all(ordered[1:-1:2]!=ordered[2::2]));oriented=bool(np.all(e[order[::2],0]==e[order[1::2],1]))
    p=surface(q);cross=np.cross(p[tri[:,1]]-p[tri[:,0]],p[tri[:,2]]-p[tri[:,0]]);signed_volume=float(np.einsum('ij,ij->',p[tri[:,0]],cross)/6)
    return {'closed_two_incident_edges':closed,'opposite_oriented_incidence':oriented,'euler_characteristic':int(len(q)-len(e)//2+len(tri)),'minimum_double_triangle_area_m2':float(np.linalg.norm(cross,axis=-1).min()),'positive_volume_m3':signed_volume,'bounds_m':[p.min(0).tolist(),p.max(0).tolist()],'floor_z_m':-1e-6,'minimum_floor_clearance_m':float(p[:,2].min()+1e-6)}

if __name__=='__main__':
    start=time.monotonic();q,tri,chart=seed_mesh();history=[]
    for iteration in range(16):
        e,r=errors(q,tri);row={'iteration':iteration,'vertices':len(q),'triangles':len(tri),**r};history.append(row);print(json.dumps(row),flush=True)
        mark=e>REFINE_PIXELS
        if not mark.any():break
        # Conservative upper estimate before allocating the next refinement.
        if len(tri)+3*mark.sum()>MAX_TRIANGLES:raise RuntimeError('Refinement estimate exceeds2M triangle cap; no build admitted.')
        q,tri,chart=refine(q,tri,chart,mark)
        if len(tri)>MAX_TRIANGLES:raise RuntimeError('Conforming refinement exceeds2M triangle cap; no build admitted.')
    for dense_iteration in range(6):
        e,check=errors(q,tri,dense=True)
        print(json.dumps({'dense_iteration':dense_iteration,'triangles':len(tri),**check}),flush=True)
        if check['passes_sampled_camera_budget']:break
        q,tri,chart=refine(q,tri,chart,e>REFINE_PIXELS)
        assert len(tri)<=MAX_TRIANGLES
    topo=topology(q,tri)
    res={'reconstruction':True,'representation':'Camera-dependent explicit retained surface approximation; exact analytic native field planned for shading','not_global20nm_qualification':True,'position_budget_pixels':POSITION_PIXELS,'refine_target_pixels':REFINE_PIXELS,'resolution':RESOLUTION,'absolute_camera_budgets_m':{k:v['data']['ortho_scale']/RESOLUTION*POSITION_PIXELS for k,v in CAMERAS.items()},'history':history,'final_dense':check,'topology':topo,'vertices':len(q),'triangles':len(tri),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'surface_contract_sha256':hashlib.sha256((ROOT/'source/surface_contract.py').read_bytes()).hexdigest(),'peak_RSS_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'seconds':time.monotonic()-start,'audit_limit':'Deterministic barycentric numerical sampling, not a supremum bound. Native map extrema between probes remain an explicit risk until further texel-based audit.'}
    (ROOT/'receipts/adaptive_plan.json').write_text(json.dumps(res,indent=2)+'\n')
    if check['passes_sampled_camera_budget'] and topo['closed_two_incident_edges'] and topo['opposite_oriented_incidence']:
        path=ROOT/'receipts/camera_lod_mesh.npz';np.savez_compressed(path,q_m=q,position_m=surface(q).astype('f4'),triangles=tri,chart_ids=chart);res['arrays_path']=str(path);res['arrays_sha256']=hashlib.sha256(path.read_bytes()).hexdigest();res['arrays_bytes']=path.stat().st_size
        (ROOT/'receipts/adaptive_plan.json').write_text(json.dumps(res,indent=2)+'\n')
    print(json.dumps(res),flush=True)
