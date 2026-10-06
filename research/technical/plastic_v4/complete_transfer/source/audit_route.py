"""Bounded planning audit only: no mesh artifact, Blender build, or renderer."""
import json, hashlib, resource, time
import numpy as np
from surface_contract import *

def check_derivatives():
    rng=np.random.default_rng(20261006);q=rng.uniform([-.04,-.04,-.004],[.04,.04,.004],(6000,3));chart=np.arange(len(q))%6
    for k in range(6):
        ids=chart==k;fixed={0:(2,.004),1:(2,-.004),2:(0,.04),3:(0,-.04),4:(1,-.04),5:(1,.04)}[k];q[ids,fixed[0]]=fixed[1]
    # Include the retained parting witness and ejector annulus deliberately.
    q=np.r_[q,np.column_stack([np.full(1001,.04),np.zeros(1001),np.linspace(-.00016,.00016,1001)]),np.column_stack([.023+.0022+np.linspace(-.0003,.0003,1001),np.full(1001,.023),np.full(1001,-.004)])]
    chart=np.r_[chart,np.full(1001,2),np.full(1001,1)];rows=[]
    for step in [1e-8,2e-9]:
        errors=[];derr=[]
        for axis in [0,1]:
            t=TANGENTS[chart,axis];db,dd=macro_differential(q,t);bp,dp=base_and_direction(q+step*t);bm,dm=base_and_direction(q-step*t)
            errors.append(np.linalg.norm((bp-bm)/(2*step)-db,axis=-1));derr.append(np.linalg.norm((dp-dm)/(2*step)-dd,axis=-1))
        rows.append({'step_m':step,'base_derivative_max':float(np.max(errors)),'direction_derivative_max_per_m':float(np.max(derr))})
    import sys
    sys.path.insert(0,str(EXPERIMENTS/'plastic_explicit_curve/height_sampling_audit/seam_revision/source'))
    from cpu_surface_contract import expected_world_normal
    qx=rng.uniform([.0398,-.0381,.0038],[.03999,-.0379,.00399],(2000,3));qx[:1000,2]=.004;qx[1000:,0]=.04
    ids=np.r_[np.zeros(1000,dtype='i1'),np.ones(1000,dtype='i1')];old=expected_world_normal(qx,ids);new=normal(qx,np.where(ids==0,0,2));angle=np.degrees(np.arctan2(np.linalg.norm(np.cross(old,new),axis=-1),np.sum(old*new,axis=-1)))
    return {'finite_difference_checks':rows,'frozen_outer_corner_normal_crosscheck_max_degrees':float(angle.max()),'frozen_outer_corner_sample_count':len(qx)}

def scan_lod(n,subdivision=8):
    # Every top facet receives an independent interior point. Densely probe a
    # deterministic stratified subset, plus every sampled native field extremum.
    ax=np.unique(np.r_[np.linspace(-.04,.04,n+1),-.037,.037]);ny=len(ax)-1;rng=np.random.default_rng(8675309)
    maxerr=0.;sumsq=0.;count=0;maxmacro=0.;quant=[]
    for first in range(0,ny,8):
        X,Y=np.meshgrid(ax[:-1],ax[first:min(first+8,ny)]);DX,DY=np.meshgrid(np.diff(ax),np.diff(ax)[first:first+8]);uv=np.column_stack([X.ravel(),Y.ravel()]);du=np.column_stack([DX.ravel(),np.zeros(DX.size)]);dv=np.column_stack([np.zeros(DY.size),DY.ravel()])
        for off in [(np.zeros(2),du,du+dv),(np.zeros(2),du+dv,dv)]:
            q=np.empty((len(uv),3,3));q[:,:,:2]=np.stack([uv+o for o in off],1);q[:,:,2]=.004;v=surface(q).astype('f4').astype('f8');b=base_and_direction(q)[0]
            # Random barycentric sample plus centroid; no claimed supremum.
            for mode in [0,1]:
                w=rng.exponential(size=(len(q),3)) if mode else np.ones((len(q),3));w/=w.sum(1,keepdims=True);query=np.einsum('tij,ti->tj',q,w);expected=surface(query);actual=np.einsum('tij,ti->tj',v,w);e=np.linalg.norm(actual-expected,axis=-1);m=np.linalg.norm(np.einsum('tij,ti->tj',b,w)-base_and_direction(query)[0],axis=-1)
                maxerr=max(maxerr,float(e.max()));maxmacro=max(maxmacro,float(m.max()));sumsq+=float(e@e);count+=len(e);quant.append(e[::16])
    z=np.unique(np.r_[np.linspace(-.004,-.001,31),np.linspace(-.001,.001,21),np.linspace(.001,.004,31),np.linspace(-.00016,.00016,65)])
    N=len(ax)-1;K=len(z)-1;verts=2*(N+1)**2+4*N*(K-1);tri=4*N*N+8*N*K
    return {'regular_axis_cells':n,'actual_axis_cells':N,'side_z_cells':K,'closed_mesh_vertices':int(verts),'closed_mesh_triangles':int(tri),'top_samples':count,'sampled_top_position_max_m':maxerr,'sampled_top_position_rms_m':float(np.sqrt(sumsq/count)),'sampled_top_position_p99_m':float(np.quantile(np.concatenate(quant),.99)),'sampled_top_macro_max_m':maxmacro,'sampled_projection_upper_bound_pixels':{'full_960':maxerr/.12960000336170197*960,'macro_960':maxerr/.024000000208616257*960},'packed_mesh_and_original_q_MiB':float((verts*24+tri*12)/2**20),'estimated_blender_build_RSS_MiB':[int(300+(verts*160+tri*120)/2**20),int(500+(verts*240+tri*200)/2**20)],'limitation':'Two barycentric probes per top triangle; sides/underside and worst native boundary positions still require separate checks. Not a bound or20nm qualification.'}

if __name__=='__main__':
    start=time.monotonic();res={'reconstruction':True,'not_historical_binary_recovery':True,'height_sha256':MAP_SHA,'macro_sha256':MACRO_SHA,'derivatives':check_derivatives(),'lods':[]}
    for n in [256,512,1024]:
        row=scan_lod(n);res['lods'].append(row);print(json.dumps(row),flush=True)
    res.update(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),surface_contract_sha256=hashlib.sha256((ROOT/'source/surface_contract.py').read_bytes()).hexdigest(),peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,seconds=time.monotonic()-start)
    (ROOT/'receipts/route_audit.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res),flush=True)
