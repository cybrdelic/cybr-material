"""Audit the shared analytic-expression graph before Blender source creation."""
from surface_contract import *
from analytic_shader import NumericOps,field_expressions
import json,resource
rng=np.random.default_rng(5179);allrows=[]
for chart in range(6):
    q=rng.uniform([-.04,-.04,-.004],[.04,.04,.004],(4000,3));axis,value={0:(2,.004),1:(2,-.004),2:(0,.04),3:(0,-.04),4:(1,-.04),5:(1,.04)}[chart];q[:,axis]=value
    expected=surface(q);en=normal(q,chart);rows=[]
    for dtype in ['f8','f4']:
        r=field_expressions(NumericOps(HEIGHT,dtype),list(q.T),TANGENTS[chart,0],TANGENTS[chart,1]);p=np.stack(r['position'],-1);n=np.stack(r['normal'],-1);n=n/np.linalg.norm(n,axis=-1,keepdims=True);angle=np.degrees(np.arctan2(np.linalg.norm(np.cross(n,en),axis=-1),np.sum(n*en,axis=-1)));b=base_and_direction(q)[0];cells=np.floor(np.clip((b[:,:2]/.08+.5)*4096-.5,0,4095));cg=np.stack(r['cell'],-1);mismatch=np.any(cg!=cells,axis=-1)
        rows.append({'dtype':dtype,'finite':bool(np.isfinite(p).all() and np.isfinite(n).all()),'position_max_m':float(np.linalg.norm(p-expected,axis=-1).max()),'normal_max_deg':float(angle.max()),'normal_rms_deg':float(np.sqrt(np.mean(angle**2))),'normal_max_same_native_cell_deg':float(angle[~mismatch].max()),'native_cell_mismatch':int(mismatch.sum())})
    allrows.append({'chart':chart,'samples':len(q),'rows':rows})
r={'graph':'Same retained symbolic expression builder for CPU and Blender nodes','rows':allrows,'source_sha256':hashlib.sha256((ROOT/'source/analytic_shader.py').read_bytes()).hexdigest(),'peak_RSS_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'not_actual_renderer_verification':True};(ROOT/'receipts/shader_math.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
