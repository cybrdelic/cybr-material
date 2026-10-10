import resource,os
resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json,numpy as np
from newton_force import MatrixFreeBundle
P=Path(__file__).resolve().parents[1];c=MatrixFreeBundle(19,128);c.controls(.00238,0);c.w[:]=np.load(P/'data/yarn_19_direct_r2.npz')['w'];c.checkpoint_path=P/'data/yarn_19_force_accepted_checkpoint.npz';r=c.solve(.00238,0,wall_s=60,maxiter=140);r['peak_rss_MiB']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024;r['scope']='Same 19-fibre stock, intrinsic crimp and root controls; optimizer resumed with exact sparse physical Newton. No finite forming friction history.';np.savez_compressed(P/'data/yarn_19_force.npz',w=c.w,points=c.points(),material_arc_m=c.rod.s,endpoint_fixed=c.endpoint_fixed,intrinsic_curvature_per_m=np.array([r.k0 for r in c.rods]),roots=np.stack([c.origins,c.targets],axis=1));(P/'receipts/yarn_19_force.json').write_text(json.dumps(r,indent=2));print('RESULT',json.dumps({k:v for k,v in r.items() if k not in ['solver_iteration_trace','root_positions_m','root_forces_world_N','clamp_torques_world_Nm']}),flush=True)
