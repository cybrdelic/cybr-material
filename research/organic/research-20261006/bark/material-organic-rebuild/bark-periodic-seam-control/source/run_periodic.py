import argparse,hashlib,json,os,resource,time
from pathlib import Path
import numpy as np
from periodic_sector import PeriodicSector,ROOT
from checkpoint_io import write_checkpoint,read_checkpoint

def run(label,multiplier=1,steps=16,end_load=1.):
    start=time.monotonic();dest=ROOT/'state'/label;dest.mkdir(parents=True,exist_ok=True)
    files=list((ROOT/'source').glob('*.py'))+[ROOT/'config/experiment.json']
    params=dict(label=label,multiplier=multiplier,pitch_m=.006,steps=steps,end_load=end_load,source_hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    m=PeriodicSector(sector_multiplier=multiplier);u=np.zeros_like(m.X);_,_,initial=m.evaluate(u,work=True)
    trace=[dict(station=0,iterations=0,relative_energy_residual=0.,eigenstrain_work_J=0.,energy_work_error_J=0.,**initial)]
    file=dest/'restart.npz'
    if file.exists():
        a,metadata=read_checkpoint(file);assert metadata['parameters']==params;u=a['displacement'];trace=metadata['trace'];m.set_load(trace[-1]['load'])
        assert abs(m.evaluate(u)[0]-trace[-1]['stored_energy_J'])<1e-12
    def save(status,error=None):
        meta=dict(status=status,error=error,parameters=params,trace=trace,reference_mass_kg=m.mass,reference_volume_m3=m.volume,
                  sector_angle_rad=m.alpha,full_nodes=len(m.base.X),periodic_nodes=len(m.X),cells=len(m.base.cells),
                  wall_s=time.monotonic()-start,peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
        meta['state_sha256']=write_checkpoint(file,dict(displacement=u,full_displacement=m.expand(u),reference_nodes=m.base.X,node_map=m.base.map,
            reference_cell_nodes=m.base.Xcell,points=m.base.points,tri=m.base.tri,periodic_groups=m.groups,rotations=m.rotations),meta)
        temp=dest/'restart.tmp.json';temp.write_text(json.dumps(meta,indent=2));os.replace(temp,dest/'restart.json')
        return meta
    save('running')
    for station in range(trace[-1]['station']+1,steps+1):
        load=end_load*station/steps
        try:m.set_load(load);candidate,result=m.equilibrate(u)
        except (ValueError,RuntimeError) as exc:
            m.set_load(trace[-1]['load']);meta=save('blocked',str(exc));print('BLOCKED',str(exc),flush=True);return
        old=trace[-1];work=old['eigenstrain_work_J']+.5*(result['eigenstrain_conjugate_J']+old['eigenstrain_conjugate_J'])*(load-old['load'])
        u=candidate;trace.append(dict(station=station,eigenstrain_work_J=work,energy_work_error_J=result['stored_energy_J']-work,**result));save('running')
        print('STATION',station,load,result['stored_energy_J'],result['iterations'],flush=True)
        if time.monotonic()-start>90:print('PAUSED',flush=True);return
    meta=save('completed');(ROOT/'receipts'/f'{label}.json').write_text(json.dumps(meta,indent=2));print('COMPLETED',json.dumps({k:meta[k] for k in ('wall_s','peak_rss_MiB','reference_mass_kg','cells')}),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('label');p.add_argument('--multiplier',type=int,default=1);p.add_argument('--steps',type=int,default=16);p.add_argument('--end-load',type=float,default=1.);a=p.parse_args();run(a.label,a.multiplier,a.steps,a.end_load)
