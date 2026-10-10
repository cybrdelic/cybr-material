"""Bounded checkpointed quasistatic loading with explicit eigenstrain ledger."""
import argparse,hashlib,json,os,resource,time
from pathlib import Path
import numpy as np
from model import Coupon,ROOT,ORGANIC,DomainError


def run(label,pitch,steps,end_load=1.):
    start=time.monotonic(); destination=ROOT/'state'/label; destination.mkdir(parents=True,exist_ok=True)
    paths=[Path(__file__),Path(__file__).with_name('model.py'),ROOT/'config/experiment.json',
           ORGANIC/'bark-directional-native-port/source/native_directional.py',
           ORGANIC/'bark-directional-native-port/source/libdirectional_native.so',
           ORGANIC/'bark-constitutive-upgrade/source/directional_candidate.py',
           ORGANIC/'causal-bark/volume-core/source/quadratic_prism.py',
           ORGANIC/'causal-bark/volume-core/source/birth_cell.py']
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    params=dict(label=label,pitch_m=pitch,steps=steps,end_load=end_load,source_hashes=hashes)
    m=Coupon(pitch); u=np.zeros_like(m.X); _,_,initial=m.evaluate(u,work=True)
    trace=[dict(station=0,iterations=0,relative_energy_residual=0.,eigenstrain_work_J=0.,energy_work_error_J=0.,**initial)]
    restart=destination/'restart.npz'; meta_path=destination/'restart.json'
    if restart.exists():
        old=json.loads(meta_path.read_text()); assert old['parameters']==params,'Restart parameters/source mismatch'
        a=np.load(restart); u=a['displacement'].copy(); trace=old['trace']; m.set_load(trace[-1]['load'])
        assert np.array_equal(a['node_map'],m.map)
        energy,_,_=m.evaluate(u)
        assert abs(energy-trace[-1]['stored_energy_J'])<1e-12
    def save(status,blocker=None):
        target=destination/'restart.tmp.npz'
        np.savez_compressed(target,displacement=u,reference_nodes=m.X,node_map=m.map,tri=m.tri,points=m.points,
                            reference_cell_nodes=m.Xcell,cohort_thickness_m=m.config['cohort_thickness_m'])
        os.replace(target,restart)
        record=dict(parameters=params,status=status,blocker=blocker,trace=trace,state_sha256=hashlib.sha256(restart.read_bytes()).hexdigest(),
                    wall_s_this_process=time.monotonic()-start,peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
                    build_s=m.build_s,reference_mass_kg=m.mass,reference_volume_m3=m.volume,nodes=len(m.X),cells=len(m.cells),
                    dependencies_unchanged=all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==v for p,v in hashes.items()))
        temp=destination/'restart.tmp.json'; temp.write_text(json.dumps(record,indent=2)); os.replace(temp,meta_path)
        return record
    save('running')
    for station in range(trace[-1]['station']+1,steps+1):
        load=end_load*station/steps
        try:
            m.set_load(load)
            candidate,result=m.equilibrate(u)
        except (ValueError,RuntimeError) as error:
            m.set_load(trace[-1]['load'])
            record=save('blocked',str(error))
            (ROOT/'receipts'/f'{label}.json').write_text(json.dumps(record,indent=2))
            print('BLOCKED',json.dumps(dict(label=label,station=station,load=load,error=str(error))),flush=True)
            return
        old=trace[-1]
        work=old['eigenstrain_work_J']+.5*(result['eigenstrain_conjugate_J']+old['eigenstrain_conjugate_J'])*(load-old['load'])
        u=candidate
        trace.append(dict(station=station,eigenstrain_work_J=work,energy_work_error_J=result['stored_energy_J']-work,**result))
        save('running')
        print('STATION',json.dumps({k:trace[-1][k] for k in ('station','load','iterations','stored_energy_J','energy_work_error_J','maximum_tangential_strain_by_cohort')}),flush=True)
        if time.monotonic()-start>m.config['process_wall_budget_s']:
            print('PAUSED',json.dumps(dict(label=label,station=station)),flush=True);return
    record=save('completed')
    (ROOT/'receipts'/f'{label}.json').write_text(json.dumps(record,indent=2))
    print('COMPLETED',json.dumps({k:record[k] for k in ('status','wall_s_this_process','peak_rss_MiB','reference_mass_kg','nodes','cells')}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('label');p.add_argument('--pitch',type=float,default=.006);p.add_argument('--steps',type=int,default=16);p.add_argument('--end-load',type=float,default=1.)
    a=p.parse_args();run(a.label,a.pitch,a.steps,a.end_load)
