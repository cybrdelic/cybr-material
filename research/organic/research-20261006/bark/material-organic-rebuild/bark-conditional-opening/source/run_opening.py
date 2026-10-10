"""Checkpointed conditional quasistatic evolution, rejecting unresolved work."""
import argparse,hashlib,json,os,resource,time
from pathlib import Path
import numpy as np
from opening_model import OpeningModel,History,ROOT
from checkpoint_io import write_checkpoint,read_checkpoint


def run(label,end_load=1.,increment=.05,radial_pitch=.001,quadrature=4,penalty=None,Gc=None):
    start=time.monotonic();dest=ROOT/'state'/label;dest.mkdir(parents=True,exist_ok=True)
    dependencies=[ROOT/'source'/name for name in ('opening_model.py','volume_coupon.py','periodic_sector.py','single_seam.py','checkpoint_io.py','run_opening.py')]+list((ROOT/'config').glob('*.json'))
    params=dict(label=label,end_load=end_load,initial_increment=increment,radial_pitch=radial_pitch,quadrature=quadrature,penalty=penalty,Gc=Gc,
                source_hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in dependencies})
    m=OpeningModel(radial_pitch,quadrature,penalty,Gc);u=np.zeros_like(m.X);accepted=m.initial_history;_,_,_,initial,_=m.evaluate(u,accepted,work=True)
    trace=[dict(station=0,load=0.,work_J=0.,energy_work_error_J=0.,**{k:v for k,v in initial.items() if k!='load'})];next_increment=increment;events=[]
    file=dest/'restart.npz'
    if file.exists():
        a,old=read_checkpoint(file);assert old['parameters']==params,'Frozen source/parameter mismatch'
        u=a['displacement'];accepted=History.make(a['maximum'],a['damage']);trace=old['trace'];events=old['events'];next_increment=old['next_increment'];m.set_load(trace[-1]['load'])
        assert abs(m.evaluate(u,accepted)[3]['stored_energy_J']-trace[-1]['stored_energy_J'])<1e-12
    def save(status,error=None):
        meta=dict(status=status,error=error,parameters=params,trace=trace,events=events,next_increment=next_increment,
                  reference_mass_kg=m.seam.mass,prepared_slit_phi_J=m.initial_phi,
                  wall_s_this_process=time.monotonic()-start,peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
                  cells=len(m.base.cells),full_nodes=len(m.base.X),periodic_nodes=len(m.X))
        arrays=dict(displacement=u,maximum=accepted.maximum,damage=accepted.damage,full_displacement=m.seam.expand(u),
                    reference_nodes=m.base.X,node_map=m.seam.cell_map,reference_cell_nodes=m.base.Xcell,periodic_groups=m.seam.groups,rotations=m.seam.rotations,
                    interface_reference=m.interface_reference,interface_area=m.area,starter=m.starter,initial_maximum=m.initial_history.maximum,
                    tri=m.base.tri,points=m.base.points)
        meta['state_sha256']=write_checkpoint(file,arrays,meta)
        temp=dest/'restart.tmp.json';temp.write_text(json.dumps(meta,indent=2));os.replace(temp,dest/'restart.json')
        return meta
    save('running')
    while trace[-1]['load']<end_load-1e-12:
        old=trace[-1];h=min(next_increment,end_load-old['load']);load=old['load']+h;begin=time.monotonic()
        try:
            m.set_load(load);candidate,new_history,result,opening=m.equilibrate(u,accepted,deadline=start+90)
        except (ValueError,RuntimeError,TimeoutError) as exc:
            if m.last_trial is not None:
                t=m.last_trial;write_checkpoint(dest/'unaccepted_trial.npz',{k:v for k,v in t.items() if k!='report'},dict(report=t['report'],reason=str(exc),attempted_load=load,accepted_load=old['load'],accepted=False))
            m.set_load(old['load']);status='paused' if isinstance(exc,TimeoutError) else 'blocked'
            events.append(dict(load=load,event=status,error=str(exc)))
            meta=save(status,str(exc));(ROOT/'receipts'/f'{label}.json').write_text(json.dumps(meta,indent=2));print(status.upper(),json.dumps(events[-1]),flush=True);return
        dw=.5*(old['eigenstrain_conjugate_J']+result['eigenstrain_conjugate_J'])*h
        de=result['stored_energy_J']-old['stored_energy_J'];dd=result['new_dissipation_from_initial_J']-old['new_dissipation_from_initial_J']
        error=de+dd-dw;relative=abs(error)/max(abs(dw),1e-12)
        if relative>m.parameters['per_step_ledger_relative_tolerance']:
            events.append(dict(event='work_rejected',attempted_load=load,increment=h,error_J=error,error_relative=relative,
                               stored_energy_change_J=de,new_dissipation_J=dd,work_J=dw))
            write_checkpoint(dest/'unaccepted_trial.npz',dict(displacement=candidate,maximum=new_history.maximum,damage=new_history.damage,opening=opening),dict(report=result,accepted=False,event=events[-1]))
            next_increment=h/2;m.set_load(old['load']);save('running')
            print('WORK_REJECTED',json.dumps(events[-1]),flush=True)
            if next_increment<m.parameters['minimum_load_increment']:
                meta=save('blocked','Unresolved work balance at minimum load increment');(ROOT/'receipts'/f'{label}.json').write_text(json.dumps(meta,indent=2));return
        else:
            write_checkpoint(dest/'previous_accepted.npz',dict(displacement=u,maximum=accepted.maximum,damage=accepted.damage),dict(load=old['load'],work_J=old['work_J'],purpose='Prior immutable history for the final station incremental-potential/stability audit'))
            u=candidate;accepted=new_history;work=old['work_J']+dw
            balance=result['stored_energy_J']+result['new_dissipation_from_initial_J']-work
            trace.append(dict(station=old['station']+1,work_J=work,energy_work_error_J=balance,step_work_error_J=error,step_work_relative_error=relative,
                              increment=h,station_wall_s=time.monotonic()-begin,**result))
            save('running');print('ACCEPTED',json.dumps({k:trace[-1][k] for k in ('station','load','iterations','maximum_opening_m','maximum_new_damage','stored_energy_J','new_dissipation_from_initial_J','energy_work_error_J','station_wall_s')}),flush=True)
        if time.monotonic()-start>88:save('paused');print('PAUSED',flush=True);return
    meta=save('completed');(ROOT/'receipts'/f'{label}.json').write_text(json.dumps(meta,indent=2));print('COMPLETED',json.dumps({k:meta[k] for k in ('wall_s_this_process','peak_rss_MiB','cells','reference_mass_kg')}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('label');p.add_argument('--end-load',type=float,default=1.);p.add_argument('--increment',type=float,default=.05)
    p.add_argument('--radial-pitch',type=float,default=.001);p.add_argument('--quadrature',type=int,default=4);p.add_argument('--penalty',type=float);p.add_argument('--Gc',type=float)
    a=p.parse_args();run(a.label,a.end_load,a.increment,a.radial_pitch,a.quadrature,a.penalty,a.Gc)
