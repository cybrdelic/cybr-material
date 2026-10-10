# Explicit external inputs replace machine-specific historical file locations.
from pathlib import Path as _ResearchInputPath
import sys as _research_input_sys
for _research_parent in _ResearchInputPath(__file__).resolve().parents:
    if (_research_parent / "research_inputs.py").is_file():
        _research_input_sys.path.insert(0, str(_research_parent))
        break
from research_inputs import required_input as _required_research_input
import bpy,numpy as np,resource,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
rss=lambda:resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
current=lambda:int(next(x.split()[1] for x in Path('/proc/self/status').read_text().splitlines() if x.startswith('VmRSS:')))/1024
bpy.ops.wm.open_mainfile(filepath=str(_required_research_input("baseline-scene")))
before=rss();current_before=current();p=np.load(ROOT/'receipts/probe_positions.npy');curves=p.shape[0]*p.shape[1];points=p.shape[2]
data=bpy.data.hair_curves.new('Memory probe packed yarn');data.add_curves([points]*curves);data.position_data.foreach_set('vector',p.ravel())
rr=data.attributes.new('radius','FLOAT','POINT');rr.data.foreach_set('value',np.full(curves*points,14e-6,dtype='f4'))
typ=data.attributes.new('curve_type','INT8','CURVE');typ.data.foreach_set('value',np.zeros(curves,dtype='i4'))
res=data.attributes.new('resolution','INT','CURVE');res.data.foreach_set('value',np.full(curves,3,dtype='i4'))
ob=bpy.data.objects.new('Memory probe packed yarn',data);bpy.context.scene.collection.objects.link(ob)
data.materials.append(bpy.data.materials['Carpet / fine wool yarn 0'])
bpy.context.view_layer.update();evaluated=ob.evaluated_get(bpy.context.evaluated_depsgraph_get());n=len(evaluated.data.points)
report={'baseline_load_peak_mib':before,'probe_peak_mib':rss(),'native_curves':curves,'native_points':n,'point_payload_bytes':curves*points*16,'point_budget_full':1764*187*192,'full_payload_mib':1764*187*192*16/2**20,'evaluated_representation':'native Catmull-Rom Curves, no conversion to mesh','cycles_bvh_memory':'not measured by builder; rendering belongs to central renderer'}
report['rss_before_mib']=current_before;report['rss_after_mib']=current();report['rss_added_mib']=current()-current_before
(ROOT/'receipts/native_curve_memory_probe.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
