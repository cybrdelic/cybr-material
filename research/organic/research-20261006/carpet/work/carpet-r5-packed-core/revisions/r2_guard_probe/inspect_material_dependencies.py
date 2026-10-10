# Explicit external inputs replace machine-specific historical file locations.
from pathlib import Path as _ResearchInputPath
import sys as _research_input_sys
for _research_parent in _ResearchInputPath(__file__).resolve().parents:
    if (_research_parent / "research_inputs.py").is_file():
        _research_input_sys.path.insert(0, str(_research_parent))
        break
from research_inputs import required_input as _required_research_input
import bpy,json,sys,resource
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'revisions/r1_frozen/source'))
from inspect_baseline import material_sha
names=[f'Carpet / fine wool yarn {i}' for i in range(4)]
with bpy.data.libraries.load(str(_required_research_input("baseline-scene")),link=False) as (src,dst):dst.materials=list(names)
h=json.loads((ROOT/'revisions/r1_frozen/receipts/packed_core_handoff.json').read_text());materials=[]
for name in names:
    m=bpy.data.materials[name];digest=material_sha(m);assert digest==h['retained_materials_sha256'][name]
    nodes=[{'name':n.name,'type':n.bl_idname} for n in m.node_tree.nodes]
    assert not any('Hair' in n['type'] or n['type']=='ShaderNodeAttribute' for n in nodes)
    links=[[l.from_node.bl_idname,l.from_socket.name,l.to_node.bl_idname,l.to_socket.name] for l in m.node_tree.links]
    materials.append({'name':name,'node_inventory':nodes,'links':links,'sha256_matches_r1':True,'sha256':digest})
report={'source_materials':'exact selected baseline material data; digests match frozen r1 candidate','materials':materials,'peak_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'no_scene_saved':True,'no_render':True}
(HERE/'material_dependency_inventory.json').write_text(json.dumps(report,indent=2));print(json.dumps({'materials':len(materials),'all_digests_match':True,'node_types':sorted(set(n['type'] for m in materials for n in m['node_inventory'])),'peak_rss_mib':report['peak_rss_mib']}),flush=True)
