"""Emit an additive plan fragment matching CATALOG_CAPTURE_CONTRACT.md.

It does not touch the active room plan or renderer. Review before integration.
"""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
targets=[];receipts={}
for candidate in sorted((ROOT/'expansion/receipts').glob('*_hero_soft*.json')):
    data=json.loads(candidate.read_text());id=data['recipe']
    if id not in receipts or candidate.name>receipts[id].name:receipts[id]=candidate
for receipt in receipts.values():
    data=json.loads(receipt.read_text());s=data['settings'];w=data['camera']['ortho_scale']/1.42
    target=[0,0,.012 if s['family'] in ('PETG','Grass') else .005]
    views=[dict(id=id,label=label,position_m=[w*x,w*y,w*z],target_m=target,lens_mm=55)
           for id,label,x,y,z in [('front','Front thickness',0,-1.8,.55),('three_quarter','Three-quarter structure',1.1,-1.5,1.25),('macro','Close structure',.35,-.55,.45)]]
    lights=[]
    for id,label,key_height,key_size in [('soft','Soft studio',1.4,1.05),('grazing','Low grazing',.20,.65),('hard','Hard structure',1.4,.09)]:
        lamps=[]
        for name,pos,power,size in [('key',[-w*.7,-w*.65,w*key_height],36*w*w,w*key_size),('fill',[w*.8,-w*.55,w*.55],5*w*w,w*1.2),('rim',[w*.55,w*.25,w*.6],22*w*w,w*.75)]:
            lamps.append(dict(name=f'Catalog / {name}',position_m=pos,target_m=target,energy=power,energy_unit='W',type='AREA',size_m=size,size_y_m=size,color=[1,1,1]))
        lights.append(dict(id=id,label=label,world_strength=.22,world_color=[.2,.2,.2],lights=lamps))
    targets.append(dict(id=s['id'],name=s['name'],category='material',family=s['family'],
                       source_scene_file=f"expansion/scenes/{receipt.stem}.blend",source_material_id=s['id'],
                       inspection_width_m=w,aspect=[1,1],exposure_ev=0,views=views,lighting=lights))
result=dict(contract='CATALOG_CAPTURE_CONTRACT.md, inspected read-only in active checkout',
            baseline='de4146ac759164e29323b152d37f64417c00c3fb',status='fragment only, not installed',
            map_migration='After review, copy expansion/maps/<id>/ to materials/<id>/ and rewrite FILE paths without changing original IDs',targets=targets)
path=ROOT/'expansion/capture_plan_fragment.json';path.write_text(json.dumps(result,indent=2)+'\n')
print('CAPTURE_FRAGMENT',len(targets),path)
