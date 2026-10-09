from pathlib import Path
import json,hashlib
R=Path(__file__).resolve().parents[1]
e=json.loads((R/'expected_native.json').read_text())
for key,value in e['png_sha256'].items():
 name='Normal_Object_RGB16.png' if key=='Normal_Object' else key+'.png'
 assert hashlib.sha256((R/'native4096'/name).read_bytes()).hexdigest()==value,name
assert hashlib.sha256((R/'native4096/GeometryHeight.npy').read_bytes()).hexdigest()==e['GeometryHeight.npy_sha256']
print('All four native4096 PNG files and macro geometry match the reviewed source exactly.')
