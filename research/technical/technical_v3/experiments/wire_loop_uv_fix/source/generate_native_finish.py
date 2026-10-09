"""Same archived finish recipe, sequential channel evaluation, native4096."""
import sys,json,hashlib,gc,time,resource
from pathlib import Path
import numpy as np
from PIL import Image
ROOT=Path(__file__).resolve().parents[1]
SOURCE=Path("/workspace/scratch/b4387906eb93/technical-recovery-20261006/restored/CYBR_technical_cloud_review/source")
sys.path.insert(0,str(SOURCE))
from surface_math import field,clamp,color
from expansion.catalog import get
from expansion.maps import generate,normals
S=get("12_lattice_wire")
def channels(n):
    fine=field(n,64,S["seed"]+1)
    yield "Roughness",clamp(S["roughness"]+.008*fine)
    h=(.5+.01*fine).astype("f4")
    del fine
    yield "Normal_OpenGL",normals(h,S)
    del h
    broad=field(n,5,S["seed"]+2,2)
    yield "BaseColor",clamp(color(S["color"],.005*broad)).astype("f4")
    del broad
    yield "Metallic",np.full((n,n),S["metallic"],np.float32)
def main():
    t=time.monotonic()
    reference=generate(S,512);reference["Normal_OpenGL"]=normals(reference["Height"],S)
    equivalence={}
    for name,a in channels(512):
        actual=np.array(Image.open(SOURCE.parents[0]/"expansion/maps/12_lattice_wire"/(name+".png")))
        encoded=np.rint(a*255).astype("u1")
        assert np.array_equal(encoded,actual),(name,"archive mismatch")
        assert np.array_equal(encoded,np.rint(reference[name]*255).astype("u1"))
        equivalence[name]=True
    del reference;gc.collect()
    out=ROOT/"maps";out.mkdir(exist_ok=True)
    saved=[]
    for name,a in channels(4096):
        p=out/(name+".png");Image.fromarray(np.rint(a*255).astype("u1")).save(p,compress_level=4)
        saved.append(dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),resolution=[4096,4096],colorspace="sRGB" if name=="BaseColor" else "Non-Color",bits=8))
        del a;gc.collect()
    report=dict(recipe="Exact archived authored finish. No measured formation simulation.",same_source512_quantized_pixels=equivalence,maps=saved,elapsed_s=time.monotonic()-t,peak_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,physical_mapping_limit="Normalized closed chart. Source60mm finish period is not reinterpreted as a calibrated wire-scale finish.")
    (ROOT/"receipts/native_maps.json").write_text(json.dumps(report,indent=2))
    print(json.dumps(report),flush=True)
if __name__=="__main__": main()

