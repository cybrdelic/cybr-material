"""Encode the unchanged analytic finish-normal field natively as RGB16 PNG."""
import sys,hashlib,json,struct,zlib,resource,time,gc
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
SOURCE=Path("/workspace/scratch/b4387906eb93/technical-recovery-20261006/restored/CYBR_technical_cloud_review/source")
sys.path.insert(0,str(SOURCE))
from surface_math import field
from expansion.catalog import get
from expansion.maps import normals

def chunk(f,name,data):
    f.write(struct.pack(">I",len(data)));f.write(name);f.write(data)
    f.write(struct.pack(">I",zlib.crc32(name+data)&0xffffffff))
def encode(p,a):
    n,m,c=a.shape;assert c==3
    with p.open("wb") as f:
        f.write(b"\x89PNG\r\n\x1a\n");chunk(f,b"IHDR",struct.pack(">IIBBBBB",m,n,16,2,0,0,0));enc=zlib.compressobj(6)
        for row in a:
            data=enc.compress(b"\0"+np.rint(row*65535).astype(">u2").tobytes())
            if data:chunk(f,b"IDAT",data)
        chunk(f,b"IDAT",enc.flush());chunk(f,b"IEND",b"")
def verify(p,a):
    data=p.read_bytes();pos=8;stream=zlib.decompressobj();raw=bytearray()
    while pos<len(data):
        length=struct.unpack(">I",data[pos:pos+4])[0];name=data[pos+4:pos+8];body=data[pos+8:pos+8+length]
        assert zlib.crc32(name+body)&0xffffffff==struct.unpack(">I",data[pos+8+length:pos+12+length])[0]
        if name==b"IHDR":assert struct.unpack(">IIBBBBB",body)==(a.shape[1],a.shape[0],16,2,0,0,0)
        if name==b"IDAT":raw.extend(stream.decompress(body))
        pos+=length+12
    raw.extend(stream.flush());stride=a.shape[1]*6+1
    assert len(raw)==a.shape[0]*stride
    for i,row in enumerate(a):
        assert raw[i*stride]==0
        q=np.frombuffer(raw,dtype=">u2",count=a.shape[1]*3,offset=i*stride+1).reshape(-1,3)
        assert np.array_equal(q,np.rint(row*65535).astype(np.uint16))
    return len(data)
def angular_stats(a,bits):
    total=0.;maximum=0.;count=0
    for row0 in range(0,len(a),64):
        f=a[row0:row0+64].astype(float);u=f*2-1;u/=np.linalg.norm(u,axis=-1,keepdims=True)
        q=np.rint(f*(2**bits-1))/(2**bits-1)*2-1;q/=np.linalg.norm(q,axis=-1,keepdims=True)
        ang=np.degrees(np.arccos(np.clip(np.sum(u*q,axis=-1),-1,1)))
        total+=float(np.sum(ang**2));count+=ang.size;maximum=max(maximum,float(ang.max()))
    return dict(rms_degrees=(total/count)**.5,max_degrees=maximum)
rows=[]
for cid in ["14_chainmail","12_lattice_wire"]:
    t=time.monotonic();s=get(cid)
    h=(.5+.01*field(4096,64,s["seed"]+1)).astype("f4");a=normals(h,s);del h;gc.collect()
    path=ROOT/"maps"/(cid+"_Normal_OpenGL16.png");encode(path,a);size=verify(path,a)
    stats={str(b):angular_stats(a,b) for b in [8,16]}
    row=dict(material_id=cid,path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),resolution=[4096,4096],bits=16,color_space="Non-Color",normal_convention="OpenGL tangent RGB, original signed derivatives retained",same_analytic_normal_field=True,PNG_integer_roundtrip_exact=True,angular_quantization_error=stats,bytes=size,elapsed_s=time.monotonic()-t,peak_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
    rows.append(row);print(json.dumps(row),flush=True)
    (ROOT/"receipts/native16_maps.json").write_text(json.dumps(rows,indent=2))
    del a;gc.collect()

