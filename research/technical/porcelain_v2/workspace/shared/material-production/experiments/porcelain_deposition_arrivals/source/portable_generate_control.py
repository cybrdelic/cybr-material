"""Native4K volume-matched finite-packet deposition layout intervention.
The only process intervention is the wet arrival point law. Conditional on the
fixed count, homogeneous Poisson locations are independent uniform positions.
Droplet profiles remain authored footprints; impact/fluid spreading and firing
are not solved. Existing linear capillary evolution and optical law are fixed.
"""
from pathlib import Path
import numpy as np,json,hashlib,resource,time,struct,zlib,gc,io
from scipy.fft import rfft2,irfft2,rfftfreq,fftfreq
from PIL import Image
R=Path(__file__).resolve().parents[1];BASE=R/'inputs/original_porcelain'
old=json.loads((BASE/'receipts/coating_state.json').read_text());P=old['parameters'];N=P['resolution'];L=P['tile_m'];assert N==4096;step_m=L/N;pixel_area=step_m**2;start=time.monotonic();rng=np.random.default_rng(P['seed']);packet=[]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def kernel(cx,cy,rx,ry,h):
 xx=np.arange(int(cx-rx*1.8),int(cx+rx*1.8+1));yy=np.arange(int(cy-ry*1.8),int(cy+ry*1.8+1));u=(xx[None,:]-cx)/rx;v=(yy[:,None]-cy)/ry;p=np.exp(-2*(u*u+v*v)).astype('f4')*h
 return xx,yy,p

def deposit(cells,amplitude):
 a=np.zeros((N,N),'f4');step=N/cells;rows=[]
 for j in range(cells):
  for i in range(cells):
   cx=(i+.5+rng.uniform(-.16,.16))*step;cy=(j+.5+rng.uniform(-.16,.16))*step;rx=step*rng.uniform(.52,.73);ry=rx*rng.uniform(.85,1.15);h=amplitude*rng.uniform(.7,1.25)
   xx,yy,k=kernel(cx,cy,rx,ry,h);a[np.ix_(yy%N,xx%N)]+=k;rows.append([cx,cy,rx,ry,h,float(k.sum(dtype='f8')*pixel_area)])
 return a,np.asarray(rows,'f8')

def chunk(f,n,d):f.write(struct.pack('>I',len(d))+n+d+struct.pack('>I',zlib.crc32(n+d)&0xffffffff))

def normal16(path,dx,dy):
 with path.open('xb') as f:
  f.write(b'\x89PNG\r\n\x1a\n');chunk(f,b'IHDR',struct.pack('>IIBBBBB',N,N,16,2,0,0,0));enc=zlib.compressobj(6)
  for j in range(N):
   a=np.stack((-dx[j],dy[j],np.ones(N,'f4')),axis=-1);a/=np.linalg.norm(a,axis=-1,keepdims=True);q=np.rint((a*.5+.5)*65535).astype('>u2');b=enc.compress(b'\0'+q.tobytes())
   if b:chunk(f,b'IDAT',b)
  chunk(f,b'IDAT',enc.flush());chunk(f,b'IEND',b'')
 # Independent CRC/deflate/integer decoding, compared to analytical row normals.
 data=path.read_bytes();pos=8;enc=zlib.decompressobj();buf=bytearray();row=0;stride=N*6+1
 while pos<len(data):
  le=struct.unpack('>I',data[pos:pos+4])[0];tag=data[pos+4:pos+8];b=data[pos+8:pos+8+le];assert zlib.crc32(tag+b)&0xffffffff==struct.unpack('>I',data[pos+8+le:pos+12+le])[0];pos+=le+12
  if tag==b'IDAT':buf.extend(enc.decompress(b))
  while len(buf)>=stride:
   assert buf[0]==0;q=np.frombuffer(buf[1:stride],dtype='>u2').reshape(N,3);a=np.stack((-dx[row],dy[row],np.ones(N,'f4')),axis=-1);a/=np.linalg.norm(a,axis=-1,keepdims=True);assert np.array_equal(q,np.rint((a*.5+.5)*65535).astype('u2'));del buf[:stride];row+=1
 assert row==N and not buf
 return {'path':str(path),'sha256':sha(path),'bytes':path.stat().st_size,'bits':16,'resolution':[N,N],'PNG_integer_roundtrip_exact':True,'colorspace':'Non-Color'}

def gray(name,a,bits):
 assert a.min()>=0 and a.max()<=1 and np.isfinite(a).all();p=R/'maps'/name;p.parent.mkdir(exist_ok=True);q=np.rint(a*(65535 if bits==16 else 255)).astype('u2' if bits==16 else 'u1');Image.fromarray(q).save(p,compress_level=4)
 return {'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size,'bits':bits,'resolution':[N,N],'colorspace':'Non-Color'}

def structure(points):
 xy=points[:,:2]/N;counts=np.histogram2d(xy[:,0],xy[:,1],bins=80,range=[[0,1],[0,1]])[0];bragg={}
 for k in [(80,0),(0,80),(80,80),(160,0),(0,160)]:bragg[str(k)]=float(abs(np.exp(2j*np.pi*(xy@np.array(k))).sum())**2/len(xy))
 return {'count':len(xy),'quadrat80_mean':float(counts.mean()),'quadrat80_variance':float(counts.var()),'conditional_Poisson_expected_variance':1-1/6400,'structure_factor':bragg,'nonzero_integer_mode_Poisson_expectation':1.0}

body,body_packets=deposit(160,P['powder_skin_amplitude_m']);wet,wet_packets=deposit(80,P['wet_deposit_amplitude_m']);body_mean=float(body.mean(dtype='f8'));wet_mean=float(wet.mean(dtype='f8'));volume=float(wet.sum(dtype='f8')*pixel_area)
body_map=gray('Substrate_Height.png',.5+(body-np.mean(body))/20e-6,16);assert body_map['sha256']==old['maps']['Substrate_Height']['sha256'],'Original pressed body replay differs'
kx=(2*np.pi*rfftfreq(N,L/N)).astype('f4');ky=(2*np.pi*fftfreq(N,L/N)).astype('f4');decay=np.exp(-P['integrated_capillary_mobility_m4']*(kx[None,:]**2+ky[:,None]**2)**2);del kx,ky
records=[]
def finish(name,wet_field):
 h0=body+wet_field+P['nominal_glaze_thickness_m'];initial_mean=float(h0.mean(dtype='f8'));initial_rms=float(np.std(h0));F=rfft2(h0);H=irfft2(F*decay,s=(N,N)).astype('f4');hot=irfft2(F*decay*decay,s=(N,N)).astype('f4');del F,h0
 err=float(abs(H.mean(dtype='f8')-initial_mean));rms=float(np.std(H));hot_rms=float(np.std(hot));del hot;assert hot_rms<rms<initial_rms and err<1e-10
 H-=np.mean(H);maps={};height_scale=max(40e-6,float(2*np.max(abs(H)))*1.000001);maps['Height']=gray(name+'/Glaze_Height.png',.5+H/height_scale,16)
 dx=(np.roll(H,-1,1)-np.roll(H,1,1))*N/(2*L);dy=(np.roll(H,-1,0)-np.roll(H,1,0))*N/(2*L);slope=np.sqrt(dx*dx+dy*dy);rough=np.clip(.17984575+.008*(slope/(np.mean(slope)+1e-12)-1),.17184575,.19584575).astype('f4');maps['Roughness']=gray(name+'/Glaze_Roughness.png',rough,8);del rough,slope
 normal_path=R/'maps'/name/'Glaze_Normal_OpenGL_RGB16.png';maps['Normal']=normal16(normal_path,dx,dy)
 if name=='control':
  assert maps['Height']['sha256']==old['maps']['Glaze_Height']['sha256'],'Original glaze-height replay differs'
  assert maps['Roughness']['sha256']==old['maps']['Glaze_Roughness']['sha256'],'Original roughness replay differs'
  samples=json.loads((BASE/'receipts/normal16_decode_samples.json').read_text())['samples'];errors=[]
  for p in samples:
   y,x=p['y_top'],p['x'];v=np.array([-dx[y,x],dy[y,x],1.],'f4');v/=np.linalg.norm(v);errors.append(int(np.max(abs(np.rint((v*.5+.5)*65535).astype('i8')-np.asarray(p['RGB16'],'i8')))))
  assert max(errors)==0,errors
  maps['Normal']['original_RGB16_decode_samples_exact']=len(samples)
 rec={'case':name,'initial_rms_m':initial_rms,'leveled_rms_m':rms,'double_mobility_rms_m':hot_rms,'mean_leveling_error_m':err,'wet_volume_m3':float(wet_field.sum(dtype='f8')*pixel_area),'height_scale_m':height_scale,'normal_and_roughness_share_same_H':True,'maps':maps}
 records.append(rec);(R/'receipts/progress.json').write_text(json.dumps(records,indent=2));del H,dx,dy;gc.collect();print('CASE_READY',name,rms,flush=True)
finish('control',wet);del wet;gc.collect()
# Same packet marks/volume, one independent spatial flux process with fixed count.
centers=np.random.default_rng(P['seed']).uniform(0,N,(len(wet_packets),2));new_packets=wet_packets.copy();new_packets[:,:2]=centers;newwet=np.zeros((N,N),'f4');per_packet=[]
for i,(cx,cy,rx,ry,h,target) in enumerate(new_packets):
 xx,yy,k=kernel(cx,cy,rx,ry,h);mass=float(k.sum(dtype='f8')*pixel_area);correction=target/mass;k=(k.astype('f8')*correction).astype('f4');actual=float(k.sum(dtype='f8')*pixel_area);per_packet.append(abs(actual-target)/target);newwet[np.ix_(yy%N,xx%N)]+=k
mass_error=float(abs(newwet.sum(dtype='f8')*pixel_area-volume)/volume);assert mass_error<1e-6 and max(per_packet)<1e-6
finish('conditional_poisson',newwet);del newwet
np.savez_compressed(R/'state/arrival_packets.npz',pressed_packets=body_packets,control_wet_packets=wet_packets,candidate_wet_packets=new_packets)
report={'parameters':P,'source_original':str(BASE/'source/coating_state.py'),'source_original_sha256':sha(BASE/'source/coating_state.py'),'original_pressed_height_and_glaze_height_and_roughness_bytes_exact':True,'intervention':'Uniform spatial Poisson deposition conditioned on the same6400 arrivals and unchanged packet marks/total delivered volume; no arbitrary height/noise field added','event_law_reference':'https://www.sciencedirect.com/science/article/pii/S0032591017309063','profile_limit':'Existing positive Gaussian footprints are retained authored deposition kernels. Impact, spreading, drying, firing and curved-edge flow remain unmodeled. Homogeneous flux is an assumed controlled-process idealization, not measured glaze statistics.','same_physical_parameters':['packet count/radii/profile/target volume','pressed substrate','linear capillary mobility','nominal film thickness','optical constants and roughness-response formula'],'packet_quadrature_mass_max_relative_error':max(per_packet),'total_wet_volume_relative_error':mass_error,'control_wet_volume_m3':volume,'retained_mean_film_m':P['nominal_glaze_thickness_m']+wet_mean,'control_statistics':structure(wet_packets),'candidate_statistics':structure(new_packets),'conditional_count_explanation':'Counts in disjoint cells are weakly anticorrelated after conditioning on a fixed total. This is not an unconditioned variable-dose Poisson sample.','cases':records,'shared_body_map':body_map,'elapsed_seconds':time.monotonic()-start,'peak_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'visual_acceptance':False};(R/'receipts/generation.json').write_text(json.dumps(report,indent=2));print('GENERATION_READY',json.dumps({k:report[k] for k in ['total_wet_volume_relative_error','retained_mean_film_m','elapsed_seconds','peak_rss_mib']}),flush=True)
