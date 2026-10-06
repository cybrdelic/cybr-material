"""Regenerate original packets at native4096 and evolve one imposed reference history.
No Blender invocation, image render, upload, network operation or appearance fit.
"""
from pathlib import Path
import gc
import hashlib
import json
import resource
import struct
import sys
import time
import zlib
import numpy as np
from PIL import Image
from scipy.fft import rfft2, irfft2, rfftfreq, fftfreq
from model import HERE, ROOT, ARRIVALS, GAMMAS, R0, finite_depth, setup, sha

if '--admitted' not in sys.argv:
    raise SystemExit('Generation requires the separately admitted resource guard.')
start = time.monotonic()
report = setup()
original = json.loads((ARRIVALS/'receipts/generation.json').read_text())
p = original['parameters']
N, L = p['resolution'], p['tile_m']
assert N == 4096
step_m, pixel_area = L/N, (L/N)**2
maps = HERE/'maps'
state = HERE/'state'
assert not (HERE/'receipts/generation.json').exists(), 'Never overwrite a completed candidate'
rng = np.random.default_rng(p['seed'])
saved = np.load(ARRIVALS/'state/arrival_packets.npz')

def deposit(cells, amplitude, key):
    a = np.zeros((N, N), dtype=np.float32)
    spacing = N/cells
    rows = []
    for j in range(cells):
        for i in range(cells):
            cx=(i+.5+rng.uniform(-.16,.16))*spacing
            cy=(j+.5+rng.uniform(-.16,.16))*spacing
            rx=spacing*rng.uniform(.52,.73)
            ry=rx*rng.uniform(.85,1.15)
            h=amplitude*rng.uniform(.7,1.25)
            xx=np.arange(int(cx-rx*1.8),int(cx+rx*1.8+1))
            yy=np.arange(int(cy-ry*1.8),int(cy+ry*1.8+1))
            u=(xx[None,:]-cx)/rx
            v=(yy[:,None]-cy)/ry
            k=np.exp(-2*(u*u+v*v)).astype('f4')*h
            a[np.ix_(yy%N,xx%N)]+=k
            rows.append([cx,cy,rx,ry,h,float(k.sum(dtype='f8')*pixel_area)])
    records=np.asarray(rows, dtype=np.float64)
    assert np.array_equal(records, saved[key]), key+' packet replay mismatch'
    return a, {'count': len(records), 'saved_packet_array_exact': True,
               'raw_packet_array_sha256': hashlib.sha256(records.tobytes()).hexdigest(),
               'field_sha256': hashlib.sha256(memoryview(a)).hexdigest(),
               'field_min_m': float(a.min()), 'field_max_m': float(a.max()),
               'field_mean_m': float(a.mean(dtype='f8'))}

def stats(a):
    mean=float(a.mean(dtype='f8'))
    variance=0.
    for j in range(0,N,64):
        variance+=float(np.sum((a[j:j+64].astype('f8')-mean)**2, dtype='f8'))
    return {'min_m': float(a.min()), 'max_m': float(a.max()),
            'mean_m': mean, 'rms_about_mean_m': (variance/(N*N))**.5}

def derivatives(a):
    dx=np.empty((N,N), dtype='f4')
    dy=np.empty((N,N), dtype='f4')
    factor=N/(2*L)
    for j in range(0,N,64):
        rows=np.arange(j,min(j+64,N))
        slab=a[rows]
        dx[rows]=(np.roll(slab,-1,axis=1)-np.roll(slab,1,axis=1))*factor
        dy[rows]=(a[(rows+1)%N]-a[(rows-1)%N])*factor
    return dx,dy

def slope_stats(dx,dy):
    minimum,maximum,total,total2=1e99,0.,0.,0.
    for j in range(0,N,64):
        s=np.hypot(dx[j:j+64].astype('f8'),dy[j:j+64].astype('f8'))
        minimum=min(minimum,float(s.min()))
        maximum=max(maximum,float(s.max()))
        total+=float(s.sum())
        total2+=float(np.sum(s*s))
    return {'min':minimum,'max':maximum,'mean':total/(N*N),'rms':(total2/(N*N))**.5}

def old_roughness_stats(dx,dy,mean_slope,multiplier=1.):
    minimum,maximum,total,total2=1.,0.,0.,0.
    for j in range(0,N,64):
        s=multiplier*np.hypot(dx[j:j+64].astype('f8'),dy[j:j+64].astype('f8'))
        r=np.clip(R0+.008*(s/(multiplier*mean_slope+1e-12)-1),R0-.008,R0+.016)
        minimum=min(minimum,float(r.min()))
        maximum=max(maximum,float(r.max()))
        total+=float(r.sum())
        total2+=float(np.sum(r*r))
    mean=total/(N*N)
    return {'height_amplitude_multiplier':multiplier,'minimum':minimum,'maximum':maximum,
            'mean':mean,'std':max(0.,total2/(N*N)-mean*mean)**.5}

def chunk(f,tag,data):
    f.write(struct.pack('>I',len(data))+tag+data+struct.pack('>I',zlib.crc32(tag+data)&0xffffffff))

def write_png(path,row_bytes,channels):
    with path.open('xb') as f:
        f.write(b'\x89PNG\r\n\x1a\n')
        chunk(f,b'IHDR',struct.pack('>IIBBBBB',N,N,16,0 if channels==1 else 2,0,0,0))
        encoder=zlib.compressobj(6)
        for row in range(N):
            encoded=encoder.compress(b'\0'+row_bytes(row))
            if encoded: chunk(f,b'IDAT',encoded)
        chunk(f,b'IDAT',encoder.flush())
        chunk(f,b'IEND',b'')
    decoder=zlib.decompressobj()
    buf=bytearray()
    row=0
    stride=N*channels*2+1
    with path.open('rb') as f:
        assert f.read(8)==b'\x89PNG\r\n\x1a\n'
        while True:
            header=f.read(8)
            if not header: break
            length,tag=struct.unpack('>I4s',header)
            data=f.read(length)
            assert struct.unpack('>I',f.read(4))[0]==zlib.crc32(tag+data)&0xffffffff
            if tag==b'IDAT': buf.extend(decoder.decompress(data))
            while len(buf)>=stride:
                assert buf[0]==0 and buf[1:stride]==row_bytes(row)
                del buf[:stride]
                row+=1
    assert row==N and not buf and decoder.eof
    return {'path':str(path),'sha256':sha(path),'bytes':path.stat().st_size,
            'resolution':[N,N],'channels':channels,'bits':16,'colorspace':'Non-Color',
            'CRC_and_integer_roundtrip_exact':True}

body,body_rec=deposit(160,p['powder_skin_amplitude_m'],'pressed_packets')
wet,wet_rec=deposit(80,p['wet_deposit_amplitude_m'],'control_wet_packets')
assert abs(wet_rec['field_mean_m']+p['nominal_glaze_thickness_m']-report['retained_mean_film_m'])<1e-18
wet_volume=float(wet.sum(dtype='f8')*pixel_area)
assert wet_volume==original['control_wet_volume_m3']
# Check the original native substrate encoding, without saving a second copy.
retained_body=np.asarray(Image.open(ARRIVALS/'maps/Substrate_Height.png'))
body_mean_old=np.mean(body)
for j in range(0,N,64):
    body_encoding=np.rint((.5+(body[j:j+64]-body_mean_old)/20e-6)*65535).astype('u2')
    assert np.array_equal(body_encoding,retained_body[j:j+64])
del retained_body, body_encoding
body_mean=body_rec['field_mean_m']
film_mean=p['nominal_glaze_thickness_m']+wet_rec['field_mean_m']
top_mean=body_mean+film_mean
# Same sampled body/wet fields. Float64 addition avoids the old f4 sum's ~1.8 pm DC bias.
initial=body.astype('f8')
initial+=wet
initial+=p['nominal_glaze_thickness_m']
initial_stats=stats(initial)
initial_film={'min_m':float(wet.min())+p['nominal_glaze_thickness_m'],
              'max_m':float(wet.max())+p['nominal_glaze_thickness_m'],
              'mean_m':film_mean,'volume_m3':film_mean*L*L}
old_rms=original['cases'][0]['initial_rms_m']
assert abs(initial_stats['rms_about_mean_m']/old_rms-1)<2e-6
initial-=top_mean
dx,dy=derivatives(initial)
initial_slope=slope_stats(dx,dy)
del dx,dy,wet
F=rfft2(initial,workers=1)
zero_initial=complex(F[0,0])
F[0,0]=0.+0.j  # Mathematical zero-mean field; total physical mean carried separately.
initial_stats['floating_residual_removed_m']=zero_initial.real/(N*N)
assert abs(initial_stats['floating_residual_removed_m'])<1e-18
del initial
gc.collect()
# Store k*F(kh)/2 in float64. In-place operations bound peak memory.
kx=2*np.pi*rfftfreq(N,step_m)
ky=2*np.pi*fftfreq(N,step_m)
rate=np.hypot(ky[:,None],kx[None,:])
for j in range(0,N,64):
    rate[j:j+64]*=finite_depth(rate[j:j+64]*film_mean)/2
assert rate[0,0]==0
J=report['history']['J_per_Pa']
responses=[]
map_records={}
central_raw=None
for gamma in GAMMAS:
    transfer=np.exp(-gamma*J*rate)
    transfer[0,0]=1.
    assert float(transfer.min())>=0 and float(transfer.max())==1.
    transformed=F*transfer
    assert transformed[0,0]==F[0,0]
    del transfer
    H=irfft2(transformed,s=(N,N),workers=1)
    del transformed
    hstat=stats(H)
    film_min,film_max,film_sum=1.,0.,0.
    for j in range(0,N,64):
        film=top_mean+H[j:j+64]-body[j:j+64]
        film_min=min(film_min,float(film.min()))
        film_max=max(film_max,float(film.max()))
        film_sum+=float(film.sum(dtype='f8'))
    film_actual=film_sum/(N*N)
    del film
    assert film_min>0 and abs(film_actual-film_mean)<1e-17
    dx,dy=derivatives(H)
    slope=slope_stats(dx,dy)
    response={'gamma_N_per_m':gamma,'exposure_m':gamma*J,'top_surface_deviation':hstat,
              'top_surface_absolute_range_m':[top_mean+hstat['min_m'],top_mean+hstat['max_m']],
              'film_min_m':film_min,'film_max_m':film_max,'film_mean_m':film_actual,
              'film_volume_m3':film_actual*L*L,'volume_relative_error':abs(film_actual/film_mean-1),
              'slope':slope,'rms_survival':hstat['rms_about_mean_m']/initial_stats['rms_about_mean_m'],
              'old_roughness_if_recomputed':old_roughness_stats(dx,dy,slope['mean'])}
    assert hstat['rms_about_mean_m']<initial_stats['rms_about_mean_m']
    if gamma==.3:
        scale=4e-5 # Retained physical encoding range, never normalization to candidate extrema.
        assert np.max(np.abs(H))<scale/2
        raw_path=state/'top_surface_deviation_m.npy'
        raw=H.astype('f4')
        cast_error=max(float(np.max(np.abs(raw[j:j+64].astype('f8')-H[j:j+64]))) for j in range(0,N,64))
        np.save(raw_path,raw,allow_pickle=False)
        del raw
        central_raw={'path':str(raw_path),'sha256':sha(raw_path),'bytes':raw_path.stat().st_size,
                     'dtype':'float32','resolution':[N,N], 'units':'m',
                     'top_surface_mean_m':top_mean,'max_float32_storage_error_m':cast_error}
        height_quant_error=0.
        normal_max_error=0.
        for j in range(0,N,64):
            decoded=(np.rint((.5+H[j:j+64]/scale)*65535)/65535-.5)*scale
            height_quant_error=max(height_quant_error,float(np.max(np.abs(decoded-H[j:j+64]))))
        def height_row(j):
            return np.rint((.5+H[j]/scale)*65535).astype('>u2').tobytes()
        def normal_row(j):
            a=np.stack((-dx[j].astype('f8'),dy[j].astype('f8'),np.ones(N)),axis=-1)
            a/=np.linalg.norm(a,axis=-1,keepdims=True)
            return np.rint((a*.5+.5)*65535).astype('>u2').tobytes()
        map_records['Height']=write_png(maps/'Glaze_Height.png',height_row,1)
        map_records['Height'].update({'scale_m':scale,'midlevel':.5,
                                     'quantization_max_abs_error_m':height_quant_error})
        map_records['Normal']=write_png(maps/'Glaze_Normal_OpenGL_RGB16.png',normal_row,3)
        map_records['Normal'].update({'convention':'OpenGL; (-dH/dx,+dH/dy,1), top-origin rows',
                                     'slope_method':'Periodic centered native-grid derivative in SI',
                                     'component_quantization_bound':1/65535})
        response['old_normalized_slope_amplitude_audit']=[old_roughness_stats(dx,dy,slope['mean'],m)
                                                        for m in (1.,1e-3,1e-6,1e-9,0.)]
    responses.append(response)
    (HERE/'receipts/progress.json').write_text(json.dumps(responses,indent=2,allow_nan=False)+'\n')
    print('NUMERIC_CASE_READY',gamma,json.dumps({'rms_m':hstat['rms_about_mean_m'],
          'max_slope':slope['max'],'film_min_m':film_min,'volume_relative_error':response['volume_relative_error']}),flush=True)
    del H,dx,dy
    gc.collect()
assert responses[0]['top_surface_deviation']['rms_about_mean_m']>responses[1]['top_surface_deviation']['rms_about_mean_m']>responses[2]['top_surface_deviation']['rms_about_mean_m']
report.update({'initial_replay':{'body':body_rec,'wet':wet_rec,'packet_marks_and_positions_exact':True,
               'original_substrate_map_pixels_exact':True,'wet_volume_m3':wet_volume,
               'initial_film':initial_film,'initial_top_surface':initial_stats,'initial_slope':initial_slope,
               'construction':'Exact f4 sampled original-control body and wet fields, combined/evolved in f8; no old leveled map used as initial state'},
               'responses':responses,'central_state':central_raw,'maps':map_records,
               'physical_encoding':{'height_scale_m':4e-5,'tile_m':L,'mean_film_m':film_mean,
                                    'mean_top_m':top_mean,'substrate_mean_m':body_mean},
               'first_binding_scope':{'history_variant':'Evolved height and normal only',
                                      'coat_roughness':'Exact existing control Glaze_Roughness image, retained for isolated attribution',
                                      'roughness_source_path':original['cases'][0]['maps']['Roughness']['path'],
                                      'roughness_source_sha256':original['cases'][0]['maps']['Roughness']['sha256'],
                                      'shared_state_PBR_qualification':False},
               'held_optical_variant':{'status':'HELD; not bound or rendered',
                                       'coat_roughness_scalar':R0,
                                       'purpose':'Separate intrinsic-only roughness accounting, controlled against the new-history frame'},
               'all_physical_fields_finite':True,'all_sensitivity_films_positive':True,
               'mean_and_volume_conservation_pass':True,
               'elapsed_seconds':time.monotonic()-start,
               'peak_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
               'source_sha256':{q.name:sha(q) for q in (HERE/'source').glob('*.py')},
               'rendered':False,'selected':False,'visual_acceptance':False})
for rel,expected in report['input_sha256'].items():
    assert sha(ROOT/rel)==expected, 'Input changed: '+rel
(HERE/'receipts/generation.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
print('GENERATION_READY',json.dumps({'elapsed_seconds':report['elapsed_seconds'],
      'peak_rss_mib':report['peak_rss_mib'],'maps':map_records}),flush=True)
