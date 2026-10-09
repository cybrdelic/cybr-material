"""Reproducible numerical audit of frozen mesh, chart and native height/normal maps."""
import hashlib, json, math, struct, sys, zlib
from pathlib import Path
import numpy as np
from PIL import Image

R = Path(__file__).resolve().parent
E = Path(sys.argv[1]).resolve()
D = json.loads((R / 'scene_data.json').read_text())
M = D['objects']['07 new / finite60micron mean glaze shell']
V = np.asarray(M['vertices'])
LO = M['loops']
F = M['polygons']
period, edge, ix, zc = .0192, .0006, .0474, .0012
outer = [p for p in F if p['normal'][2] > .02]
n = len(V) // 2
# Solidify pairs each original vertex with the corresponding offset vertex.
assert all(LO[k]['vertex'] < n for p in outer for k in p['loops'])
assert np.max(np.linalg.norm(V[:n] - V[n:], axis=1)) < .00007

def unit(a):
    a = np.asarray(a, dtype=float)
    return a / np.linalg.norm(a, axis=-1, keepdims=True)

def angle(a, b):
    return float(np.degrees(np.arccos(np.clip(np.dot(unit(a), unit(b)), -1, 1))))

def stats(x):
    x = np.asarray(x, dtype=float)
    return {'count': len(x), 'min': float(x.min()), 'median': float(np.median(x)),
            'p95': float(np.quantile(x, .95)), 'max': float(x.max())}

def nominal_normal(p):
    # Inner positions are the original 600 um bevel. Use those rather than
    # interpreting numerical Solidify offset tangential drift as curvature.
    q = p - np.array([np.clip(p[0], -ix, ix), np.clip(p[1], -ix, ix), zc])
    q[np.abs(q) < 1e-8] = 0
    return unit(q)

expected = np.asarray([nominal_normal(p) for p in V[n:]])
errors, end_errors, edge_errors, corner_errors, top_errors = [], [], [], [], []
examples = []
for p in outer:
    for k in p['loops']:
        vid = LO[k]['vertex']
        error = angle(LO[k]['normal'], expected[vid])
        errors.append(error)
        ref = expected[vid]
        if abs(ref[2]) < 1e-5:
            end_errors.append(error)
        elif np.count_nonzero(np.abs(ref[:2]) > 1e-5) == 1:
            edge_errors.append(error)
        elif np.count_nonzero(np.abs(ref[:2]) > 1e-5) == 2:
            corner_errors.append(error)
        else:
            top_errors.append(error)
        if p['index'] in [0, 209, 216] and k == p['loops'][0]:
            examples.append({'face': p['index'], 'loop': k, 'vertex': vid,
                             'inner_position_m': V[vid+n].tolist(),
                             'shading_normal': LO[k]['normal'],
                             'nominal_normal': ref.tolist(), 'error_degrees': error})

edges = {}
for p in outer:
    ll = p['loops']
    for a, b in zip(ll, ll[1:] + ll[:1]):
        va, vb = LO[a]['vertex'], LO[b]['vertex']
        edges.setdefault(tuple(sorted([va, vb])), []).append((p['index'], {va:a, vb:b}))
seams = []
normal_jumps, tangent_jumps, uv_jumps = [], [], []
for key, pair in edges.items():
    if len(pair) != 2:
        continue
    for vertex in key:
        a, b = pair[0][1][vertex], pair[1][1][vertex]
        ne, te = angle(LO[a]['normal'], LO[b]['normal']), angle(LO[a]['tangent'], LO[b]['tangent'])
        delta = np.asarray(LO[a]['uv']) - LO[b]['uv']
        uv_jumps.append(float(np.max(np.abs(delta - np.rint(delta)))))
        normal_jumps.append(ne)
        tangent_jumps.append(te)
        if ne > .1:
            seams.append({'faces': [pair[0][0], pair[1][0]], 'vertex': vertex,
                          'normal_jump_degrees': ne, 'tangent_jump_degrees': te})

metrics = {'flat': [], 'straight': [], 'corner': []}
for t in M['triangles']:
    p = F[t['polygon']]
    if p['normal'][2] <= .02:
        continue
    ll = t['loops']
    positions = V[[LO[k]['vertex'] for k in ll]]
    uv = np.asarray([LO[k]['uv'] for k in ll]) * period
    jacobian = (positions[1:] - positions[0]).T @ np.linalg.inv((uv[1:] - uv[0]).T)
    s = np.linalg.svd(jacobian, compute_uv=False)
    kind = 'flat' if p['normal'][2] > .999 else ('straight' if min(abs(p['normal'][0]), abs(p['normal'][1])) < 1e-6 else 'corner')
    metrics[kind].extend(s.tolist())

# Face-center deviation compares smooth interpolation with the nominal parent
# cylinder/sphere. The long strips have no intermediate axial vertices.
center_errors = []
for p in outer:
    ids = [LO[k]['vertex'] for k in p['loops']]
    if np.max(np.ptp(V[ids, :2], axis=0)) < .09 or p['normal'][2] > .999:
        continue
    center_errors.append({'face': p['index'],
        'nominal_center_normal': nominal_normal(np.mean(V[np.asarray(ids)+n], axis=0)).tolist(),
        'smooth_center_normal': unit(np.mean([LO[k]['normal'] for k in p['loops']], axis=0)).tolist(),
        'error_degrees': angle(np.mean([LO[k]['normal'] for k in p['loops']], axis=0),
                               nominal_normal(np.mean(V[np.asarray(ids)+n], axis=0)))})

def normal_samples(path):
    samples = []
    with path.open('rb') as f:
        assert f.read(8) == b'\x89PNG\r\n\x1a\n'
        dec, buf, row = zlib.decompressobj(), bytearray(), 0
        while True:
            size = struct.unpack('>I', f.read(4))[0]
            tag, data = f.read(4), f.read(size)
            crc = struct.unpack('>I', f.read(4))[0]
            assert zlib.crc32(tag + data) & 0xffffffff == crc
            if tag == b'IHDR':
                width, height, bits, mode, *_ = struct.unpack('>IIBBBBB', data)
                assert (width, height, bits, mode) == (4096, 4096, 16, 2)
                stride = width * 6 + 1
            if tag == b'IDAT':
                buf.extend(dec.decompress(data))
                while len(buf) >= stride:
                    assert buf[0] == 0
                    if row % 16 == 0:
                        samples.append(np.frombuffer(buf[1:stride], dtype='>u2').reshape(width, 3)[::16].copy())
                    del buf[:stride]
                    row += 1
            if tag == b'IEND':
                break
        assert row == 4096 and not buf
    return np.asarray(samples, dtype=float).reshape(-1, 3) / 65535 * 2 - 1

generation = json.loads((E/'porcelain_deposition_arrivals/receipts/generation.json').read_text())
fields = {}
for case in generation['cases']:
    name = case['case']
    base = E/'porcelain_deposition_arrivals/maps'/name
    sample = normal_samples(base/'Glaze_Normal_OpenGL_RGB16.png')
    q = np.asarray(Image.open(base/'Glaze_Height.png'))
    y, x = np.meshgrid(np.arange(0, 4096, 16), np.arange(0, 4096, 16), indexing='ij')
    factor = case['height_scale_m'] / 65535 / (2*period/4096)
    dx = (q[y, (x+1)%4096].astype(float) - q[y, (x-1)%4096]) * factor
    dy = (q[(y+1)%4096, x].astype(float) - q[(y-1)%4096, x]) * factor
    height_normal = unit(np.stack([-dx, dy, np.ones_like(dx)], axis=-1).reshape(-1, 3))
    discrepancies = np.degrees(np.arccos(np.clip(np.sum(unit(sample)*height_normal, axis=1), -1, 1)))
    inclinations = np.degrees(np.arctan2(np.linalg.norm(sample[:, :2], axis=1), sample[:, 2]))
    fields[name] = {'sample_count': len(sample), 'sample_stride': 16,
        'normal_tilt_degrees': stats(inclinations),
        'decoded_normal_vs_quantized_height_degrees': stats(discrepancies),
        'height_scale_m': case['height_scale_m'],
        'normal_sha256': hashlib.sha256((base/'Glaze_Normal_OpenGL_RGB16.png').read_bytes()).hexdigest()}
    del q

result = {'scene_sha256': D['sha256'], 'outer_face_count': len(outer),
    'macro_loop_error_vs_nominal_parent_degrees': {
        'all': stats(errors), 'termination_boundary': stats(end_errors),
        'cylindrical_interior': stats(edge_errors), 'corner_interior': stats(corner_errors),
        'top_boundary_and_flat': stats(top_errors)},
    'macro_examples': examples, 'normal_discontinuities': seams,
    'maximum_outer_uv_phase_jump': max(uv_jumps),
    'maximum_outer_normal_jump_degrees': max(normal_jumps),
    'maximum_outer_tangent_jump_degrees': max(tangent_jumps),
    'chart_metric_surface_length_per_chart_length': {k: stats(v) for k, v in metrics.items()},
    'straight_strip_center_normal_errors': center_errors,
    'native_fields': fields,
    'source_files_unmodified': True,
    'scope': 'No render. Geometry and atlas are frozen. Nominal macro normal is evaluated on paired original inner bevel; corner normals assume the intended rounded-box sphere.'}
(R/'measurements.json').write_text(json.dumps(result, indent=2))
print(json.dumps({k:v for k,v in result.items() if k not in ['straight_strip_center_normal_errors','normal_discontinuities']}, indent=2))
