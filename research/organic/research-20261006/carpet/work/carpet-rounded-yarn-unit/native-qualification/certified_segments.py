"""Self-contained r4 projection-interval segment lower bound.

The approximate closest-point calculation only chooses a projection direction.
It is never used as a distance certificate. This reproduces the conservative
projection construction in revisions/r4_bounded_loop/source/certified_segments.py
without importing its solver dependency tree. Numerical pads are engineering
allowances, not formal outward-rounded interval arithmetic.
"""
import numpy as np


def projection_distance_lower(a, b, c, d):
    """Return a continuous lower bound for batches of finite chord pairs."""
    a, b, c, d = (np.asarray(x, dtype=np.float64) for x in (a, b, c, d))
    u, v, w = b-a, d-c, a-c
    aa = np.sum(u*u, axis=-1)
    bb = np.sum(u*v, axis=-1)
    cc = np.sum(v*v, axis=-1)
    dd = np.sum(u*w, axis=-1)
    ee = np.sum(v*w, axis=-1)
    determinant = aa*cc-bb*bb
    product = aa*cc
    relative = np.divide(determinant, product, out=np.zeros_like(product), where=product>0)
    regular = relative > 64*np.finfo(float).eps
    ss, tt = [], []
    for s in (0., 1.):
        ss.append(np.full_like(aa, s))
        tt.append(np.clip(np.divide(ee+s*bb, cc, out=np.zeros_like(cc), where=cc>0), 0, 1))
    for t in (0., 1.):
        ss.append(np.clip(np.divide(t*bb-dd, aa, out=np.zeros_like(aa), where=aa>0), 0, 1))
        tt.append(np.full_like(cc, t))
    si = np.divide(bb*ee-cc*dd, determinant, out=np.zeros_like(aa), where=regular)
    ti = np.divide(aa*ee-bb*dd, determinant, out=np.zeros_like(aa), where=regular)
    ss.append(si)
    tt.append(ti)
    s, t = np.stack(ss), np.stack(tt)
    delta = w[None]+s[..., None]*u[None]-t[..., None]*v[None]
    square = np.sum(delta*delta, axis=-1)
    square[4] = np.where(regular & (si>=0) & (si<=1) & (ti>=0) & (ti<=1), square[4], np.inf)
    best = np.argmin(square, axis=0)
    selected = delta.reshape(5, -1, 3)[best.ravel(), np.arange(aa.size)].reshape(aa.shape+(3,))
    norm = np.linalg.norm(selected, axis=-1)
    normal = np.divide(selected, norm[..., None], out=np.zeros_like(selected), where=norm[..., None]>0)
    # Shrink slightly so roundoff in normalization cannot create a norm > 1.
    normal *= 1-32*np.finfo(float).eps
    zero = np.zeros_like(aa)
    p1 = np.sum(u*normal, axis=-1)
    q0 = np.sum((c-a)*normal, axis=-1)
    q1 = np.sum((d-a)*normal, axis=-1)
    gap = np.maximum.reduce((zero, np.minimum(q0, q1)-np.maximum(zero, p1),
                             np.minimum(zero, p1)-np.maximum(q0, q1)))
    pad = 128*np.finfo(float).eps*(np.linalg.norm(u, axis=-1)+
                                 np.linalg.norm(c-a, axis=-1)+np.linalg.norm(d-a, axis=-1))
    return np.maximum(gap-pad, 0.)
