"""Preserve signed pass roundoff outside both independently filtered components."""
import numpy as np
from two_component_math import tests as original_tests

def partition(t, r, delta):
    assert t.dtype == r.dtype == np.float32 and delta.dtype == np.float64
    assert t.shape == r.shape == delta.shape and np.isfinite(t).all() and np.isfinite(r).all() and np.isfinite(delta).all()
    assert np.all(r >= 0)
    beauty = t.astype('f8') + r.astype('f8') + delta
    scale = np.maximum(np.maximum(abs(beauty), abs(t.astype('f8')) + r), np.finfo(np.float32).tiny)
    assert np.all(abs(delta) <= scale * 2.**-16), 'Original physical residual exceeds guard'
    assert np.all(np.maximum(-t.astype('f8'), 0) <= scale * 2.**-16), 'Signed transmission exceeds existing guard'
    positive = np.where(t >= 0, t, np.float32(0))
    signed = np.where(t < 0, t, np.float32(0))
    assert np.array_equal(positive.astype('f8') + signed.astype('f8'), t.astype('f8'))
    return positive, signed

def combine(t, r, delta, signed):
    assert t.shape == r.shape == delta.shape == signed.shape
    assert delta.dtype == np.float64 and all(np.isfinite(a).all() for a in [t,r,delta,signed])
    assert np.all(t >= 0) and np.all(r >= 0) and np.all(signed <= 0)
    result = (t.astype('f8') + signed.astype('f8') + r.astype('f8') + delta).astype('f4')
    assert np.isfinite(result).all() and np.all(result >= 0), 'No clipping of negative final radiance'
    return result

def tests():
    report = original_tests()
    t = np.array([[[.3,-3e-9,.1]]],dtype='f4');r=np.array([[[.2,.4,.2]]],dtype='f4')
    beauty=(t.astype('f8')+r).astype('f4');delta=beauty.astype('f8')-t.astype('f8')-r
    positive,signed=partition(t,r,delta)
    assert np.array_equal(combine(positive,r,delta,signed),beauty)
    assert np.array_equal(signed[t<0],t[t<0])
    for bad in [np.array([[[-.1,.1,.1]]],dtype='f4'),np.full_like(t,np.nan)]:
        try:partition(bad,r,delta)
        except AssertionError:pass
        else:raise AssertionError('Invalid signed term admitted')
    report.update(signed_identity_exact=True,signed_term_preserved=True,excess_signed_term_rejected=True)
    return report
