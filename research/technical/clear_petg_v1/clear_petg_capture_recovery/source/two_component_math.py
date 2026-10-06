"""Explicit two-component filtering accounting. Never clips radiance."""
import numpy as np
def combine(filtered_transmission,filtered_remainder,delta64):
    t=np.asarray(filtered_transmission);r=np.asarray(filtered_remainder);d=np.asarray(delta64)
    assert t.shape==r.shape==d.shape and t.ndim==3 and t.shape[2]==3
    assert d.dtype==np.float64 and np.isfinite(t).all() and np.isfinite(r).all() and np.isfinite(d).all()
    assert np.all(t>=0) and np.all(r>=0),'Filtered physical radiance must be nonnegative; no clipping'
    out=(t.astype(np.float64)+r.astype(np.float64)+d).astype(np.float32)
    assert np.isfinite(out).all() and np.all(out>=0),'Final radiance must be nonnegative; no clipping'
    return out

def tests():
    rng=np.random.default_rng(2306);t=rng.uniform(.001,2,(13,17,3)).astype('f4');r=rng.uniform(.001,2,t.shape).astype('f4');b=(t.astype('f8')+r+rng.uniform(-1e-7,1e-7,t.shape)).astype('f4');delta=b.astype('f8')-t.astype('f8')-r.astype('f8')
    assert np.array_equal(combine(t,r,delta),b)
    assert np.array_equal(combine(t,r,delta),combine(r,t,delta))
    for bad in [-np.ones_like(t),np.full_like(t,np.nan)]:
        try:combine(bad,r,delta)
        except AssertionError:pass
        else:raise AssertionError('Invalid filtered radiance admitted')
    try:combine(np.zeros_like(t),np.zeros_like(r),-np.ones(t.shape,'f8'))
    except AssertionError:pass
    else:raise AssertionError('Negative final radiance admitted')
    return {'identity_float32_exact':True,'component_exchange_symmetric':True,'negative_filter_rejected':True,'nonfinite_filter_rejected':True,'negative_final_rejected':True,'clipping':False}
