"""Small independent reference for the audited self-clearance condition.

This is an audit/test utility, not a production or interval-arithmetic certificate.
The theorem is exact; floating computations here use explicit engineering margins.
Input leaves must be in curve order and have valid continuous error/arclength bounds.
"""
from dataclasses import dataclass
from pathlib import Path
import importlib.util
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
R4 = ROOT / 'carpet-rounded-yarn-unit/revisions/r4_bounded_loop/source'
sys.path.insert(0, str(R4))
from certified_segments import segment_closest


@dataclass
class Leaves:
    a: np.ndarray
    b: np.ndarray
    error: np.ndarray
    arc_upper: np.ndarray


def projection_distance_lower(a, b, c, d):
    """Reuse the independently reviewed r4 projection-separation bound.

    Candidate closest-point error can weaken the bound but not invent geometric
    separation. Its engineering rounding allowance is not interval arithmetic.
    """
    distance, _, _, _, error, _ = segment_closest(a, b, c, d)
    return np.maximum(0., distance-error)


def assess_self(leaves, curvature_upper, radius, *, closed=False,
                numerical_distance_margin=1e-12):
    """Sufficient physical normal-disc injectivity; zero construction-gap target.

    C1,1 regularity and the supplied curvature/error/arclength bounds are assumed.
    Failed inequalities mean UNPROVEN, not an intersection witness.
    The all-pairs reference is intentionally small; production can use a complete
    midpoint-tree broad phase with horizon 2*r+e_i+e_j+half_i+half_j+margin.
    """
    a, b, error, upper = (np.asarray(x, dtype=float) for x in
                          (leaves.a, leaves.b, leaves.error, leaves.arc_upper))
    K, r = float(curvature_upper), float(radius)
    valid = (len(a)>0 and K>=0 and r>0 and np.isfinite(K) and
             np.all(np.isfinite(a)) and np.all(np.isfinite(b)) and
             np.all(np.isfinite(error)) and np.all(error>=0) and
             np.all(np.isfinite(upper)) and np.all(upper>0))
    if not valid:
        return {'pass':False, 'status':'unproven_invalid_inputs'}
    if not K*r < 1:
        return {'pass':False, 'status':'unproven_nonfold', 'kappa_radius':K*r}
    local = np.inf if K == 0 else np.pi/K
    # Self pairs within a leaf are not enumerated; this checks their whole span.
    # This extra condition also guarantees any consecutive pair is local.
    if np.any(upper > local/2):
        return {'pass':False, 'status':'unproven_refine_long_leaves'}
    i, j = np.triu_indices(len(a), 1)
    prefix = np.r_[0., np.cumsum(upper)]
    # Inflate accumulated-sum/subtraction error. Formal certification needs
    # outward-rounded intervals for these bounds and for pi/K as well.
    arc_pad = 128*np.finfo(float).eps*float(prefix[-1])*max(1, len(a))
    span = prefix[j+1]-prefix[i]+arc_pad
    if closed:
        reverse = prefix[-1]-(prefix[j]-prefix[i+1])+arc_pad
        span = np.minimum(span, reverse)
    nearby = span <= local*(1-128*np.finfo(float).eps)
    far_i, far_j = i[~nearby], j[~nearby]
    lower = projection_distance_lower(a[far_i],b[far_i],a[far_j],b[far_j])
    clearance = lower-error[far_i]-error[far_j]-2*r-numerical_distance_margin
    passed = bool(np.all(clearance > 0))
    bad = np.flatnonzero(clearance<=0)
    return {'pass':passed, 'status':'sufficient_condition_pass' if passed else
            'unproven_distant_rectangle', 'kappa_radius':K*r,
            'local_arc_cutoff_m':None if not np.isfinite(local) else local,
            'leaf_count':len(a), 'local_pair_count':int(nearby.sum()),
            'distant_pair_count':len(far_i),
            'minimum_tested_surface_clearance_lower_m':float(clearance.min())
            if len(clearance) else None,
            'first_unresolved_leaf_pair': [int(far_i[bad[0]]),int(far_j[bad[0]])]
            if len(bad) else None, 'self_clearance_target_m':0.,
            'regularity_and_input_bounds_assumed':True,
            'formal_outward_rounded_arithmetic':False}


def join(*pieces):
    return Leaves(*(np.concatenate([getattr(x, name) for x in pieces])
                    for name in ('a','b','error','arc_upper')))


def line(a, b, count=4):
    points = np.linspace(a,b,count+1)
    return Leaves(points[:-1], points[1:], np.zeros(count),
                  np.linalg.norm(np.diff(points,axis=0),axis=1))


def circular_arc(center, radius, theta0, theta1, count=64):
    theta = np.linspace(theta0, theta1, count+1)
    points = np.asarray(center)+radius*np.c_[np.cos(theta),np.sin(theta),
                                            np.zeros_like(theta)]
    angle = abs(theta1-theta0)/count
    assert angle <= np.pi
    # Exact real-arithmetic Hausdorff error and arclength for circular arcs.
    return Leaves(points[:-1],points[1:],
                  np.full(count, radius*(1-np.cos(angle/2))),
                  np.full(count, radius*angle))


def cubic_leaves(control, subdivisions=64):
    """de Casteljau dyadic leaves with control-polygon length upper bounds."""
    assert subdivisions > 0 and subdivisions & (subdivisions-1) == 0
    bez = np.asarray(control,dtype=float)[None]
    while len(bez)<subdivisions:
        first=(bez[:,:-1]+bez[:,1:])/2
        second=(first[:,:-1]+first[:,1:])/2
        mid=(second[:,0]+second[:,1])/2
        left=np.stack((bez[:,0],first[:,0],second[:,0],mid),axis=1)
        right=np.stack((mid,second[:,1],first[:,2],bez[:,3]),axis=1)
        # Keep parameter order; the production bounded_chords depth batches do not.
        bez=np.stack((left,right),axis=1).reshape(-1,4,3)
    chord=bez[:,3]-bez[:,0]
    error=np.maximum(np.linalg.norm(bez[:,1]-bez[:,0]-chord/3,axis=1),
                     np.linalg.norm(bez[:,2]-bez[:,0]-2*chord/3,axis=1))
    upper=np.linalg.norm(np.diff(bez,axis=1),axis=2).sum(axis=1)
    return Leaves(bez[:,0],bez[:,3],error,upper)
