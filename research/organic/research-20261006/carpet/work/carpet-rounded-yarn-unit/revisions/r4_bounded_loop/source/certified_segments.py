"""Projection-interval lower bound independent of closest-point conditioning."""
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];R3=ROOT.parent/'r3_swept_feasibility'
sys.path.insert(0,str(R3/'source'))
from swept_feasibility import segment_closest as candidate_closest

def segment_closest(a,b,c,d):
    distance,s,t,normal,_,feature=candidate_closest(a,b,c,d)
    # Translate first: this reduces cancellation for far-from-origin scenes.
    p0=np.zeros_like(distance);p1=np.sum((b-a)*normal,axis=-1);q0=np.sum((c-a)*normal,axis=-1);q1=np.sum((d-a)*normal,axis=-1)
    pmin=np.minimum(p0,p1);pmax=np.maximum(p0,p1);qmin=np.minimum(q0,q1);qmax=np.maximum(q0,q1)
    gap=np.maximum.reduce((np.zeros_like(distance),qmin-pmax,pmin-qmax))
    numerical=128*np.finfo(float).eps*(np.linalg.norm(b-a,axis=-1)+np.linalg.norm(c-a,axis=-1)+np.linalg.norm(d-a,axis=-1))
    lower=np.maximum(gap-numerical,0);error=np.maximum(distance-lower,0)
    return distance,s,t,normal,error,feature
