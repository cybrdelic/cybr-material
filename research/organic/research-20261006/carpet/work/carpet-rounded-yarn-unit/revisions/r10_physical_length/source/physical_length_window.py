"""Selective physical-length window with original disks and complete histories.

The optimization uses straight body segments. Only the separate native gate can
qualify the actual stored Catmull-Rom geometry. Same-fibre nonlocal clearance is
also reserved for that gate. Every frozen fixed body segment stays in the search.
"""
from pathlib import Path
import sys
import math
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.integrate import quad

ROOT=Path(__file__).resolve().parents[1]
UNIT=ROOT.parents[1]
WORK=UNIT.parent
sys.path.insert(0,str(WORK))
sys.path.insert(0,str(UNIT/'revisions/r6_c2_bridge/source'))
from bridge_field import BridgeField
# Keep the solver distance implementation distinct from the native gate module.
from advance_loop import segment_closest

DISK_M=30e-6
RESERVE_M=1.95e-6
CONSTRUCTION_TARGET_M=2e-6
ROUNDING_INSET_M=1e-9
BATCH=8192
MAX_PAIRS=1_000_000

def make_mask(stations,contact_fibres,physical_width,physical_stations=None):
    stations=np.asarray(stations,dtype=float)
    assert len(stations)>=4 and np.all(np.diff(stations)>0) and physical_width>0
    # Anchor is at/before end-width: subdivisions cannot shorten the physical
    # window. Root (index0) remains fixed even if the requested width reaches it.
    metric=stations if physical_stations is None else np.asarray(physical_stations)
    assert metric.shape==stations.shape and np.all(np.diff(metric)>0)
    anchor=max(0,int(np.searchsorted(metric,metric[-1]-physical_width,side='right')-1))
    mask=np.zeros((len(stations),187),dtype=bool)
    mask[anchor+1:,np.asarray(contact_fibres,dtype=int)]=True
    mask[-2:]=True
    mask[0]=False
    return anchor,mask

def pair_candidates(a,b,affected,owner,move,threshold):
    """Complete midpoint-sphere search at original disk centers, then pruning.

Each segment shifts by at most the largest allowed displacement of its two
endpoints. Thus minimum distance can fall by at most move_i+move_j, and a pair
excluded at those centers can never reach the threshold anywhere in the disks.
No k-nearest-neighbor cap is used. Same-owner pairs go to the native gate later.
"""
    middle=(a+b)*.5;half=np.linalg.norm(b-a,axis=1)*.5
    tree=cKDTree(middle);ii=[];jj=[];total=0
    ids=np.flatnonzero(affected)
    for start in range(0,len(ids),16):
        selected=ids[start:start+16]
        radii=threshold+half[selected]+half.max()+move[selected]+move.max()+1e-11
        counts=tree.query_ball_point(middle[selected],radii,workers=1,return_length=True)
        if total+int(counts.sum())>MAX_PAIRS:raise RuntimeError('complete pair budget exceeded')
        lists=tree.query_ball_point(middle[selected],radii,workers=1)
        i=np.repeat(selected,[len(x) for x in lists])
        j=np.concatenate([np.asarray(x,dtype=int) for x in lists])
        valid=(owner[i]!=owner[j]) & (~affected[j] | (i<j))
        i=i[valid];j=j[valid];total+=len(i)
        for first in range(0,len(i),BATCH):
            sl=slice(first,first+BATCH)
            d,_,_,_,e,_=segment_closest(a[i[sl]],b[i[sl]],a[j[sl]],b[j[sl]])
            keep=d-e<=threshold+move[i[sl]]+move[j[sl]]+1e-11
            ii.append(i[sl][keep]);jj.append(j[sl][keep])
    return np.concatenate(ii),np.concatenate(jj)

def wrapper_candidates(a,b,affected,move,wa,wb,error,threshold):
    middle=(a+b)*.5;half=np.linalg.norm(b-a,axis=1)*.5
    wm=(wa+wb)*.5;wh=np.linalg.norm(wb-wa,axis=1)*.5
    tree=cKDTree(wm);ii=[];jj=[];total=0
    ids=np.flatnonzero(affected)
    for start in range(0,len(ids),16):
        selected=ids[start:start+16]
        radii=threshold+half[selected]+wh.max()+move[selected]+error.max()+1e-11
        counts=tree.query_ball_point(middle[selected],radii,workers=1,return_length=True)
        if total+int(counts.sum())>MAX_PAIRS:raise RuntimeError('complete wrapper pair budget exceeded')
        lists=tree.query_ball_point(middle[selected],radii,workers=1)
        i=np.repeat(selected,[len(x) for x in lists]);j=np.concatenate([np.asarray(x,dtype=int) for x in lists]);total+=len(i)
        for first in range(0,len(i),BATCH):
            sl=slice(first,first+BATCH)
            d,_,_,_,e,_=segment_closest(a[i[sl]],b[i[sl]],wa[j[sl]],wb[j[sl]])
            keep=d-e<=threshold+move[i[sl]]+error[j[sl]]+1e-11
            ii.append(i[sl][keep]);jj.append(j[sl][keep])
    return np.concatenate(ii),np.concatenate(jj)

def guide_pieces(field,s0,s1):
    """Exact polynomial pieces of the original cubic / local C2 quintic guide."""
    split=[s0,s1]
    original=field.original_spline
    split.extend(original.x[(original.x>s0)&(original.x<s1)].tolist())
    split.extend(v for lo,hi,_ in field.spline.parts for v in (lo,hi) if s0<v<s1)
    split=np.unique(split)
    for lo,hi in zip(split[:-1],split[1:]):
        midpoint=(lo+hi)*.5
        bridge=next((part for part in field.spline.parts if part[0]<=midpoint<=part[1]),None)
        if bridge:
            origin,end,power=bridge;ta=(lo-origin)/(end-origin);h=(hi-lo)/(end-origin)
        else:
            index=np.searchsorted(original.x,midpoint,side='right')-1
            power=original.c[:,index][::-1];ta=lo-original.x[index];h=hi-lo
        restricted=np.zeros_like(power)
        degree=len(power)-1
        for k in range(degree+1):
            for j in range(k,degree+1):restricted[k]+=power[j]*math.comb(j,k)*ta**(j-k)*h**k
        controls=np.zeros_like(power)
        for i in range(degree+1):
            for k in range(i+1):controls[i]+=restricted[k]*math.comb(i,k)/math.comb(degree,k)
        fractions=(lo-s0+(hi-lo)*np.arange(degree+1)/degree)/(s1-s0)
        yield fractions,controls

class PhysicalLengthWindow:
    def __init__(self,field,stations,xy,world,references,contact_fibres,physical_width,wrapper_leaves):
        self.field=field;self.stations=np.asarray(stations);self.original_world=np.asarray(world).copy()
        self.original_xy=np.asarray(xy).copy();self.references=np.asarray(references).copy()
        # The C2 repair retains its old parameter; that parameter is not exactly
        # physical arclength. Integrate speed before choosing the anchor.
        arc=[];errors=[]
        for s0,s1 in zip(self.stations[:-1],self.stations[1:]):
            split=field.original_spline.x[(field.original_spline.x>s0)&(field.original_spline.x<s1)]
            points=sorted(set(split.tolist()+[v for lo,hi,_ in field.spline.parts for v in (lo,hi) if s0<v<s1]))
            value,error=quad(lambda s:float(np.linalg.norm(field.spline(s,1))),s0,s1,
                             points=points,epsabs=1e-13,epsrel=1e-9,limit=max(100,len(points)+20))
            arc.append(value);errors.append(error)
        self.physical_stations=np.r_[0,np.cumsum(arc)]
        self.arc_integration_error_m=float(sum(errors))
        self.anchor,self.mask=make_mask(stations,contact_fibres,physical_width,self.physical_stations)
        self.physical_window_width_m=float(self.physical_stations[-1]-self.physical_stations[self.anchor])
        assert self.physical_window_width_m-self.arc_integration_error_m>=physical_width
        self.node_station,self.node_owner=np.where(self.mask);self.n=len(self.node_station)
        self.mapping=np.full(self.mask.shape,-1,dtype=int);self.mapping[self.mask]=np.arange(self.n)
        basis=[field.basis(s) for s in stations]
        self.centers=np.stack([p[0] for p in basis])/field.R
        self.basis=np.stack([np.stack((p[1],p[2]),axis=1) for p in basis])
        self.refs=self.references[self.mask]/field.R
        self.initial=self.original_xy[self.mask]/field.R
        self.fixed=self.original_world/field.R
        self.li=self.mapping[:-1].ravel();self.ri=self.mapping[1:].ravel()
        self.affected=(self.li>=0)|(self.ri>=0)
        self.owner=np.tile(np.arange(187),len(stations)-1)
        disk_world=self.world(self.refs.ravel())
        a,b=disk_world[:-1].reshape(-1,3),disk_world[1:].reshape(-1,3)
        self.move=np.where(self.affected,DISK_M/field.R,0.)
        self.i,self.j=pair_candidates(a,b,self.affected,self.owner,self.move,
                                     (2*field.r+CONSTRUCTION_TARGET_M)/field.R)
        self.wa=wrapper_leaves['a']/field.R;self.wb=wrapper_leaves['b']/field.R
        self.we=wrapper_leaves['error']/field.R
        self.wi,self.wj=wrapper_candidates(a,b,self.affected,self.move,self.wa,self.wb,self.we,
                                          (field.r+16e-6+CONSTRUCTION_TARGET_M)/field.R)
        # Bernstein hull controls give continuous straight-body containment
        # over every guide piece, including the window's unchanged left seam.
        es=[];fractions=[];controls=[]
        for station in range(self.anchor,len(stations)-1):
            segments=station*187+np.flatnonzero(self.affected[station*187:(station+1)*187])
            for t,p in guide_pieces(field,stations[station],stations[station+1]):
                es.append(np.repeat(segments,len(t)));fractions.append(np.tile(t,len(segments)))
                controls.append(np.tile(p,(len(segments),1))/field.R)
        self.es=np.concatenate(es);self.et=np.concatenate(fractions);self.ec=np.concatenate(controls)

    def world(self,flat):
        x=np.asarray(flat).reshape(self.n,2)
        result=self.fixed.copy()
        result[self.mask]=self.centers[self.node_station]+np.einsum('nij,nj->ni',self.basis[self.node_station],x)
        return result

    def segments(self,flat):
        world=self.world(flat)
        return world[:-1].reshape(-1,3),world[1:].reshape(-1,3)

    def evaluate(self,flat,jacobian=True,actual_world=None):
        x=np.asarray(flat).reshape(self.n,2);R=self.field.R
        if actual_world is None:a,b=self.segments(flat)
        else:
            world=actual_world/R;a,b=world[:-1].reshape(-1,3),world[1:].reshape(-1,3)
            x=np.einsum('ni,nij->nj',world[self.mask]-self.centers[self.node_station],self.basis[self.node_station])
        values=[];rows=[];cols=[];entries=[];offset=0
        def add(raw,terms):
            nonlocal offset
            values.append(np.maximum(raw,0))
            if jacobian:
                active=raw>0
                for ids,gradient in terms:
                    live=(ids>=0)&active
                    if not live.any():continue
                    ids0=ids[live]
                    projected=np.einsum('ni,nij->nj',gradient[live],self.basis[self.node_station[ids0]])
                    for k in range(2):
                        rows.append(offset+np.flatnonzero(live));cols.append(2*ids0+k);entries.append(projected[:,k])
            offset+=len(raw)
        delta=x-self.refs;length=np.linalg.norm(delta,axis=1)
        # A1nm solver inset covers world quantization while the authoritative
        # acceptance still uses the original30um disks and original references.
        move_raw=length-(DISK_M-ROUNDING_INSET_M)/R
        values.append(np.maximum(move_raw,0));offset+=self.n
        if jacobian:
            gradient=delta/np.maximum(length[:,None],1e-30)
            active=move_raw>0
            for k in range(2):
                rows.append(np.flatnonzero(active));cols.append(2*np.flatnonzero(active)+k);entries.append(gradient[active,k])
        self.last_body_lower=np.inf
        for first in range(0,len(self.i),BATCH):
            i=self.i[first:first+BATCH];j=self.j[first:first+BATCH]
            d,s,t,n,e,_=segment_closest(a[i],b[i],a[j],b[j])
            self.last_body_lower=min(self.last_body_lower,float((d-e-2*self.field.r/R).min()))
            raw=(2*self.field.r+CONSTRUCTION_TARGET_M)/R-(d-e)
            add(raw,[(self.li[i],-(1-s)[:,None]*n),(self.ri[i],-s[:,None]*n),
                     (self.li[j],(1-t)[:,None]*n),(self.ri[j],t[:,None]*n)])
        self.last_wrap_lower=np.inf
        for first in range(0,len(self.wi),BATCH):
            i=self.wi[first:first+BATCH];j=self.wj[first:first+BATCH]
            d,s,t,n,e,_=segment_closest(a[i],b[i],self.wa[j],self.wb[j])
            lower=d-e-self.we[j]-(self.field.r+16e-6)/R
            self.last_wrap_lower=min(self.last_wrap_lower,float(lower.min()))
            add(CONSTRUCTION_TARGET_M/R-lower,[(self.li[i],-(1-s)[:,None]*n),(self.ri[i],-s[:,None]*n)])
        diff=(1-self.et)[:,None]*a[self.es]+self.et[:,None]*b[self.es]-self.ec
        length=np.linalg.norm(diff,axis=1);normal=diff/np.maximum(length[:,None],1e-30)
        self.last_envelope=float(length.max()*R+self.field.r)
        add(length-(R-self.field.r-ROUNDING_INSET_M)/R,
            [(self.li[self.es],(1-self.et)[:,None]*normal),(self.ri[self.es],self.et[:,None]*normal)])
        result=np.concatenate(values)
        if not jacobian:return result
        matrix=coo_matrix((np.concatenate(entries) if entries else [],
                          (np.concatenate(rows) if rows else [],np.concatenate(cols) if cols else [])),
                          shape=(len(result),2*self.n)).tocsr()
        matrix.eliminate_zeros()
        return result,matrix

    def assess_actual(self,flat):
        exact=self.world(flat)*self.field.R
        actual=exact.astype('f4')
        # Crucial: quantization of moving points must never rewrite any fixed
        # historical coordinate. Their original stored bytes are authoritative.
        actual[~self.mask]=self.original_world[~self.mask].astype('f4')
        residual=self.evaluate(flat,False,actual.astype('f8'))
        reference_world=self.centers[self.node_station]*self.field.R+np.einsum(
            'nij,nj->ni',self.basis[self.node_station],self.refs*self.field.R)
        movement=np.linalg.norm(actual[self.mask].astype(float)-reference_world,axis=1)
        fixed_equal=np.array_equal(actual[~self.mask],self.original_world[~self.mask].astype('f4'))
        body=self.last_body_lower*self.field.R;wrap=self.last_wrap_lower*self.field.R
        return actual,{'max_residual_m':float(residual.max()*self.field.R),
                       'affected_body_clearance_lower_m':min(body,CONSTRUCTION_TARGET_M),
                       'affected_wrapper_clearance_lower_m':min(wrap,CONSTRUCTION_TARGET_M),
                       'positive_physical_clearance_pass':body>0 and wrap>0,
                       'construction_reserve_pass':body>=RESERVE_M and wrap>=RESERVE_M,
                       'continuous_straight_surface_envelope_upper_m':self.last_envelope+1e-12,
                       'continuous_straight_containment_pass':self.last_envelope+1e-12<=self.field.R,
                       'max_movement_from_original_reference_m':float(movement.max()),
                       'original30um_disks_pass':bool(movement.max()<=DISK_M+1e-12),
                       'fixed_world_bytes_preserved':fixed_equal,'native_curve_qualified':False}
