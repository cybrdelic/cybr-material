"""Post-solve capsule separation outside declared fused junction disks.

This detects inadmissible geometry. It supplies no contact force and cannot
continue a contacting branch. Incident centerlines are trimmed at radius 2t.
Their finite-thickness capsules can extend to radius 2.5t; that larger disk
is the unresolved fused junction region. No absence-of-overlap claim is made
inside it.
"""
import numpy as np


def cross(a,b): return a[...,0]*b[...,1]-a[...,1]*b[...,0]


def segment_distances(a,b,c,d):
    u=b-a;v=d-c;r=c-a;den=cross(u,v)
    nonparallel=abs(den)>1e-25
    safe=np.where(nonparallel,den,1.)
    ta=cross(r,v)/safe;tb=cross(r,u)/safe
    intersects=nonparallel&(ta>=0)&(ta<=1)&(tb>=0)&(tb<=1)
    def pd(p,p0,p1):
        edge=p1-p0
        t=np.clip(np.sum((p-p0)*edge,axis=-1)/np.sum(edge*edge,axis=-1),0,1)
        return np.linalg.norm(p-p0-t[...,None]*edge,axis=-1)
    dist=np.minimum.reduce([pd(a,c,d),pd(b,c,d),pd(c,a,b),pd(d,a,b)])
    return np.where(intersects,0.,dist)


def outside_disk(points,center,radius):
    pieces=[]
    for idx,(a,b) in enumerate(zip(points[:-1],points[1:])):
        u=b-a;v=a-center;aa=u@u;bb=2*(u@v);cc=v@v-radius**2
        disc=bb*bb-4*aa*cc;cuts=[0.,1.]
        if disc>0:
            cuts.extend(t for t in ((-bb-np.sqrt(disc))/(2*aa),(-bb+np.sqrt(disc))/(2*aa)) if 0<t<1)
        cuts=sorted(cuts)
        for lo,hi in zip(cuts[:-1],cuts[1:]):
            if np.linalg.norm(a+(lo+hi)/2*u-center)>=radius:
                pieces.append((a+lo*u,a+hi*u,idx))
    return pieces


def check_contact(network,q):
    n=network;p=n.p;t=p.wall_thickness_m
    x=q[:2*n.nn].reshape(-1,2)*n.L
    lines=[x[w['nodes']] for w in n.walls]
    minimum=np.inf;limiting=None;max_curvature=0.;pairs=0
    for wi,line in enumerate(lines):
        ds=np.diff(line,axis=0);ell=np.linalg.norm(ds,axis=1)
        theta=np.unwrap(np.arctan2(ds[:,1],ds[:,0]))
        max_curvature=max(max_curvature,float(np.max(abs(np.diff(theta))/((ell[:-1]+ell[1:])/2))*t/2))
        ref=n.X[n.walls[wi]['nodes']]*n.L
        arc=np.r_[0.,np.cumsum(np.linalg.norm(np.diff(ref,axis=0),axis=1))]
        i,j=np.triu_indices(len(ds),2)
        keep=arc[j]-arc[i+1]>=2*t;i=i[keep];j=j[keep]
        if len(i):
            dd=segment_distances(line[i],line[i+1],line[j],line[j+1]);k=np.argmin(dd);pairs+=len(dd)
            if dd[k]<minimum:minimum=float(dd[k]);limiting=[wi,int(i[k]),wi,int(j[k])]
        for wj in range(wi+1,len(lines)):
            other=lines[wj];common=set(n.walls[wi]['junctions'])&set(n.walls[wj]['junctions'])
            if common:
                vertex=next(iter(common));radius=p.junction_contact_radius_m
                aa=outside_disk(line,x[vertex],radius);bb=outside_disk(other,x[vertex],radius)
                if not aa or not bb:continue
                a,b,ai=zip(*aa);c,d,bi=zip(*bb)
                a=np.array(a);b=np.array(b);c=np.array(c);d=np.array(d)
            else:
                a,b=line[:-1],line[1:];c,d=other[:-1],other[1:]
                ai=np.arange(len(a));bi=np.arange(len(c))
            dd=segment_distances(a[:,None,:],b[:,None,:],c[None,:,:],d[None,:,:])
            k=np.unravel_index(np.argmin(dd),dd.shape);pairs+=dd.size
            if dd[k]<minimum:minimum=float(dd[k]);limiting=[wi,int(ai[k[0]]),wj,int(bi[k[1]])]
    areas=[]
    for polygon in n.cell_polygons(q):
        poly=polygon*n.L
        areas.append(float(np.sum(cross(poly,np.roll(poly,-1,axis=0)))/2))
    result=dict(minimum_nonlocal_centerline_distance_m=float(minimum),
                minimum_surface_gap_m=float(minimum-t),wall_thickness_m=t,
                limiting_segment_pair=limiting,checked_segment_pairs=int(pairs),
                max_half_thickness_current_curvature=float(max_curvature),
                cell_enclosed_midline_areas_m2=areas,
                incident_centerline_trim_radius_m=p.junction_contact_radius_m,
                unresolved_fused_junction_radius_m=p.junction_contact_radius_m+t/2,
                same_wall_material_interval_gap_exclusion_m=2*t,
                contact_forces_implemented=False,
                passed=bool(minimum>t and max_curvature<1 and min(areas)>0))
    return result
