"""Metric graded mineral fragments seated below geometry-derived binder contact.

No photograph, scanned surface or decorative random field is used. Three seeded
fragment bands occupy different coarse-grain channel widths. A shared highest
surface winner supplies exposed area, mineral identity, roughness, and relief.
All coefficients are authored, not measured rheology or solved particle packing.
"""
import numpy as np
from scipy.ndimage import gaussian_filter, distance_transform_edt
from fracture_field import hash01

PALETTE=np.array([[.245,.253,.25],[.37,.36,.325],[.165,.18,.178],[.47,.458,.403],[.285,.30,.29],[.52,.515,.47]],'f4')
BANDS=(
    dict(name='contact_fines',pitch=.00034,radius=(.000085,.000170),peak=(.000042,.000100),probability=.92),
    dict(name='channel_grit',pitch=.00068,radius=(.000140,.000290),peak=(.000070,.000160),probability=.78),
    dict(name='open_pocket_chips',pitch=.00119,radius=(.000240,.000465),peak=(.000110,.000230),probability=.54),
)

def smooth(x,a,b):
    t=np.clip((x-a)/(b-a),0,1);return t*t*(3-2*t)

def contact_state(mask,tile):
    """Coarse solid support and periodic matrix-channel width, in SI units."""
    n=len(mask);pixel=tile/n
    support=gaussian_filter(mask,.0009/pixel,mode='wrap').astype('f4')
    broad=gaussian_filter(mask,.0030/pixel,mode='wrap').astype('f4')
    # Support only needs the local channel width through 1.5 mm. Padding avoids
    # seam errors without a nine-tile full-resolution distance transform.
    pad=max(2,int(np.ceil(.0016/pixel)))
    free=np.pad(mask<.5,pad,mode='wrap')
    clearance=distance_transform_edt(free,sampling=pixel)[pad:-pad,pad:-pad].astype('f4')
    clearance=np.minimum(clearance,.0015)
    pressure=np.clip(.72*support+.28*broad,0,1)
    # The meniscus is a geometric binder-fill response to nearby solid contact.
    bed=(.000010+.000018*support+.000022*(1-smooth(clearance,.00006,.00034))).astype('f4')
    return pressure,clearance,bed

def graded_fragment_field(n,tile,seed,coarse_mask,strip_rows=128):
    if not isinstance(n,int) or n<4 or not np.isfinite(tile) or tile<=0 or not isinstance(seed,int) or seed<0:
        raise ValueError('Expected integer resolution>=4, finite positive tile and nonnegative seed')
    if coarse_mask.shape!=(n,n) or not np.isfinite(coarse_mask).all():
        raise ValueError('Finite square coarse mask must match resolution')
    if min(round(tile/b['pitch']) for b in BANDS)<3:
        raise ValueError('Tile must span at least three cells in every fine band')
    pressure,clearance,bed=contact_state(coarse_mask,tile)
    height=np.zeros((n,n),'f4');coverage=height.copy();exposure=height.copy()
    ids=np.zeros((n,n),'u1');band_ids=ids.copy();pixel=tile/n
    x=(np.arange(n,dtype='f4')+.5)*pixel;reports=[]
    for band_index,band in enumerate(BANDS,1):
        cells=round(tile/band['pitch']);pitch=tile/cells;ix=np.floor(x/pitch).astype('i4')
        band_seed=seed+band_index*1009
        for start in range(0,n,strip_rows):
            stop=min(start+strip_rows,n);y=x[start:stop];iy=np.floor(y/pitch).astype('i4')
            hh=height[start:stop];cc=coverage[start:stop];ee=exposure[start:stop];ii=ids[start:stop];bb=band_ids[start:stop]
            for oy in (-1,0,1):
                cy=iy[:,None]+oy
                for ox in (-1,0,1):
                    cx=ix[None,:]+ox;hx=cx%cells;hy=cy%cells
                    a=hash01(hx,hy,band_seed);b=hash01(hx,hy,band_seed+41);c=hash01(hx,hy,band_seed+79)
                    d=hash01(hx,hy,band_seed+131);q=hash01(hx,hy,band_seed+193)
                    centerx=(cx+.5+(a-.5)*.88)*pitch;centery=(cy+.5+(b-.5)*.88)*pitch
                    px=np.floor(centerx/pixel).astype('i4')%n;py=np.floor(centery/pixel).astype('i4')%n
                    load=pressure[py,px];gap=clearance[py,px]
                    lo,hi=band['radius'];radius=lo*(hi/lo)**d
                    if band_index==1:
                        probability=band['probability']*(.82+.18*load)
                    else:
                        probability=band['probability']*smooth(gap,radius*.30,radius*1.45)*(.92-.28*load)
                    present=q<probability
                    dx=x[None,:]-centerx;dy=y[:,None]-centery
                    angle=c*np.float32(2*np.pi);co=np.cos(angle);si=np.sin(angle)
                    u=dx*co+dy*si;v=(-dx*si+dy*co)/(.48+.43*b)
                    # Unequal opposing facets form finite elongated mineral shards,
                    # rather than uniformly sized symmetric hexagonal beads.
                    edge=np.minimum.reduce((radius*(.77+.22*a)-u,radius*(.74+.25*b)+u,
                        radius*(.71+.26*c)-(.40*u+.9165*v),radius*(.72+.25*d)+(.48*u+.8773*v),
                        radius*(.74+.20*b)-(-.53*u+.8480*v),radius*(.75+.20*c)+(-.58*u+.8146*v)))
                    cover=np.clip(edge/pixel+.5,0,1)*present
                    low,high=band['peak'];peak=low+(high-low)*c
                    # Compaction depth uses local support and grain aspect; the
                    # thin binder film leaves only part of each mineral exposed.
                    burial=peak*(.18+.28*load+.10*a)
                    plane=peak+u*(a-.5)*.075+v*(b-.5)*.075
                    crown=np.minimum(plane,np.maximum(edge,0)*.76)
                    above=np.maximum(crown-burial,0)
                    relief=above*cover
                    exposed=cover*smooth(above,.000006,.000047)
                    # AA fringe selects the same surface owner even if relief is
                    # zero; no default palette ID at a positive coverage edge.
                    winner=(relief+cover*1e-9)>(hh+cc*1e-9)
                    hh[:]=np.where(winner,relief,hh);cc[:]=np.where(winner,cover,cc)
                    ee[:]=np.where(winner,exposed,ee)
                    ii[:]=np.where(winner,np.minimum((hash01(hx,hy,band_seed+263)*6).astype('u1'),5),ii)
                    bb[:]=np.where(winner,band_index,bb)
        reports.append(dict(band,actual_pitch_m=pitch,cells_per_axis=cells,seed=band_seed))
    matrix=coarse_mask<.05;mass=float(matrix.sum())
    report={'bands':reports,'matrix_exposed_mean':float(exposure[matrix].mean()),
        'matrix_exposed_fraction_gt_half':float((exposure[matrix]>.5).mean()),
        'matrix_relief_mean_m':float(height[matrix].mean()),'maximum_relief_m':float(height.max()),
        'binder_bed_range_m':[float(bed.min()),float(bed.max())],
        'matrix_owner_fraction_by_band':{b['name']:float(np.count_nonzero((band_ids==k)&matrix)/mass) for k,b in enumerate(BANDS,1)},
        'coverage_pixels_without_owner':int(np.count_nonzero((coverage>0)&(band_ids==0))),
        'exposure_pixels_without_relief':int(np.count_nonzero((exposure>0)&(height<=0))),
        'state':'Three metric finite fragment bands with coarse-geometry channel grading; one surface owner supplies height, exposure and mineral identity',
        'optical_law':'Binder/mineral blend is controlled by exposed crown height above geometry-dependent burial, not independent color noise',
        'seed':seed,'packing_qualified':False}
    return dict(height_m=height,coverage=coverage,exposure=exposure,mineral_rgb=PALETTE[ids],
                binder_bed_m=bed,band_ids=band_ids,report=report)
