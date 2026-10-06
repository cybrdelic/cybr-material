"""Sample a timber growth volume and project finite anatomical openings.

The trunk axis is Y. A board samples the XZ radial cross-section at a varying
cut depth. Annual growth, earlywood vessels, ray plates and finish wear share
that volume. No color image is an input, and color is not used to infer height.
"""
import json,math
from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter, map_coordinates
from surface_math import clamp,smooth,field,color,blend,paths

CONFIG=Path(__file__).resolve().parent/'wood_profiles.json'


def growth_volume(x,y,n,p):
    rng=np.random.default_rng(p['seed'])
    width=p['extent_m'];length=p['extent_m']
    axis_y=(y-.5)*length
    # Small trunk curvature and a sloping saw plane change radial depth.
    # Radius, rather than a shared wave phase, determines every growth band.
    center_x=p['pith_x_m']+.0009*np.sin(axis_y*8.5)+.0004*np.sin(axis_y*23)
    center_z=.0016*np.sin(axis_y*7.1+.8)
    cut_z=p['cut_depth_m']+p['cut_slope']*axis_y+.004*(axis_y/length)**2
    dx=x*width-center_x;dz=cut_z-center_z
    radius=np.sqrt((dx*.985)**2+(dz*1.04)**2)
    theta=np.arctan2(dx,dz)
    # Growth widths are correlated across neighboring years, but the
    # spacing and latewood fraction are not constant stripes.
    widths=np.exp(rng.normal(np.log(p['ring_width_m']),.32,240))
    widths=np.clip(gaussian_filter(widths,.70,mode='nearest'),.0010,.009)
    levels=np.concatenate(([0],np.cumsum(widths)))
    distorted=radius+.00065*np.sin(theta*3+radius*63)+.00035*np.sin(theta*7-radius*29)
    for branch in p.get('branches',[]):
        # Distance to an oblique cylindrical branch axis in XYZ. Its
        # intersection with the saw plane is elliptical. The transition
        # turns growth and fiber coordinates around a sound branch knot;
        # it adds no dark vignette or corresponding height bulge.
        bx=x*width-branch['x_m'];by=axis_y-branch['y_m']
        tilt=branch['axis_tilt'];cosine=np.sqrt(1-tilt*tilt)
        bz=cut_z-(p['cut_depth_m']+p['cut_slope']*branch['y_m'])
        branch_r=np.sqrt(bx*bx+(bz*cosine-by*tilt)**2)
        support=branch['radius_m']
        influence=1-smooth(branch_r,support*.64,support*3.8)
        at_center=np.sqrt(((branch['x_m']-p['pith_x_m'])*.985)**2+(p['cut_depth_m']*1.04)**2)
        branch_growth=at_center-support*2.2+branch_r*2.2
        distorted=distorted*(1-influence)+branch_growth*influence
    index=np.clip(np.searchsorted(levels,distorted,side='right')-1,0,len(widths)-1)
    phase=(distorted-levels[index])/widths[index]
    late_start=rng.uniform(.54,.77,len(widths)).astype(np.float32)
    late=smooth(phase,late_start[index],.98)
    early=np.exp(-((phase-.17)/(.14 if p['ring_porous'] else .27))**2)
    annual=gaussian_filter(rng.uniform(-1,1,len(widths)),1.05,mode='nearest').astype(np.float32)[index]
    density=late*(.76+.24*rng.uniform(0,1,len(widths)).astype(np.float32)[index])
    return density.astype(np.float32),early.astype(np.float32),annual.astype(np.float32),theta.astype(np.float32),distorted.astype(np.float32)


def fiber_bundles(radius,y,n,p,seed,width_mm,length_mm):
    """Correlated cell bundles in radial/longitudinal coordinates.

    Sampling radial distance bends the bundles with the actual board cut.
    A bundle continues over many vessel openings; it is not a row of dots.
    The grid spacing is set in meters and does not depend on output resolution.
    """
    sx=max(8,math.ceil(float(np.max(radius))/(width_mm*.001))+2)
    sy=max(4,round(p['extent_m']/(length_mm*.001)))
    rng=np.random.default_rng(seed)
    grid=rng.uniform(-1,1,(sy+4,sx+4)).astype(np.float32)
    u=radius/(width_mm*.001)
    v=np.broadcast_to(y,(n,n))*p['extent_m']/(length_mm*.001)
    return map_coordinates(grid,np.stack([v+1,u+1]),order=3,mode='reflect',prefilter=True).astype(np.float32)


def vessels(n,p,probability,seed,count,radius_um,length_mm):
    """Longitudinal lumen windows with rounded floors and varied intersections.

    Diameter and depth are tied to anatomical radius. Features below a pixel
    retain fractional coverage instead of becoming full-strength pixel dots.
    """
    rng=np.random.default_rng(seed);tile=p['extent_m']
    mask=np.zeros((n,n),np.float32);depth=mask.copy()
    for _ in range(count):
        cx,cy=rng.uniform(0,n,2)
        if rng.random()>probability[int(cy),int(cx)]:continue
        r=np.exp(rng.uniform(np.log(radius_um[0]),np.log(radius_um[1])))*1e-6
        opening=np.exp(rng.uniform(np.log(length_mm[0]),np.log(length_mm[1])))*.001
        rx=r*n/tile;ry=opening*.5*n/tile
        tilt=rng.uniform(-.016,.016)
        drift=tilt*ry
        ix=np.arange(max(0,int(cx-max(rx,.65)-abs(drift)-2)),min(n,int(cx+max(rx,.65)+abs(drift)+3)))
        iy=np.arange(max(0,int(cy-ry-2)),min(n,int(cy+ry+3)))
        along=(iy[:,None]-cy)/max(ry,.65)
        bend=drift*along+rng.uniform(-.30,.30)*rx*along**2
        across=(ix[None,:]-cx-bend)/max(rx,.65)
        # Elliptical intersections taper continuously, rather than having
        # parallel sides and abruptly rounded rectangular ends.
        radius_at_y=np.sqrt(clamp(1-along**2))
        coverage=1-smooth(abs(across),radius_at_y*.54,radius_at_y+.26)
        coverage*=1-smooth(abs(along),.92,1.04)
        aa=min(1,rx/.65)*min(1,ry/.65)
        coverage*=aa
        floor=np.sqrt(clamp(1-across**2-along**2))
        # Some lumen windows are partly sealed by finish or deposits.
        filled=.28 if rng.random()<.18 else 1
        relief=min(r*.40,62e-6)*filled*rng.uniform(.50,.95)
        target=np.ix_(iy,ix)
        mask[target]=np.maximum(mask[target],coverage)
        depth[target]=np.maximum(depth[target],floor*relief*aa)
    return mask,depth


def rays(n,p,seed,count,length_mm,width_mm):
    """Finite medullary plates: broad bodies with tapered broken ends."""
    rng=np.random.default_rng(seed);tile=p['extent_m']
    out=np.zeros((n,n),np.float32)
    for _ in range(count):
        cx,cy=rng.uniform(0,n,2)
        rx=np.exp(rng.uniform(np.log(length_mm[0]),np.log(length_mm[1])))*.001*n/tile*.5
        ry=np.exp(rng.uniform(np.log(width_mm[0]),np.log(width_mm[1])))*.001*n/tile*.5
        ix=np.arange(max(0,int(cx-rx-2)),min(n,int(cx+rx+3)))
        iy=np.arange(max(0,int(cy-max(ry,.65)-rx*.035-3)),min(n,int(cy+max(ry,.65)+rx*.035+4)))
        u=(ix[None,:]-cx)/max(rx,.65)
        v=(iy[:,None]-cy-rng.uniform(-.028,.028)*(ix[None,:]-cx))/max(ry,.65)
        edge=np.maximum(.025,(.73+.11*np.sin(u*7+rng.uniform(0,6.28))+.07*np.cos(u*13))*clamp(1-u*u)**.52)
        body=(1-smooth(abs(v),edge*.65,edge))*(1-smooth(abs(u),.73,1.02))
        body*=min(1,ry/.65)*rng.uniform(.4,1)
        target=np.ix_(iy,ix);out[target]=np.maximum(out[target],body)
    return out


def wood(x,y,n,oak=False):
    name='fumed_oak' if oak else 'american_walnut'
    p=json.loads(CONFIG.read_text())[name];seed=p['seed'];tile=p['extent_m']
    density,early,annual,theta,radius=growth_volume(x,y,n,p)
    bundles=fiber_bundles(radius,y,n,p,seed+14,1.3,38)
    fibers=fiber_bundles(radius,y,n,p,seed+15,.34,14)
    cellular=fiber_bundles(radius,y,n,p,seed+19,.13,5.5)
    # The bulk color and pore distribution come from the same growth volume.
    probability=clamp(.018+.982*early**1.6) if oak else clamp(.22+.70*early)
    probability*=.36+.64*smooth(bundles,-.64,.58)
    pore,pore_m=vessels(n,p,probability,seed+11,p['vessel_count'],p['vessel_radius_um'],p['vessel_length_mm'])
    if oak:
        small,small_m=vessels(n,p,np.ones((n,n),np.float32)*.32,seed+12,38000,(25,60),(.35,2.1))
        pore=np.maximum(pore,small*.32);pore_m=np.maximum(pore_m,small_m)
    ray=rays(n,p,seed+13,p['ray_count'],p['ray_length_mm'],p['ray_width_mm'])
    # Low-amplitude coherent cell-bundle variation has a fixed finite scale.
    # It never becomes independent per-pixel color or white-noise bump.
    contact=smooth(field(n,(6,8),seed+16,2),.06,.52)
    scuff=paths(n,seed+17,380,width=(.00009,.00025),length=(.003,.025),angle=.9,jitter=.14,bend=.004)*contact
    checks=paths(n,seed+18,45,width=(.00016,.00044),length=(.004,.032),angle=math.pi/2,jitter=.025,bend=.006)*(.22+.78*early)
    rgb=color(tuple(p['base_srgb']),.033*annual+.021*bundles+.016*fibers+.006*cellular)
    late_tone=clamp(density*(.25 if oak else .22)*smooth(bundles+.5*fibers,-.7,.7))
    rgb=blend(rgb,tuple(p['late_srgb']),late_tone)
    rgb=blend(rgb,tuple(p['pore_srgb']),pore*(.48 if oak else .24)+checks*.08)
    rgb=blend(rgb,tuple(p['ray_srgb']),ray*(.24 if oak else .07))
    rgb=blend(rgb,tuple(p['abraded_srgb']),scuff*.04)
    # Relief is independently authored in meters from lumen and finish models.
    # Broad pigment bands do not become carved trenches.
    meters=2e-6*fibers+1e-6*cellular+3e-6*density+3e-6*ray
    meters-=pore_m*(1-contact*.10)+checks*65e-6+scuff*5e-6
    h=.56+meters/p['height_scale_m']
    rough=p['roughness']+.075*pore+.09*checks-.057*contact+.105*scuff+.012*ray+.012*cellular+.018*bundles
    ao=clamp(1-.12*pore-.19*checks)
    wear=clamp(contact*.38+scuff+.42*checks)
    return rgb,h,rough,ao,np.zeros((n,n),np.float32),wear
