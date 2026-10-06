"""Sequential bounded-memory procedural maps. No image inputs or lighting in color."""
import argparse,gc,hashlib,json,sys
from pathlib import Path
import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from surface_math import field,clamp,smooth,voronoi,blend,color,TAU,specks,paths
from expansion.catalog import SPECS,get
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'expansion'

def generate(s,n):
    x=(np.arange(n,dtype=np.float32)[None,:]+.5)/n
    y=(np.arange(n,dtype=np.float32)[:,None]+.5)/n
    seed=s['seed'];family=s['family'];fine=field(n,64,seed+1)
    broad=field(n,5,seed+2,2)
    h=.5+.01*fine;rough=s['roughness']+.008*fine
    rgb=color(s['color'],.005*broad);metal=np.full((n,n),s['metallic'],np.float32)
    ao=np.ones((n,n),np.float32);opacity=ao.copy();extras={}
    if family=='Bark':
        # Nonuniform staggered, vertically elongated plates plus split fiber bands.
        u=x+.023*field(n,(12,5),seed+3,2)
        v=y+.032*field(n,(8,9),seed+4,2)
        gap,dist,tone=voronoi(n,12,seed+5,u,v*.35)
        fissure=1-smooth(gap,.025,.12)
        axial=np.sin(TAU*(u*105+.42*field(n,(25,8),seed+6)))
        splinter=smooth(axial,.15,.94)*(1-fissure)
        fiber_field=field(n,(34,95),seed+7,2)
        fine_splits=smooth(-fiber_field,.18,.58)*(1-fissure)
        finite_fibers=paths(n,seed+18,1800,width=(.00022,.00048),length=(.007,.12),angle=np.pi/2,jitter=.12,bend=.045)
        h=.65-.43*fissure+.055*tone+.008*splinter+.055*fiber_field-.025*fine_splits-.024*finite_fibers
        rgb=color(s['color'],.04*tone+.025*broad+.018*fiber_field)
        # intrinsic fresh plate lip color, no painted cavity shadow.
        rgb=blend(rgb,(.37,.235,.115),splinter*.045+finite_fibers*.07)
        rough=.85+.05*fissure+.02*fine;ao=1-.42*fissure
    elif family in ('Concrete','Cement paste','Asphalt') or s['id']=='26_tile_terrazzo':
        gap,d,tone=voronoi(n,36 if family=='Asphalt' else 25,seed+8)
        inclusion=smooth(d,.08,.17)*(1-smooth(d,.20,.42))
        # Stable mineral grains: angular shapes, grain-sized correlated color.
        aggregate=smooth(gap,.015,.06)*smooth(tone,-.12,.04)
        fine_grain=field(n,170,seed+9)
        paste=color(s['color'],.012*broad+.006*fine_grain)
        mineral=color((.33,.325,.29) if family=='Asphalt' else (.56,.54,.47),.095*tone)
        visible=aggregate
        if family=='Concrete' and s.get('finish') not in ('polished','exposed_aggregate'):visible=aggregate*.12
        rgb=paste*(1-visible[...,None])+mineral*visible[...,None]
        h=.5+.035*aggregate+.012*fine_grain;rough=s['roughness']+.05*aggregate+.01*broad
        finish=s.get('finish','')
        if family=='Asphalt':h=.37+.36*aggregate+.04*tone*aggregate;rough=.86+.07*(1-aggregate)
        elif finish=='board_formed':
            board=(y*8)%1;seam=1-smooth(np.minimum(board,1-board),.007,.025)
            grain=np.sin(TAU*(x*60+.55*field(n,(4,30),seed+10)))
            h+=.04*grain+.17*seam;rgb+=.006*grain[...,None]
        elif finish=='brushed':h+=.10*np.sin(TAU*(y*180+.16*field(n,(9,18),seed+10)))
        elif finish=='polished':rough=.25+.05*aggregate;h=.5+.006*fine_grain
        elif finish=='exposed_aggregate':h=.32+.43*aggregate+.035*tone*aggregate
        elif finish=='weathered':
            erosion=smooth(field(n,12,seed+10,2),.0,.5)
            h-=.16*erosion*(1-aggregate);rough+=.06*erosion
            rgb=blend(rgb,(.40,.42,.37),erosion*.15)
        elif finish=='paste':rgb=paste;h=.5+.085*field(n,(9,6),seed+11,2);rough=.7+.035*broad
        else:
            pores=specks(n,seed+11,90,(.001,.003),(.001,.0025))
            h-=.18*pores;ao-=.18*pores
        extras['AggregateMask']=aggregate if finish!='paste' else np.zeros((n,n),np.float32)
    elif family=='Plastic':
        if s['variant']=='Molded stipple':
            gap,d,tone=voronoi(n,140,seed+12)
            h=.50+.21*(1-smooth(d,.04,.54));rough=.49+.045*smooth(d,.1,.5)
        else:h=.5+.025*field(n,18,seed+12);rough=.28+.008*broad
    elif family=='PETG':
        phase=(y*s['tile_m']/s['layer_height_m'])%1
        bead=np.sqrt(clamp(1-((phase-.5)/.52)**2))
        h=.22+.55*bead+.008*fine
        boundary=np.broadcast_to(1-smooth(np.minimum(phase,1-phase),.015,.12),(n,n)).copy()
        rough=s['roughness']+.14*boundary+.006*fine
        extras.update(Transmission=np.full((n,n),s.get('transmission',0),np.float32),
                      Thickness=np.ones((n,n),np.float32),PrintBoundary=boundary)
    elif family in ('LCD','CRT'):
        pitch=s['pixel_pitch_m'];px=x*s['tile_m']/pitch;py=y*s['tile_m']/pitch
        sub=(px*3)%1;row=py%1
        mask=(1-smooth(np.abs(sub-.5),.31,.44))*(1-smooth(np.abs(row-.5),.30,.45))
        if family=='CRT':mask=(1-smooth(np.abs(sub-.5),.29,.43))*(1-smooth(np.abs(row-.5),.42,.49))
        stripe=np.broadcast_to((np.floor(px*3).astype(int)%3),(n,n))
        # Analytic color bars and circle, not a projected reference image.
        test=np.stack(np.broadcast_arrays(.3+.6*x,.3+.6*y,.35+.45*np.sin(TAU*x)**2),axis=-1)
        ring=1-smooth(np.abs(np.sqrt((x-.5)**2+(y-.5)**2)-.25),.009,.019)
        test=test*(1-ring[...,None])+ring[...,None]*.95
        emission=np.zeros((n,n,3),np.float32)
        for c in range(3):emission[...,c]=(stripe==c)*mask*test[...,c]
        extras['Emission']=emission;extras['PhosphorMask' if family=='CRT' else 'SubpixelMask']=mask
        extras['Transmission']=np.full((n,n),.96,np.float32)
        extras['Thickness']=np.ones((n,n),np.float32)
        h=.5+.08*mask;rough=.28+.1*(1-mask)
    elif family=='Tile':h=.5+.012*field(n,30,seed+13);rough=s['roughness']+.012*broad
    elif family=='Engineered':h=.5+.012*field(n,(6,160),seed+13);rough=s['roughness']+.012*broad
    elif family in ('Blanket','Carpet'):h=.5+.015*fine;rough=.85+.012*broad
    elif family=='Grass':rgb=color(s['color'],.025*broad);rough=.72+.02*broad
    # Geometric openings are deliberately not encoded as painted/alpha holes.
    return dict(BaseColor=clamp(rgb).astype(np.float32),Height=clamp(h).astype(np.float32),
                Roughness=clamp(rough).astype(np.float32),Metallic=metal,AO=clamp(ao),Opacity=opacity,**extras)

def normals(h,s):
    n=h.shape[0];dx=(np.roll(h,-1,1)-np.roll(h,1,1))*n*.5*s['height_scale_m']/s['tile_m']
    dy=(np.roll(h,-1,0)-np.roll(h,1,0))*n*.5*s['height_scale_m']/s['tile_m']
    v=np.stack((-dx,dy,np.ones_like(h)),axis=-1);v/=np.linalg.norm(v,axis=-1,keepdims=True)
    return (.5+.5*v).astype(np.float32)

def write(s,n,destination=None):
    out=(destination or OUT/'maps')/s['id'];out.mkdir(parents=True,exist_ok=True)
    values=generate(s,n);checks={}
    macro=gaussian_filter(values['Height'],max(.5,n/512),mode='wrap')
    values.update(Height_Macro=macro,Normal_OpenGL=normals(values['Height'],s),
                  Normal_Micro_OpenGL=normals(values['Height']-macro,s))
    direct=values['Normal_OpenGL'].copy();direct[...,1]=1-direct[...,1];values['Normal_DirectX']=direct
    values['ORM']=np.stack((values['AO'],values['Roughness'],values['Metallic']),axis=-1)
    map_info={}
    for name,a in values.items():
        assert a.shape[:2]==(n,n) and np.isfinite(a).all() and a.min()>=0 and a.max()<=1,name
        sixteen=name.startswith('Height') or name=='Thickness'
        encoded=np.rint(a*(65535 if sixteen else 255)).astype(np.uint16 if sixteen else np.uint8)
        path=out/(name+'.png');Image.fromarray(encoded).save(path,compress_level=4)
        checks[name]=dict(min=float(a.min()),max=float(a.max()),sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        map_info[name+'.png']=dict(color_space='sRGB' if name=='BaseColor' else 'linear / Non-Color',bit_depth=16 if sixteen else 8)
    metadata=dict(s,resolution=n,maps=map_info,height_midlevel=.5,
                  displacement='(height - 0.5) * height_scale_m',
                  normal_pairing='Height_Macro + Normal_Micro_OpenGL, or full Normal_OpenGL without displacement',
                  emission_encoding='linear RGB8 intensity pattern; strength in scene is separate',
                  thickness_scale_m=s.get('wall_thickness_m',s.get('glass_thickness_m',0)),
                  map_limits='Surface maps do not replace fibers, blades, links, lattice voids, panel seams or tile grout geometry.',
                  checks=checks)
    (out/'material.json').write_text(json.dumps(metadata,indent=2)+'\n')
    del values;gc.collect();return metadata

def main():
    p=argparse.ArgumentParser();p.add_argument('--resolution',type=int,default=512);p.add_argument('--only',nargs='*');p.add_argument('--parameters',type=Path);a=p.parse_args()
    for s in SPECS:
        if a.only and s['id'] not in a.only:continue
        s=get(s['id'],a.parameters);m=write(s,a.resolution);print(s['id'],a.resolution,'validated',flush=True)
    (OUT/'catalog.json').write_text(json.dumps(SPECS,indent=2)+'\n')
if __name__=='__main__':main()
