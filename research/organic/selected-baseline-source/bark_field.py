"""Authored mature fir bark: anisotropic irregular fracture islands, not rows.
No external texture pixels. Distances, dimensions and colors are authored.
"""
from pathlib import Path
import numpy as np
from scipy.ndimage import zoom,gaussian_filter
from scipy.spatial import cKDTree
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]

def noise(rng,n,shape):
 a=rng.random(shape);a=zoom(a,(n/shape[0],n/shape[1]),order=3)[:n,:n]
 return (a-a.mean())/(a.std()+1e-9)

def build(n=641,seed=604):
 rng=np.random.default_rng(seed);y,x=np.mgrid[0:1:complex(n),0:1:complex(n)]
 # Irregular, locally leaning fibrous axis and variable Voronoi cell sizes.
 wx=x+.017*noise(rng,n,(11,14))+.006*noise(rng,n,(35,28))
 wy=y+.016*noise(rng,n,(12,17))
 # Mature fir-inspired elongated ridge islands. No uniform vertical tracks.
 points=rng.uniform([-.12,-.12],[1.12,1.12],(125,2));pts=points.copy();pts[:,1]/=2.8
 q=np.stack((wx,wy/2.8),-1).reshape(-1,2)
 ds,ids=cKDTree(pts).query(q,k=2);ds=ds.reshape(n,n,2);ids=ids.reshape(n,n,2)
 edge=(ds[:,:,1]-ds[:,:,0])*0.5*.18
 # Vary depth and shoulder width independently: ridges do not share a cross-section.
 widths=rng.uniform(.0011,.0048,len(pts))[ids[:,:,0]]
 body=np.power(1-np.exp(-edge/widths),.68)
 amplitude=rng.uniform(.0046,.0115,len(pts))[ids[:,:,0]]
 broad=noise(rng,n,(19,29));chip=noise(rng,n,(81,57));grain=noise(rng,n,(240,180))
 h=.0023+amplitude*body+.0007*broad*body+.00038*chip+.00017*grain
 # Local subordinate fissures fade toward their ends and never become a row grid.
 subpts=rng.uniform([-.1,-.1],[1.1,1.1],(560,2));subpts[:,1]/=2.0
 dd,ii=cKDTree(subpts).query(np.stack((wx,wy/2.),-1).reshape(-1,2),k=2)
 dd=dd.reshape(n,n,2);ii=ii.reshape(n,n,2)
 sd=(dd[:,:,1]-dd[:,:,0])*.18*.5
 active=np.clip(noise(rng,n,(26,17))+.38,0,1)
 checks=np.exp(-(sd/.00045)**1.1)*active*body
 h-=checks*.0015
 # Small brittle chips: high-frequency angular ridges, recess-dependent amplitude.
 h+=.00042*np.maximum(chip-.25,0)*body
 h=np.maximum(h,.0013)
 shade=np.clip(body,.0,1);cellcol=rng.uniform(-.022,.022,len(pts))[ids[:,:,0]]
 base=np.array([.102,.073,.047]);top=np.array([.248,.194,.132])
 rgb=base[None,None,:]+shade[:,:,None]*(top-base)[None,None,:]
 rgb+=(.011*broad+.008*chip+cellcol)[:,:,None]*body[:,:,None]
 # Weathered gray crowns mixed with warmer newly fractured cork.
 gray=np.clip(.43+.25*noise(rng,n,(22,16)),0,.82)*body
 rgb=rgb*(1-gray[:,:,None])+np.array([.22,.21,.185])*gray[:,:,None]
 rgb=np.clip(rgb-checks[:,:,None]*.018,.014,1)
 np.savez_compressed(ROOT/'data/bark_r4.npz',height=h.astype('f4'),color=rgb.astype('f4'))
 stats={'seed':seed,'grid':n,'specimen_width_m':.18,'primary_cells':len(pts),'subordinate_cells':len(subpts),'height_range_m':[float(h.min()),float(h.max())],'model':'anisotropic branching junction network with variable ridge widths, amplitudes and terminated secondary fissures','photographic_inputs':False}
 (ROOT/'receipts/bark_r4_field.json').write_text(json.dumps(stats,indent=2))
 print(stats)
if __name__=='__main__':build()
