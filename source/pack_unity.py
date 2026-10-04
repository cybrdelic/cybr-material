"""Convert original ORM components to native Unity mask layouts.

python source/pack_unity.py --pipeline hdrp
python source/pack_unity.py --pipeline urp
Requires NumPy and Pillow. Outputs exports/unity_<pipeline>/<material>/.
"""
import argparse
from pathlib import Path
import numpy as np
from PIL import Image

ROOT=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser()
ap.add_argument('--pipeline',choices=['hdrp','urp'],default='hdrp')
args=ap.parse_args()
for folder in sorted(p for p in (ROOT/'materials').iterdir() if p.is_dir()):
    metal=np.asarray(Image.open(folder/'Metallic.png').convert('L'))
    ao=np.asarray(Image.open(folder/'AO.png').convert('L'))
    smooth=255-np.asarray(Image.open(folder/'Roughness.png').convert('L'))
    white=np.full_like(metal,255)
    target=ROOT/'exports'/('unity_'+args.pipeline)/folder.name
    target.mkdir(parents=True,exist_ok=True)
    if args.pipeline=='hdrp':
        rgba=np.stack([metal,ao,white,smooth],axis=-1)
        Image.fromarray(rgba).save(target/'HDRP_MaskMap.png')
    else:
        rgba=np.stack([metal,white,white,smooth],axis=-1)
        Image.fromarray(rgba).save(target/'URP_MetallicSmoothness.png')
        Image.fromarray(np.stack([white,ao,white],axis=-1)).save(target/'URP_Occlusion.png')
    print(target)
