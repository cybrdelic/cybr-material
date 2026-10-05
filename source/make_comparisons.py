"""Display copies only; never transforms or overwrites raw Cycles evidence."""
from pathlib import Path
from PIL import Image,ImageDraw
ROOT=Path(__file__).resolve().parents[1]

def pair(key,after,dest):
    size=(1024,768) if int(key[:2])>10 else (768,768)
    out=Image.new('RGB',(size[0]*2,size[1]+40),'#20221f');draw=ImageDraw.Draw(out)
    for i,path in enumerate((ROOT/'evidence'/f'{key}.png',after)):
        im=Image.open(path).convert('RGB');im.thumbnail(size,Image.Resampling.LANCZOS)
        out.paste(im,(i*size[0]+(size[0]-im.width)//2,40))
        draw.text((i*size[0]+15,13),'V3.1.1 BASELINE' if i==0 else 'REVISED / REAL CYCLES RENDER',fill='white')
    dest.parent.mkdir(parents=True,exist_ok=True);out.save(dest,quality=94)

if __name__=='__main__':
    import sys
    draft='--draft' in sys.argv
    for path in (ROOT/'path_traced'/('drafts' if draft else 'renders')).glob('*.png'):
        if (ROOT/'evidence'/path.name).exists():pair(path.stem,path,ROOT/'docs/images/comparisons'/f'{path.stem}.jpg')
