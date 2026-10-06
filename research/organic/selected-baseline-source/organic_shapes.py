"""Structural grass, yarn-bundle carpet, and hemmed short-nap fleece."""
import bpy,math,random
from mathutils import Vector
from expansion.geometry import mesh,plain,box,Tubes

def mat(name,rgb,rough=.85,sheen=0):
 m=plain(name,rgb,rough);p=m.node_tree.nodes.get('Principled BSDF');p.inputs['Sheen Weight'].default_value=sheen
 return m

def grass():
 rng=random.Random(2705);w=.09;v=[];f=[];mi=[]
 greens=[mat('Grass / blade '+str(i),c,.69) for i,c in enumerate(((.048,.135,.018),(.10,.245,.026),(.155,.285,.038),(.08,.19,.026),(.16,.22,.052),(.225,.25,.09)))]
 for m in greens:
  p=m.node_tree.nodes.get('Principled BSDF');col=p.inputs['Base Color'].default_value; p.inputs['Base Color'].default_value=(col[0]*.42,col[1]*.42,col[2]*.42,1);p.inputs['Roughness'].default_value=.60;p.inputs['Subsurface Weight'].default_value=.025;p.inputs['Subsurface Radius'].default_value=(.0005,.001,.0002)
 # Poa-like tuft construction with common basal crowns and blade age sequence.
 for clump in range(205):
  cx=rng.uniform(-w*.48,w*.48);cy=rng.uniform(-w*.48,w*.48);rot=rng.uniform(0,math.tau)
  for leaf in range(rng.randint(6,11)):
   age=rng.random();a=rot+leaf*2.399+rng.uniform(-.25,.25)
   h=rng.uniform(.014,.038)*(1-.28*age);bw=rng.uniform(.00065,.00155)
   x=cx+rng.uniform(-.0012,.0012);y=cy+rng.uniform(-.0012,.0012)
   lean=rng.uniform(.25,.75)+age*.36;tipz=1-age*.55;twist=rng.uniform(-.65,.65)
   q=len(v);steps=15
   for j in range(steps+1):
    t=j/steps;omt=1-t
    # Cubic arc bends continuously and old blades actually droop below their apex.
    d=h*(3*omt*t*t*lean*.50+t**3*lean)
    z=.0015+h*(3*omt*omt*t*.58+3*omt*t*t*1.1+t**3*tipz)
    width=bw*(.23+.77*math.sin(math.pi*min(1,t/.85)/2))*(1-t)**.65
    if j==steps:width=.000006
    aa=a+twist*t*t;fold=.25+.45*(1-t)
    for side in (-1,0,1):
     lateral=width*.5*side
     v.append((x+math.cos(a)*d-math.sin(aa)*lateral,y+math.sin(a)*d+math.cos(aa)*lateral,z-abs(side)*width*fold*.36))
   idx=rng.choices(range(6),weights=(1,3,2,2,1,.4))[0]
   for j in range(steps):
    for k in range(2):b=q+j*3+k;f.append((b,b+1,b+4,b+3));mi.append(idx)
 blades=mesh('Grass / folded arcing leaves in basal rosettes',v,f)
 for m in greens:blades.data.materials.append(m)
 for p,i in zip(blades.data.polygons,mi):p.material_index=i;p.use_smooth=True
 blades['clump_count']=205;blades['blade_cross_section']='folded three-point V with torsion';blades['segments_per_blade']=15
 # Basal litter is short, curled and embedded, without straight projecting straw.
 thatch=Tubes()
 for i in range(600):
  x=rng.uniform(-w*.49,w*.49);y=rng.uniform(-w*.49,w*.49);a=rng.uniform(0,math.tau);length=rng.uniform(.002,.007)
  pts=[]
  for j in range(8):
   t=j/7;aa=a+t*rng.uniform(.25,.9);pts.append((x+math.cos(aa)*length*t,y+math.sin(aa)*length*t,.0002+.0007*math.sin(math.pi*t)))
  thatch.path(pts,.00006,4)
 litter=thatch.object('Grass / embedded curled basal litter',mat('Grass / decomposing litter',(.095,.068,.026),.98))
 # Rough mineral/root plug with uneven crown instead of a smooth brown board.
 n=80;verts=[];faces=[]
 for j in range(n+1):
  for i in range(n+1):verts.append(((i/n-.5)*w,(j/n-.5)*w,-.0004+rng.uniform(-.0005,.0005)))
 for j in range(n):
  for i in range(n):q=j*(n+1)+i;faces.append((q,q+1,q+n+2,q+n+1))
 soil=mesh('Grass / root-zone granular soil',verts,faces,mat('Grass / moist organic soil',(.025,.016,.008),.98));solid=soil.modifiers.new('Root plug depth','SOLIDIFY');solid.thickness=.003
 return [blades,litter,soil],w*1.4,{'id':'27_grass','name':'Tufted curved folded grass','family':'Grass','seed':2705},'r3'

def carpet_old():
 rng=random.Random(2405);w=.06;objs=[]
 backing=mat('Carpet / woven backing',(.13,.115,.078),.98)
 objs.append(box('Carpet / flexible primary backing',(0,0,.00055),(w,w,.0011),backing,.00015))
 yarns=[mat('Carpet / heather yarn '+str(i),c,.96,.22) for i,c in enumerate(((.25,.285,.215),(.29,.315,.25),(.22,.25,.19),(.34,.35,.27)))]
 groups=[Tubes() for _ in yarns];fuzz=Tubes();count=42;pitch=w/count
 # Each loop is a three-ply twisted yarn bundle; local pressure lays the loops over.
 for j in range(count):
  for i in range(count):
   x=-w/2+(i+.5)*pitch+rng.uniform(-.00011,.00011);y=-w/2+(j+.5)*pitch+rng.uniform(-.0001,.0001)
   compression=.5*math.exp(-((x-.008)/.011)**2-((y+.008)/.020)**2)
   height=rng.uniform(.0021,.0033)*(1-compression);yaw=rng.uniform(-.28,.28);lean=.0003+.0014*compression
   path=[];segments=15
   for k in range(segments+1):
    t=k/segments;u=(t-.5)*.00088;z=.0011+height*math.sin(math.pi*t)**.76
    path.append(Vector((x+u+lean*math.sin(math.pi*t),y+math.sin(math.pi*t)*rng.uniform(-.00018,.00018),z)))
   color=rng.choices(range(4),weights=(4,3,2,1))[0]
   for ply in range(3):
    pts=[]
    for k,p in enumerate(path):
     t=k/segments;phi=math.tau*(t*2.6+ply/3);tan=(path[min(k+1,segments)]-path[max(k-1,0)]).normalized();basis=tan.cross(Vector((0,1,0))).normalized();other=tan.cross(basis).normalized()
     pts.append(tuple(p+.00016*(basis*math.cos(phi)+other*math.sin(phi))))
    groups[color].path(pts,rng.uniform(.00016,.000185),6)
   # Fine short flyaways stay on the yarn, not a field of upright spikes.
   for k in range(2):
    tt=rng.uniform(.16,.84);idx=round(tt*segments);base=path[idx];angle=rng.uniform(0,math.tau);ln=rng.uniform(.00035,.00085)
    pts=[tuple(base+Vector((math.cos(angle)*ln*t,math.sin(angle)*ln*t,ln*.42*math.sin(t*math.pi)))) for t in (0,.25,.5,.75,1)]
    fuzz.path(pts,[.000019,.000018,.000016,.000013,.000005],3)
 for i,g in enumerate(groups):objs.append(g.object('Carpet / three-ply wool loops '+str(i),yarns[i]))
 objs.append(fuzz.object('Carpet / short yarn-bound flyaways',yarns[1]))
 # Genuine crosshatched edge yarns in the backing support close examination.
 weave=Tubes()
 for i in range(90):
  d=(i/89-.5)*w
  weave.path([(-w/2,d,.0006),(w/2,d,.0006)],.0001,5)
  weave.path([(d,-w/2,.0008),(d,w/2,.0008)],.0001,5)
 objs.append(weave.object('Carpet / crosswoven primary backing strands',backing))
 for o in objs:o['construction']='three-ply loop wool; pressure-compressed region; crosswoven backing'
 return objs,w*1.25,{'id':'24_carpet','name':'Three-ply loop wool carpet','family':'Carpet','seed':2405},'r2'

def carpet():
 from transport_tubes import TransportTubes
 Tubes=TransportTubes
 rng=random.Random(2406);w=.06;objs=[]
 backing=mat('Carpet / woven backing',(.13,.115,.078),.98)
 objs.append(box('Carpet / flexible primary backing',(0,0,-.00045),(w,w,.0011),backing,.00015))
 yarns=[mat('Carpet / fine wool yarn '+str(i),c,.96,.32) for i,c in enumerate(((.25,.285,.215),(.29,.315,.25),(.22,.25,.19),(.34,.35,.27)))]
 for m in yarns:
  nodes=m.node_tree.nodes;links=m.node_tree.links
  uv=nodes.new('ShaderNodeTexCoord');mul=nodes.new('ShaderNodeVectorMath');mul.operation='MULTIPLY';mul.inputs[1].default_value=(32,2,1)
  no=nodes.new('ShaderNodeTexNoise');no.inputs['Scale'].default_value=1.5;no.inputs['Detail'].default_value=3
  bump=nodes.new('ShaderNodeBump');bump.inputs['Strength'].default_value=.15;bump.inputs['Distance'].default_value=.000008
  links.new(uv.outputs['UV'],mul.inputs[0]);links.new(mul.outputs[0],no.inputs['Vector']);links.new(no.outputs['Fac'],bump.inputs['Height']);links.new(bump.outputs[0],nodes.get('Principled BSDF').inputs['Normal'])
 groups=[Tubes() for _ in yarns];fuzz=Tubes();count=42;pitch=w/count
 for j in range(count):
  for i in range(count):
   x=-w/2+(i+.5)*pitch+rng.uniform(-.00030,.00030);y=-w/2+(j+.5)*pitch+rng.uniform(-.00030,.00030)
   compression=.45*math.exp(-((x-.008)/.011)**2-((y+.008)/.020)**2)
   height=rng.uniform(.00175,.00255)*(1-compression);yaw=rng.uniform(-.55,.55);lean=.00022+.0011*compression
   radius=rng.uniform(.00029,.000335);sway=rng.uniform(-.00009,.00009);span=rng.uniform(.00105,.00125);segments=48;path=[]
   for k in range(segments+1):
    t=k/segments;u=(t-.5)*span;v=math.sin(math.pi*t)*sway
    path.append(Vector((x+u*math.cos(yaw)-v*math.sin(yaw)+lean*math.sin(math.pi*t),y+u*math.sin(yaw)+v*math.cos(yaw),.0001+height*math.sin(math.pi*t)**.82)))
   color=rng.choices(range(4),weights=(4,3,2,1))[0];groups[color].path([tuple(p) for p in path],radius,6)
   # Fine outer filaments around a continuous spun-yarn core, not three thick ropes.
   twist=rng.uniform(2.3,3.5);phase=rng.uniform(0,math.tau)
   for fiber in range(6):
    pts=[]
    for k,pnt in enumerate(path):
     t=k/segments;phi=phase+math.tau*(t*twist+fiber/6);tan=(path[min(k+1,segments)]-path[max(k-1,0)]).normalized();basis=tan.cross(Vector((0,1,0))).normalized();other=tan.cross(basis).normalized()
     pts.append(tuple(pnt+(radius*.99)*(basis*math.cos(phi)+other*math.sin(phi))))
    groups[color].path(pts,.000016,4)
   for k in range(18):
    idx=rng.randrange(5,44);base=path[idx];angle=rng.uniform(0,math.tau);ln=rng.uniform(.00045,.00095)
    pts=[tuple(base+Vector((math.cos(angle)*ln*t,math.sin(angle)*ln*t,radius+ln*.33*math.sin(t*math.pi)))) for t in (0,.25,.5,.75,1)]
    fuzz.path(pts,[.000016,.000016,.000014,.000010,.000003],3)
 for i,g in enumerate(groups):objs.append(g.object('Carpet / spun wool loop cores and fine filaments '+str(i),yarns[i]))
 objs.append(fuzz.object('Carpet / wool flyaway halo',yarns[1]))
 weave=Tubes()
 for i in range(90):
  d=(i/89-.5)*w;weave.path([(-w/2,d,-.0004),(w/2,d,-.0004)],.0001,5);weave.path([(d,-w/2,-.0002),(d,w/2,-.0002)],.0001,5)
 objs.append(weave.object('Carpet / crosswoven backing strands',backing))
 for ob in objs:ob['construction']='smooth yarn core + six fine surface filaments; low broad loop; wool flyaways; compressed nap'
 return objs,w*1.25,{'id':'24_carpet','name':'Fine spun wool loop carpet','family':'Carpet','seed':2406},'r5'

def fleece():
 rng=random.Random(1705);w=.08;objs=[]
 wool=mat('Fleece / warm oatmeal brushed knit',(.49,.365,.24),.98,.45)
 pale=mat('Fleece / raised nap fibers',(.57,.435,.30),.97,.55)
 edge=mat('Fleece / overlock binding yarn',(.41,.28,.165),.93,.22)
 def height(x,y):
  # Gravity-supported broad fold, soft secondary ripples and a settling edge.
  return .0025+.010*math.exp(-((x+.012+.022*math.sin(y/w*2.8))/.022)**2)+.0032*math.sin((y/w+.5)*math.pi)*math.cos(x/w*5.4)**2
 n=120;v=[];f=[]
 for j in range(n+1):
  for i in range(n+1):
   x=(i/n-.5)*w;y=(j/n-.5)*w;v.append((x,y,height(x,y)))
 for j in range(n):
  for i in range(n):q=j*(n+1)+i;f.append((q,q+1,q+n+2,q+n+1))
 body=mesh('Fleece / draped interlock knit body',v,f,wool)
 for p in body.data.polygons:p.use_smooth=True
 sol=body.modifiers.new('Fleece body thickness','SOLIDIFY');sol.thickness=.00125
 objs.append(body)
 # Actual interlock-knit construction sits beneath the brushed face. It becomes
 # legible along edge turns and through the compressed nap patch.
 knit=Tubes();knitcount=100;pitch=w/knitcount
 for row in range(knitcount):
  for col in range(knitcount):
   cx=-w/2+(col+.5)*pitch;cy=-w/2+(row+.5)*pitch
   pts=[]
   for k in range(11):
    t=k/10;ang=t*math.tau
    xx=cx+math.sin(ang)*pitch*.38
    yy=cy+math.cos(ang)*pitch*.43+pitch*.20*math.sin(ang*2)
    zz=height(xx,yy)-.00130+.000055*math.cos(ang*2)
    pts.append((xx,yy,zz))
   knit.path(pts,.000075,4)
 objs.append(knit.object('Fleece / underside looped backing yarn',wool))
 # Dense low nap grows in correlated microtufts; each fiber is curved, never 5mm straw.
 nap=Tubes();count=92000
 tufts=[(rng.uniform(-w*.495,w*.495),rng.uniform(-w*.495,w*.495)) for _ in range(7000)]
 for i in range(count):
  cx,cy=tufts[i%len(tufts)]
  if i%10:cx,cy=rng.uniform(-w*.495,w*.495),rng.uniform(-w*.495,w*.495)
  x=max(-w*.498,min(w*.498,cx+rng.gauss(0,.00018)));y=max(-w*.498,min(w*.498,cy+rng.gauss(0,.00018)))
  ang=.75+.42*math.sin(x*140)+rng.gauss(0,.85);length=rng.uniform(.00155,.00225);curl=rng.uniform(.8,2.4)
  pressure=.72*math.exp(-((x-.022)/.018)**2-((y+.020)/.043)**2)
  length*=1-pressure
  z=height(x,y);pts=[]
  for k in range(6):
   t=k/5;a=ang+curl*t
   dx=math.cos(a)*length*t*.72;dy=math.sin(a)*length*t*.72
   pts.append((x+dx,y+dy,height(x+dx,y+dy)+length*.55*math.sin(t*2.2)))
  nap.path(pts,[.000029,.000032,.000031,.000027,.000020,.000006],3)
 fiber=nap.object('Fleece / short curled brushed nap',pale);fiber['fiber_count']=count;fiber['fiber_length_range_m']=[.00155,.00225];objs.append(fiber)
 # Hem has a physical rounded cord and wrapping overlock stitches on all four edges.
 hem=Tubes();stitch=Tubes()
 for side in range(4):
  pts=[]
  for i in range(121):
   t=i/120;x=(t-.5)*w;y=-w/2
   if side==1:x,y=w/2,(t-.5)*w
   elif side==2:x,y=(.5-t)*w,w/2
   elif side==3:x,y=-w/2,(.5-t)*w
   pts.append((x,y,height(x,y)-.00015))
  hem.path(pts,.00055,8)
  for i in range(90):
   t=(i+.5)/90;x=(t-.5)*w;y=-w/2
   if side==1:x,y=w/2,(t-.5)*w
   elif side==2:x,y=(.5-t)*w,w/2
   elif side==3:x,y=-w/2,(.5-t)*w
   pp=[]
   for k in range(9):
    a=k/8*math.tau;d=.00065*math.cos(a);z=height(x,y)-.0001+.0007*math.sin(a)
    pp.append((x+(d if side%2 else k/8*.0005),y+(d if side%2==0 else k/8*.0005),z))
   stitch.path(pp,.000065,5)
 objs.append(hem.object('Fleece / rounded self-fabric bound edge',wool));objs.append(stitch.object('Fleece / physical wrapping overlock stitches',edge))
 # Material microstructure is generated noise only, below resolved nap scale.
 nodes=wool.node_tree.nodes;links=wool.node_tree.links;p=nodes.get('Principled BSDF');tex=nodes.new('ShaderNodeTexNoise');tex.inputs['Scale'].default_value=2500;tex.inputs['Detail'].default_value=3
 bump=nodes.new('ShaderNodeBump');bump.inputs['Strength'].default_value=.16;bump.inputs['Distance'].default_value=.00016;links.new(tex.outputs['Fac'],bump.inputs['Height']);links.new(bump.outputs['Normal'],p.inputs['Normal'])
 for ob in objs:ob.location.z-=.0024
 return objs,w*1.25,{'id':'17_blanket','name':'Bound draped brushed fleece','family':'Blanket','seed':1705},'r5'
