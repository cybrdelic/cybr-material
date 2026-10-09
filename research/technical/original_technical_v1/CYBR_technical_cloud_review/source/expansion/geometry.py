"""Deterministic editable geometry. Dimensions in meters; no image projection."""
import math,random
import bpy
from mathutils import Vector
from pathlib import Path
from .catalog import get
ROOT=Path(__file__).resolve().parents[2]

def mesh(name,verts,faces,mat=None,uv=None):
    data=bpy.data.meshes.new(name);data.from_pydata(verts,[],faces);data.update()
    ob=bpy.data.objects.new(name,data);bpy.context.collection.objects.link(ob)
    if mat:data.materials.append(mat)
    if uv:
        layer=data.uv_layers.new(name='UVMap')
        for p in data.polygons:
            for i in p.loop_indices:layer.data[i].uv=uv[data.loops[i].vertex_index]
    return ob

def plain(name,rgb,rough=.5,metal=0,transmission=0,ior=1.5):
    m=bpy.data.materials.new(name);m.use_nodes=True;p=m.node_tree.nodes.get('Principled BSDF')
    p.inputs['Base Color'].default_value=(*rgb,1);p.inputs['Roughness'].default_value=rough
    p.inputs['Metallic'].default_value=metal;p.inputs['Transmission Weight'].default_value=transmission
    p.inputs['IOR'].default_value=ior
    return m

def mapped(s):
    m=plain(s['name'],s['color'],s['roughness'],s['metallic'],s.get('transmission',0),s['ior'])
    p=m.node_tree.nodes.get('Principled BSDF');nodes=m.node_tree.nodes;links=m.node_tree.links
    for name,socket in [('BaseColor','Base Color'),('Roughness','Roughness'),('Metallic','Metallic')]:
        path=ROOT/'expansion/maps'/s['id']/(name+'.png')
        if not path.exists():continue
        t=nodes.new('ShaderNodeTexImage');t.image=bpy.data.images.load(str(path),check_existing=True)
        t.image.colorspace_settings.name='sRGB' if name=='BaseColor' else 'Non-Color';t.label=name
        links.new(t.outputs['Color'],p.inputs[socket])
    # Finite geometry supplies PETG beads/bark macro relief. Avoid double relief.
    if s['family'] not in ('PETG','LCD','CRT'):
        name='Normal_Micro_OpenGL.png' if s['family'] in ('Bark','Concrete','Cement paste','Asphalt','Plastic') else 'Normal_OpenGL.png'
        path=ROOT/'expansion/maps'/s['id']/name
        if path.exists():
            t=nodes.new('ShaderNodeTexImage');t.image=bpy.data.images.load(str(path),check_existing=True);t.image.colorspace_settings.name='Non-Color'
            normal=nodes.new('ShaderNodeNormalMap');links.new(t.outputs['Color'],normal.inputs['Color']);links.new(normal.outputs[0],p.inputs['Normal'])
    p.inputs['Anisotropic'].default_value=s.get('anisotropy',0)
    p.inputs['Coat Weight'].default_value=s.get('coat',0)
    p.inputs['Sheen Weight'].default_value=s.get('sheen',0)
    if s['family']=='PETG':
        # UV tangent follows the extruded circumference.
        tangent=nodes.new('ShaderNodeTangent');tangent.direction_type='UV_MAP';tangent.uv_map='UVMap'
        links.new(tangent.outputs[0],p.inputs['Tangent'])
        if s.get('transmission',0):
            volume=nodes.new('ShaderNodeVolumeAbsorption');volume.inputs['Color'].default_value=(.8,.95,.98,1)
            volume.inputs['Density'].default_value=8
            links.new(volume.outputs[0],nodes.get('Material Output').inputs['Volume'])
    return m

def box(name,pos,size,mat,bevel=0):
    bpy.ops.mesh.primitive_cube_add(size=1,location=pos);o=bpy.context.object;o.name=name
    o.dimensions=size;bpy.ops.object.transform_apply(location=False,rotation=False,scale=True)
    if mat:o.data.materials.append(mat)
    if bevel:
        mod=o.modifiers.new('Physical eased edges','BEVEL');mod.width=bevel;mod.segments=3
        o.modifiers.new('Weighted corner normals','WEIGHTED_NORMAL')
    return o

class Tubes:
    """Closed circular tubes with indexed quads; one mesh, many disjoint loops."""
    def __init__(self):self.v=[];self.f=[];self.uv=[];self.components=0
    def path(self,points,radius,sides=6,closed=False):
        start=len(self.v);n=len(points);self.components+=1
        for i,p in enumerate(points):
            p=Vector(p);before=Vector(points[(i-1)%n] if closed or i else points[i]);after=Vector(points[(i+1)%n] if closed or i<n-1 else points[i])
            tangent=(after-before).normalized();reference=Vector((0,0,1)) if abs(tangent.z)<.92 else Vector((1,0,0))
            a=tangent.cross(reference).normalized();b=tangent.cross(a).normalized()
            r=radius[i] if isinstance(radius,list) else radius
            for j in range(sides):
                self.v.append(tuple(p+r*(a*math.cos(j*math.tau/sides)+b*math.sin(j*math.tau/sides))))
                self.uv.append((j/sides,i/max(1,n-1)))
        for i in range(n if closed else n-1):
            for j in range(sides):self.f.append((start+i*sides+j,start+i*sides+(j+1)%sides,start+((i+1)%n)*sides+(j+1)%sides,start+((i+1)%n)*sides+j))
        if not closed:
            self.f.append(tuple(start+j for j in reversed(range(sides))))
            self.f.append(tuple(start+(n-1)*sides+j for j in range(sides)))
    def object(self,name,mat):
        o=mesh(name,self.v,self.f,mat,self.uv)
        for p in o.data.polygons:p.use_smooth=True
        o['procedural_components']=self.components;o['geometry_role']='physical tubes, real open voids'
        return o

def ring(t,center,radius,wire,axis='z',segments=48):
    cx,cy,cz=center;points=[]
    for i in range(segments):
        a=math.tau*i/segments;c=radius*math.cos(a);s=radius*math.sin(a)
        points.append((cx+c,cy+s,cz) if axis=='z' else (cx+c,cy,cz+s) if axis=='y' else (cx,cy+c,cz+s))
    t.path(points,wire,8,True)

def chainmail(s,mat):
    t=Tubes();pitch=s['pitch_m'];r=s['ring_radius_m'];w=s['wire_radius_m'];count=7
    # Planar large rings; perpendicular smaller links pass through both holes.
    for j in range(count):
        for i in range(count):
            x=(i-(count-1)/2)*pitch;y=(j-(count-1)/2)*pitch
            ring(t,(x,y,.004),r,w)
            if i<count-1:ring(t,(x+pitch/2,y,.004),pitch*.35,w,'y')
            if j<count-1:ring(t,(x,y+pitch/2,.004),pitch*.35,w,'x')
    o=t.object('Chainmail / closed interlocking rings',mat)
    o['construction']='Japanese 4-in-1: horizontal main rings, perpendicular bridge rings'
    o['nominal_centerline_clearance_m']=min(r-.15*pitch,.85*pitch-r)-2*w
    return [o],pitch*count*1.2

def lattice(s,mat):
    t=Tubes();pitch=s['pitch_m'];w=s['wire_radius_m'];extent=.04 if s['id']=='12_lattice_wire' else .045
    count=int(extent*2/pitch)
    if s['id']=='12_lattice_wire':
        for direction in (0,1):
            for i in range(count+1):
                fixed=-extent+i*pitch;points=[]
                for j in range(count*12+1):
                    q=-extent+j/(count*12)*extent*2
                    z=.003+(-1 if direction else 1)*w*1.08*math.cos(math.pi*(q+extent)/pitch+i*math.pi)
                    points.append((q,fixed,z) if direction else (fixed,q,z))
                t.path(points,w,8)
    else:
        # Diamond strands have free voids; nodes are joined at manufacturing junctions.
        for j in range(-4,5):
            for i in range(-4,5):
                x=i*pitch;y=j*pitch
                for dx,dy in ((pitch/2,pitch/2),(-pitch/2,pitch/2)):
                    t.path([(x,y,.003),(x+dx,y+dy,.003+w*.4)],w,4)
    return [t.object('Lattice / actual strand network',mat)],extent*2.3

def shell(s,mat):
    # Watertight printed cup shell: outer/inner wall + closed bottom + top rim.
    r=.014;wall=s['wall_thickness_m'];height=.026;pitch=s['layer_height_m'];nz=round(height/pitch)*8;nt=160
    verts=[];uv=[];faces=[]
    for side in (0,1):
        for j in range(nz+1):
            z=height*j/nz;phase=(z/pitch)%1
            bead=math.sqrt(max(0,1-((phase-.5)/.52)**2));radius=r+(bead-.65)*s['bead_width_m']*.155-side*wall
            for i in range(nt):
                a=math.tau*i/nt;verts.append((radius*math.cos(a),radius*math.sin(a),z+.0012))
                uv.append((a*r/s['tile_m'],z/s['tile_m']))
    stride=(nz+1)*nt
    for side in (0,1):
        for j in range(nz):
            for i in range(nt):
                a=side*stride+j*nt+i;b=side*stride+j*nt+(i+1)%nt
                face=(a,b,b+nt,a+nt);faces.append(face if side==0 else tuple(reversed(face)))
    for j in (0,nz):
        for i in range(nt):
            a=j*nt+i;b=j*nt+(i+1)%nt
            faces.append((a,a+stride,b+stride,b) if j==nz else (a,b,b+stride,a+stride))
    ob=mesh('PETG / fused bead wall, closed thickness',verts,faces,mat,uv)
    for p in ob.data.polygons:p.use_smooth=True
    ob['layer_height_m']=pitch;ob['wall_thickness_m']=wall;ob['layer_count']=round(height/pitch)
    ob['approximation']='Continuous fused shell with layer cross sections; internal bead voids excluded'
    # Separate printed bottom disc, correctly overlapping fused wall region.
    bpy.ops.mesh.primitive_cylinder_add(vertices=128,radius=r-wall/2,depth=.0012,location=(0,0,.0013));base=bpy.context.object
    base.name='PETG / printed solid base';base.data.materials.append(mat)
    return [ob,base],.041

def slab(s,mat):
    import numpy as np
    # Modest grid carries actual macro displacement, aligned with source UVs.
    path=ROOT/'expansion/maps'/s['id']/'Height_Macro.png'
    n=200 if s['family']=='Bark' else 140
    image=bpy.data.images.load(str(path),check_existing=False);image.colorspace_settings.name='Non-Color'
    image.scale(n+1,n+1);pixels=np.empty((n+1)*(n+1)*4,np.float32);image.pixels.foreach_get(pixels)
    a=pixels.reshape(n+1,n+1,4)[...,0];bpy.data.images.remove(image)
    width=s['tile_m'];verts=[];uv=[];faces=[]
    for j in range(n+1):
        for i in range(n+1):
            x=(i/n-.5)*width;y=(j/n-.5)*width
            verts.append((x,y,.008+(a[j,i]-.5)*s['height_scale_m']));uv.append((i/n,j/n))
    for j in range(n):
        for i in range(n):
            q=j*(n+1)+i;faces.append((q,q+1,q+n+2,q+n+1))
    o=mesh(s['family']+' / real displaced face',verts,faces,mat,uv)
    for p in o.data.polygons:p.use_smooth=True
    solid=o.modifiers.new('Solid specimen edge','SOLIDIFY');solid.thickness=.008
    o['height_geometry_band']='Height_Macro; micro normal omitted for bark proof'
    return [o],width*1.2

def textile(s,mat):
    rng=random.Random(s['seed']);t=Tubes();width=.06;fleece=s['family']=='Blanket'
    def z(x,y):return .006+.003*math.sin(x/width*math.tau)*math.cos(y/width*math.pi) if fleece else .003
    n=50;verts=[];faces=[];uv=[]
    for j in range(n+1):
        for i in range(n+1):
            x=(i/n-.5)*width;y=(j/n-.5)*width;verts.append((x,y,z(x,y)));uv.append((i/n*width/s['tile_m'],j/n*width/s['tile_m']))
    for j in range(n):
        for i in range(n):q=j*(n+1)+i;faces.append((q,q+1,q+n+2,q+n+1))
    backing=mesh('Textile / woven backing',verts,faces,mat,uv);sol=backing.modifiers.new('Fabric thickness','SOLIDIFY');sol.thickness=.0008
    # Correlated tufts, continuous curved fibers, tapered ends. No alpha cards.
    for i in range(s['fiber_count']):
        x=rng.uniform(-width/2,width/2);y=rng.uniform(-width/2,width/2)
        length=s['fiber_length_m']*rng.uniform(.55,1.2);direction=rng.uniform(0,math.tau)
        lean=.5 if fleece else .28;points=[];radii=[]
        for j in range(6):
            a=j/5
            if fleece:
                dx=math.cos(direction)*length*lean*a*a;dy=math.sin(direction)*length*lean*a*a;dz=length*a*(1-.24*a)
            else:
                dx=math.cos(direction)*length*.5*a;dy=math.sin(direction)*length*.5*a;dz=length*math.sin(math.pi*a)
            points.append((x+dx,y+dy,z(x,y)+dz));radii.append(s['fiber_radius_m']*(1-.7*a if fleece else 1))
        t.path(points,radii,4)
    fibers=t.object('Textile / continuous individual pile fibers',mat)
    fibers['fiber_radius_m']=s['fiber_radius_m'];fibers['pile_length_m']=s['fiber_length_m']
    return [backing,fibers],width*1.25

def grass(s,mat):
    rng=random.Random(s['seed']);verts=[];faces=[];uv=[];width=.09
    for i in range(s['blade_count']):
        x=rng.uniform(-width/2,width/2);y=rng.uniform(-width/2,width/2);angle=rng.uniform(0,math.tau)
        h=s['blade_length_m']*rng.uniform(.35,1.05);w=s['blade_width_m']*rng.uniform(.5,1.1);bend=rng.uniform(.15,.65)*h
        start=len(verts)
        for j in range(7):
            a=j/6;dx=math.cos(angle)*bend*a*a;dy=math.sin(angle)*bend*a*a
            for k in (-1,1):verts.append((x+dx-math.sin(angle)*w/2*k*(1-a*.99),y+dy+math.cos(angle)*w/2*k*(1-a*.99),.003+h*a));uv.append((.5+k*.5,a))
        for j in range(6):q=start+j*2;faces.append((q,q+1,q+3,q+2))
    o=mesh('Grass / curved tapered blades',verts,faces,mat,uv)
    sol=o.modifiers.new('Blade thickness','SOLIDIFY');sol.thickness=.00009
    thatch=Tubes()
    for i in range(600):
        x=rng.uniform(-width/2,width/2);y=rng.uniform(-width/2,width/2);angle=rng.uniform(0,math.tau);length=rng.uniform(.008,.025)
        thatch.path([(x+math.cos(angle)*length*j/4,y+math.sin(angle)*length*j/4,.002+.001*math.sin(math.pi*j/4)) for j in range(5)],.00018,4)
    dry=plain('Grass / dry basal thatch',(.26,.20,.065),.89)
    base=box('Grass / soil sample',(0,0,0),(width,width,.004),plain('Soil / binder',(.06,.038,.019),.94))
    return [o,thatch.object('Grass / actual dry thatch',dry),base],width*1.4

def tiles(s,mat):
    pitch=s['tile_pitch_m'];width=pitch*3;objects=[]
    grout=plain('Tile / mineral grout',(.32,.32,.29),.87)
    objects.append(box('Tile / continuous grout bed',(0,0,.003),(width,width,.006),grout))
    for j in range(3):
        for i in range(3):
            o=box(f'Tile / fired slab {i},{j}',((i-1)*pitch,(j-1)*pitch,.008),(pitch-s['grout_width_m'],pitch-s['grout_width_m'],.008),mat,.0007)
            layer=o.data.uv_layers.active
            for polygon in o.data.polygons:
                for index in polygon.loop_indices:
                    vertex=o.data.vertices[o.data.loops[index].vertex_index].co+o.location
                    layer.data[index].uv=((vertex.x+width/2)/s['tile_m'],(vertex.y+width/2)/s['tile_m'])
            objects.append(o)
    return objects,width*1.18

def engineered(s,mat):
    objects=[];width=.11
    base=box('Engineered / machined backing',(0,0,.002),(width,width,.004),mat,.0003);objects.append(base)
    if s['id']=='16_future_ribbed':
        count=int(width/s['fin_pitch_m'])
        for i in range(count):objects.append(box('Engineered / cooling fin',((i-(count-1)/2)*s['fin_pitch_m'],0,.006),(s['fin_pitch_m']*.30,width,.007),mat,.00016))
    else:
        p=s['panel_pitch_m'];radius=p/1.8
        for j in range(-3,4):
            for i in range(-3,4):
                x=i*p+(j%2)*p/2;y=j*p*math.sqrt(3)/2
                bpy.ops.mesh.primitive_cylinder_add(vertices=6,radius=radius,depth=.003,location=(x,y,.006),rotation=(0,0,math.pi/6))
                ob=bpy.context.object;ob.name='Engineered / isolated ceramic hex plate';ob.data.materials.append(mat)
                mod=ob.modifiers.new('Eased ceramic panel edges','BEVEL');mod.width=.0003;mod.segments=3
                objects.append(ob)
    return objects,width*1.25

def display(s,mat):
    width=.035;objects=[];frame=plain('Display / black dielectric mask',(.012,.012,.012),.38)
    objects.append(box('Display / opaque rear substrate',(0,0,.001),(width,width,.002),frame))
    # Actual emissive subpixel islands below a separate glass solid.
    pitch=s['pixel_pitch_m'];count=int(width/pitch);v=[];f=[];uv=[]
    for j in range(count):
        for i in range(count):
            for c in range(3):
                cx=(i+.5)*pitch-width/2+(c-1)*pitch/3;cy=(j+.5)*pitch-width/2
                sx=pitch*.115;sy=pitch*(.34 if s['family']=='LCD' else .43)
                q=len(v)
                for dx,dy in ((-sx,-sy),(sx,-sy),(sx,sy),(-sx,sy)):
                    v.append((cx+dx,cy+dy,.0021));uv.append(((cx+dx+width/2)/width,(cy+dy+width/2)/width))
                f.append((q,q+1,q+2,q+3))
    e=plain('Display / procedural emissive islands',(.025,.025,.025),.38)
    node=e.node_tree.nodes.new('ShaderNodeTexImage');node.image=bpy.data.images.load(str(ROOT/'expansion/maps'/s['id']/'Emission.png'),check_existing=True);node.image.colorspace_settings.name='Non-Color';node.interpolation='Closest'
    p=e.node_tree.nodes.get('Principled BSDF');e.node_tree.links.new(node.outputs['Color'],p.inputs['Emission Color']);p.inputs['Emission Strength'].default_value=2.4
    pixels=mesh('Display / separated RGB emissive islands',v,f,e,uv);pixels['subpixel_count']=len(f);objects.append(pixels)
    glass=plain('Display / actual solid glass',(.98,.99,1),.035,transmission=1,ior=1.52)
    thickness=s['glass_thickness_m']
    cover=box('Display / optical glass layer',(0,0,.0023+thickness/2),(width,width,thickness),glass,.00015);cover['glass_thickness_m']=thickness;objects.append(cover)
    if s['family']=='CRT':
        # Physical grille bars between phosphors and thick glass. No electronic simulation.
        grille=plain('CRT / aperture grille metal',(.06,.06,.06),.4,metal=.8)
        for i in range(count):objects.append(box('CRT / aperture grille bar',((i+.0)*pitch-width/2,0,.00225),(pitch*.07,width,.00006),grille))
    return objects,width*1.25

def build(id,parameters=None):
    s=get(id,parameters);mat=mapped(s);family=s['family']
    builder={'Chainmail':chainmail,'Lattice':lattice,'PETG':shell,'Blanket':textile,'Carpet':textile,
             'Grass':grass,'Tile':tiles,'Engineered':engineered,'LCD':display,'CRT':display}.get(family,slab)
    objects,width=builder(s,mat)
    for o in objects:o['recipe_id']=id;o['seed']=s['seed']
    return s,objects,width
