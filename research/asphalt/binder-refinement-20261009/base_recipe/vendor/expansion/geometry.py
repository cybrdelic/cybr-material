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
        # One connected slit/expanded sheet. Shared junction IDs, actual diamond holes.
        verts=[];faces=[];lookup={};uv=[];nx=7;ny=9;px=pitch*1.4;py=pitch*.9
        def vertex(x,y):
            key=(round(x,9),round(y,9))
            if key not in lookup:
                lookup[key]=len(verts);verts.append((x,y,.003));uv.append((x/s['tile_m']+.5,y/s['tile_m']+.5))
            return lookup[key]
        for j in range(ny):
            for i in range(nx):
                x=(i-(nx-1)/2)*px;y=(j-(ny-1)/2)*py
                outer=[(-px/2,-py/2),(0,-py/2),(px/2,-py/2),(px/2,0),(px/2,py/2),(0,py/2),(-px/2,py/2),(-px/2,0)]
                hx=px/2-w*1.5;hy=py/2-w*1.5
                inner=[(-hx/2,-hy/2),(0,-hy),(hx/2,-hy/2),(hx,0),(hx/2,hy/2),(0,hy),(-hx/2,hy/2),(-hx,0)]
                a=[vertex(x+dx,y+dy) for dx,dy in outer];b=[vertex(x+dx,y+dy) for dx,dy in inner]
                for k in range(8):faces.append((a[k],a[(k+1)%8],b[(k+1)%8],b[k]))
        ob=mesh('Lattice / connected expanded diamond sheet',verts,faces,mat,uv)
        solid=ob.modifiers.new('Actual sheet gauge and aperture walls','SOLIDIFY');solid.thickness=.0008
        bpy.context.view_layer.objects.active=ob;ob.select_set(True);bpy.ops.object.modifier_apply(modifier=solid.name)
        bevel=ob.modifiers.new('Machined edge rounding','BEVEL');bevel.width=.00007;bevel.segments=2
        ob['diamond_apertures']=nx*ny;ob['sheet_gauge_m']=.0008;ob['expected_connected_components']=1
        return [ob],max(nx*px,ny*py)*1.2
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
                a=math.tau*i/nt;layer=int(z/pitch);seam_angle=.55+.018*math.sin(layer*1.7)
                seam_distance=(a-seam_angle+math.pi)%math.tau-math.pi
                seam=.000024*math.exp(-(seam_distance/.033)**2)*(bead**2)
                local_radius=radius+seam
                verts.append((local_radius*math.cos(a),local_radius*math.sin(a),z+.0012))
                uv.append((a*r/s['tile_m'],z/s['tile_m']))
    stride=(nz+1)*nt
    for side in (0,1):
        for j in range(nz):
            for i in range(nt):
                a=side*stride+j*nt+i;b=side*stride+j*nt+(i+1)%nt
                face=(a,b,b+nt,a+nt);faces.append(face if side==0 else tuple(reversed(face)))
    for j in (0,):
        for i in range(nt):
            a=j*nt+i;b=j*nt+(i+1)%nt
            faces.append((a,a+stride,b+stride,b) if j==nz else (a,b,b+stride,a+stride))
    # Three perimeter roads remain visible at the open top instead of an ideal smooth rim.
    nr=24;roads=max(1,round(wall/s['bead_width_m']));rows=[]
    rows.append([nz*nt+i for i in range(nt)])
    for j in range(1,nr):
        t=j/nr;row=[]
        for i in range(nt):
            outer=verts[nz*nt+i];inner=verts[stride+nz*nt+i]
            row.append(len(verts));verts.append((outer[0]*(1-t)+inner[0]*t,outer[1]*(1-t)+inner[1]*t,
                height+.0012+pitch*.42*(math.sin(math.pi*t*roads)**2)**.65))
            uv.append((math.tau*i/nt*r/s['tile_m'],height/s['tile_m']))
        rows.append(row)
    rows.append([stride+nz*nt+i for i in range(nt)])
    for j in range(nr):
        for i in range(nt):faces.append((rows[j][i],rows[j][(i+1)%nt],rows[j+1][(i+1)%nt],rows[j+1][i]))
    ob=mesh('PETG / fused bead wall, closed thickness',verts,faces,mat,uv)
    for p in ob.data.polygons:p.use_smooth=True
    ob['layer_height_m']=pitch;ob['wall_thickness_m']=wall;ob['layer_count']=round(height/pitch)
    ob['approximation']='Continuous fused shell with layer cross sections; internal bead voids excluded'
    ob['perimeter_roads']=roads;ob['start_stop_seam_peak_m']=.000024;ob['seam_angle_rad']=.55
    # Separate printed bottom disc, correctly overlapping fused wall region.
    bpy.ops.mesh.primitive_cylinder_add(vertices=128,radius=r-wall/2,depth=.0012,location=(0,0,.0013));base=bpy.context.object
    base.name='PETG / printed solid base';base.data.materials.append(mat)
    # The studio floor top is -1 mm; ground the base rather than floating the cup.
    ob.location.z-=.0017;base.location.z-=.0017
    return [ob,base],.041

def slab(s,mat):
    import numpy as np
    # Modest grid carries actual macro displacement, aligned with source UVs.
    path=ROOT/'expansion/maps'/s['id']/'Height_Macro.png'
    n=512 if s['family'] in ('Bark','Asphalt','Concrete','Cement paste') else 140
    geometric_path=path.with_name('GeometryHeight.npy')
    if geometric_path.exists():
        a=np.load(geometric_path)[::-1].copy();n=a.shape[0]-1
    else:
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
    o['height_geometry_band']='Finite macro geometry; residual normal includes unresolved sampling band'
    return [o],width*1.2

def bark_flakes_rejected(s,mat):
    from .bark_model import flakes,sample
    import numpy as np
    objects=[];extent=s['tile_m']
    inner=plain('Bark / inner fissure periderm',(.19,.082,.026),.94)
    backing=box('Bark / inner periderm substrate',(0,0,.003),(extent,extent,.006),inner)
    objects.append(backing)
    # Closed individual flake shells: actual lifted lips and heterogeneous broken edges.
    for index,f in enumerate(flakes(s['seed'],extent)):
        nu=64;nv=24;u,v=np.meshgrid(np.linspace(0,1,nu+1),np.linspace(0,1,nv+1))
        x,y,z=sample(f,u,v);verts=list(zip(x.ravel(),y.ravel(),z.ravel()));count=len(verts)
        verts+=list(zip(x.ravel(),y.ravel(),(z-f['depth']).ravel()))
        uv=[((a+extent/2)/extent,(b+extent/2)/extent) for a,b,_ in verts];faces=[]
        for j in range(nv):
            for i in range(nu):
                q=j*(nu+1)+i;face=(q,q+1,q+nu+2,q+nu+1);faces.append(face);faces.append(tuple(reversed(tuple(a+count for a in face))))
        perimeter=list(range(nu+1))+[(j+1)*(nu+1)+nu for j in range(nv)]+[nv*(nu+1)+i for i in range(nu-1,-1,-1)]+[j*(nu+1) for j in range(nv-1,0,-1)]
        for a,b in zip(perimeter,perimeter[1:]+perimeter[:1]):faces.append((a,a+count,b+count,b))
        ob=mesh(f'Bark / fractured outer flake {index:03}',verts,faces,mat,uv)
        for p in ob.data.polygons:p.use_smooth=True
        ob['flake_thickness_m']=f['depth'];ob['lifted_lip_m']=f['lift'];ob['growth_direction']='longitudinal, locally varied'
        objects.append(ob)
    backing['bark_model']='Plated pine-type outer periderm, authored approximation; not identified measured species'
    return objects,extent*1.25

def bark(s,mat):
    objects,width=slab(s,mat);ob=objects[0];ob.name='Bark / mature fractured ridge volume r3'
    # Flatten the inward cut surface: bark does not extrude its furrows through the back.
    ob.modifiers.clear();verts=[tuple(v.co) for v in ob.data.vertices];count=len(verts)
    faces=[tuple(p.vertices) for p in ob.data.polygons];uv=[((x/s['tile_m']+.5),(y/s['tile_m']+.5)) for x,y,z in verts]
    verts +=[(x,y,.0005) for x,y,z in verts];uv+=uv.copy();faces+=[tuple(a+count for a in reversed(face)) for face in faces.copy()]
    n=512;perimeter=list(range(n+1))+[(j+1)*(n+1)+n for j in range(n)]+[n*(n+1)+i for i in range(n-1,-1,-1)]+[j*(n+1) for j in range(n-1,0,-1)]
    for a,b in zip(perimeter,perimeter[1:]+perimeter[:1]):faces.append((a,a+count,b+count,b))
    data=ob.data;new=mesh('Bark / mature fractured ridge volume r3',verts,faces,mat,uv);bpy.data.objects.remove(ob,do_unlink=True);bpy.data.meshes.remove(data)
    for p in new.data.polygons:p.use_smooth=True
    new['bark_model']='Mature Douglas-fir-inspired ridge/furrow hierarchy; r3 unreviewed'
    new['shape_reference']='NC State Extension / OSU botanical descriptions, no pixel inputs'
    # Sparse outer scales detach at furrow shoulders, never in rows of whole plates.
    import numpy as np
    def read_field(name):
        im=bpy.data.images.load(str(ROOT/'expansion/maps'/s['id']/(name+'.png')),check_existing=False)
        im.colorspace_settings.name='Non-Color';n=im.size[0];pixels=np.empty(n*n*4,np.float32);im.pixels.foreach_get(pixels);bpy.data.images.remove(im)
        return pixels.reshape(n,n,4)[...,0]
    heights=read_field('Height_Macro');furrows=read_field('FurrowMask');resolution=heights.shape[0]
    candidates=np.argwhere((furrows>.12)&(furrows<.32));rng=random.Random(s['seed']+3900);objects=[new]
    for index in range(min(18,len(candidates))):
        j,i=candidates[rng.randrange(len(candidates))];cx=(i/resolution-.5)*s['tile_m'];cy=(j/resolution-.5)*s['tile_m']
        if abs(cx)>.075 or abs(cy)>.075:continue
        angle=rng.uniform(-.4,.4);w=rng.uniform(.0025,.0065);length=rng.uniform(.004,.012);lift=rng.uniform(.00025,.0009);thickness=rng.uniform(.0002,.00055)
        vertices=[];uv=[];faces=[];nx=8;ny=12
        for v in range(ny+1):
            t=v/ny
            for u in range(nx+1):
                q=u/nx;xx=(q-.5)*w*(.65+.35*math.sin(t*math.pi));yy=(t-.5)*length
                x=cx+xx*math.cos(angle)-yy*math.sin(angle);y=cy+xx*math.sin(angle)+yy*math.cos(angle)
                ix=max(0,min(resolution-1,int((x/s['tile_m']+.5)*resolution)));iy=max(0,min(resolution-1,int((y/s['tile_m']+.5)*resolution)))
                h=.008+(heights[iy,ix]-.5)*s['height_scale_m']+lift*t**4+.00004*math.sin(q*11+t*7)
                vertices.append((x,y,float(h)));uv.append((x/s['tile_m']+.5,y/s['tile_m']+.5))
        c=len(vertices);vertices +=[(x,y,z-thickness) for x,y,z in vertices];uv+=uv.copy()
        for v in range(ny):
            for u in range(nx):
                a=v*(nx+1)+u;f=(a,a+1,a+nx+2,a+nx+1);faces.append(f);faces.append(tuple(a+c for a in reversed(f)))
        boundary=list(range(nx+1))+[(j+1)*(nx+1)+nx for j in range(ny)]+[ny*(nx+1)+i for i in range(nx-1,-1,-1)]+[j*(nx+1) for j in range(ny-1,0,-1)]
        for a,b in zip(boundary,boundary[1:]+boundary[:1]):faces.append((a,a+c,b+c,b))
        scale=mesh('Bark / sparse shoulder scale '+str(index),vertices,faces,mat,uv)
        for p in scale.data.polygons:p.use_smooth=True
        scale['curl_m']=lift;scale['thickness_m']=thickness;scale['formation']='Outer scale detachment along a primary-furrow shoulder'
        objects.append(scale)
    return objects,s['tile_m']*1.25

def textile(s,mat):
    rng=random.Random(s['seed']);t=Tubes();width=.06;fleece=s['family']=='Blanket'
    def z(x,y):
        if not fleece:return -.0002
        # Cloth swatch draped over a broad fold; ends settle near the studio ground.
        return -.0002+.009*math.exp(-(y/.014)**2)+.0015*math.sin(x/width*math.tau*2)*math.exp(-(y/.023)**2)
    n=50;verts=[];faces=[];uv=[]
    for j in range(n+1):
        for i in range(n+1):
            x=(i/n-.5)*width;y=(j/n-.5)*width;verts.append((x,y,z(x,y)));uv.append((i/n*width/s['tile_m'],j/n*width/s['tile_m']))
    for j in range(n):
        for i in range(n):q=j*(n+1)+i;faces.append((q,q+1,q+n+2,q+n+1))
    backing=mesh('Textile / woven backing',verts,faces,mat,uv);sol=backing.modifiers.new('Fabric thickness','SOLIDIFY');sol.thickness=.0008
    for p in backing.data.polygons:p.use_smooth=True
    if fleece:
        hem=Tubes()
        for edge in range(4):
            points=[]
            for j in range(81):
                q=(j/80-.5)*width
                x,y=(q,-width/2) if edge==0 else (width/2,q) if edge==1 else (q,width/2) if edge==2 else (-width/2,q)
                points.append((x,y,z(x,y)-.00015))
            hem.path(points,.00040,8)
        hem_ob=hem.object('Blanket / softly folded sewn hem',mat)
        # Quiet physical support under the draped fold, with matching broad bulge.
        support_v=[];support_f=[]
        for i in range(51):
            x=(i/50-.5)*width
            for j in range(32):
                a=j/32*math.tau
                support_v.append((x,.0045*math.cos(a),max(-.001,.0035+.0015*math.sin(x/width*math.tau*2)+.0045*math.sin(a))))
        for i in range(50):
            for j in range(32):a=i*32+j;b=i*32+(j+1)%32;support_f.append((a,b,b+32,a+32))
        support_f.append(tuple(reversed(range(32))));support_f.append(tuple(range(50*32,51*32)))
        support=mesh('Study / fleece fold support',support_v,support_f,plain('Study / soft fold support',(.13,.12,.10),.92))
    # Correlated tufts, continuous curved fibers, tapered ends. No alpha cards.
    for i in range(s['fiber_count']):
        if fleece:
            x=rng.uniform(-width/2,width/2);y=rng.uniform(-width/2,width/2)
            # Brushed nap varies continuously with cloth position, with local curling.
            direction=.45+.7*math.sin(x/width*5)+rng.uniform(-.6,.6)
        else:
            tuft=i//26;row=tuft//24;column=tuft%24
            x=(column/24-.5)*width+rng.uniform(-.00022,.00022)
            y=(row/21-.5)*width+rng.uniform(-.00022,.00022)
            direction=math.pi/2+rng.uniform(-.13,.13)
        length=s['fiber_length_m']*rng.uniform(.55,1.2);curl=rng.uniform(1.2,3.0)
        if fleece:length*=1-.23*math.exp(-(y/.008)**2)
        lean=.5 if fleece else .28;points=[];radii=[]
        segments=8 if fleece else 6
        for j in range(segments):
            a=j/(segments-1)
            if fleece:
                # Low curled arcs, varying tangents and cloth-following roots; no straight stalks.
                theta=direction+curl*a
                dx=length*.55*a*math.cos(theta);dy=length*.55*a*math.sin(theta)
                dz=length*.48*math.sin(math.pi*a*.80)+.00015*a
            else:
                dx=math.cos(direction)*length*.5*a;dy=math.sin(direction)*length*.5*a;dz=length*math.sin(math.pi*a)
            points.append((x+dx,y+dy,z(x+dx,y+dy)+dz));radii.append(s['fiber_radius_m']*(1-.60*a if fleece else 1))
        t.path(points,radii,4)
    fibers=t.object('Textile / continuous individual pile fibers',mat)
    fibers['fiber_radius_m']=s['fiber_radius_m'];fibers['pile_length_m']=s['fiber_length_m']
    # Backing yarns are visible at cut edges, separate from pile.
    backing_yarns=Tubes()
    for direction in (0,1):
        for i in range(45):
            fixed=(i/44-.5)*width;points=[]
            for j in range(61):
                q=(j/60-.5)*width;x,y=(fixed,q) if direction else (q,fixed)
                points.append((x,y,z(x,y)-.00042+(-1 if direction else 1)*.00010*math.cos(j/60*math.pi*44+i*math.pi)))
            backing_yarns.path(points,.00011,4)
    yarns=backing_yarns.object('Textile / backing yarn construction',mat)
    return [backing,fibers,yarns]+([hem_ob,support] if fleece else []),width*1.25

def grass(s,mat):
    rng=random.Random(s['seed']);verts=[];faces=[];uv=[];width=.09
    clumps=[(rng.uniform(-width*.48,width*.48),rng.uniform(-width*.48,width*.48),rng.uniform(0,math.tau)) for _ in range(max(1,s['blade_count']//12))]
    for i in range(s['blade_count']):
        cx,cy,fan=clumps[i%len(clumps)];x=cx+rng.uniform(-.00055,.00055);y=cy+rng.uniform(-.00055,.00055)
        angle=fan+rng.uniform(-1.15,1.15)
        h=s['blade_length_m']*rng.uniform(.35,1.05);w=s['blade_width_m']*rng.uniform(.5,1.1);bend=rng.uniform(.15,.65)*h
        start=len(verts)
        for j in range(7):
            a=j/6;dx=math.cos(angle)*bend*a*a;dy=math.sin(angle)*bend*a*a
            for k in (-1,0,1):
                fold=w*.18*(1-a)*(1-abs(k))
                verts.append((x+dx-math.sin(angle)*w/2*k*(1-a*.99),y+dy+math.cos(angle)*w/2*k*(1-a*.99),.002+h*a+fold));uv.append((.5+k*.5,a))
        for j in range(6):
            for k in range(2):q=start+j*3+k;faces.append((q,q+1,q+4,q+3))
    o=mesh('Grass / curved tapered blades',verts,faces,mat,uv)
    sol=o.modifiers.new('Blade thickness','SOLIDIFY');sol.thickness=.00009
    thatch=Tubes()
    for i in range(600):
        x,y,fan=clumps[i%len(clumps)];angle=fan+rng.uniform(-2,2);length=rng.uniform(.008,.025)
        thatch.path([(x+math.cos(angle)*length*j/4,y+math.sin(angle)*length*j/4,.0018+.0005*math.sin(math.pi*j/4)) for j in range(5)],.00012,4)
    dry=plain('Grass / dry basal thatch',(.26,.20,.065),.89)
    base=box('Grass / explicitly cut sod root-mat section',(0,0,.0005),(width,width,.003),plain('Soil / root-mat binder',(.06,.038,.019),.94),.0007)
    base['specimen']='Finite cut sod/root mat, not an unbounded lawn surface'
    roots=Tubes()
    for i in range(700):
        x,y,_=clumps[i%len(clumps)];angle=rng.uniform(0,math.tau)
        roots.path([(x+.004*math.cos(angle)*j/4,y+.004*math.sin(angle)*j/4,.0018-.0023*j/4) for j in range(5)],.000035,4)
    rootmat=plain('Grass / fine basal roots',(.19,.13,.063),.92)
    return [o,thatch.object('Grass / actual dry thatch',dry),roots.object('Grass / coherent fine root system',rootmat),base],width*1.4

def tiles(s,mat):
    pitch=s['tile_pitch_m'];width=pitch*3;objects=[]
    grout=plain('Tile / mineral grout',(.32,.32,.29),.87)
    objects.append(box('Tile / continuous grout bed',(0,0,.00525),(width,width,.0105),grout))
    for direction in (0,1):
        for joint in (-.5,.5):
            verts=[];faces=[]
            for j in range(21):
                along=(j/20-.5)*width
                for i in range(9):
                    across=(i/8-.5)*s['grout_width_m'];h=.0105+.0003*(2*i/8-1)**2
                    verts.append((along,joint*pitch+across,h) if direction else (joint*pitch+across,along,h))
            for j in range(20):
                for i in range(8):a=j*9+i;faces.append((a,a+1,a+10,a+9))
            objects.append(mesh('Tile / recessed grout meniscus',verts,faces,grout))
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
    builder={'Bark':bark,'Chainmail':chainmail,'Lattice':lattice,'PETG':shell,'Blanket':textile,'Carpet':textile,
             'Grass':grass,'Tile':tiles,'Engineered':engineered,'LCD':display,'CRT':display}.get(family,slab)
    objects,width=builder(s,mat)
    for o in objects:o['recipe_id']=id;o['seed']=s['seed']
    return s,objects,width
