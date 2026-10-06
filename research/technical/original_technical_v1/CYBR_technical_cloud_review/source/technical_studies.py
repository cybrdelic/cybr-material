"""Technical material study improvements. SI units, analytic inputs, no render side effects."""
import bpy, math, json, sys, argparse, hashlib
from pathlib import Path
from mathutils import Vector
sys.path.insert(0,str(Path(__file__).resolve().parent))
from expansion.geometry import plain, mesh, box, Tubes, ring, mapped, build
from expansion.catalog import get
ROOT=Path(__file__).resolve().parents[1]

def aim(o,t):o.rotation_euler=(Vector(t)-o.location).to_track_quat('-Z','Y').to_euler()

def material(name,color,rough=.4,metal=0,**kw):
    m=plain(name,color,rough,metal,**kw)
    return m

def configure(id,width,height=0,view='hero'):
    sc=bpy.context.scene;sc.render.engine='CYCLES';sc.cycles.samples=96;sc.cycles.use_denoising=True
    sc.cycles.adaptive_threshold=.025;sc.cycles.adaptive_min_samples=24;sc.cycles.max_bounces=12;sc.cycles.transmission_bounces=12
    sc.render.threads_mode='FIXED';sc.render.threads=4
    sc.render.resolution_x=960;sc.render.resolution_y=800;sc.render.resolution_percentage=100
    sc.render.image_settings.file_format='PNG';sc.render.image_settings.color_mode='RGB';sc.render.image_settings.color_depth='16'
    sc.view_settings.view_transform='AgX';sc.view_settings.look='AgX - Medium High Contrast'
    sc.world=bpy.data.worlds.new('Neutral studio');sc.world.use_nodes=True
    sc.world.node_tree.nodes['Background'].inputs[0].default_value=(.15,.18,.22,1)
    sc.world.node_tree.nodes['Background'].inputs[1].default_value=.16
    sc.unit_settings.system='METRIC';sc.unit_settings.scale_length=1
    box('Studio / floor',(0,0,-.003),(width*10,width*10,.005),plain('Floor warm gray',(.16,.15,.14),.86))
    bpy.ops.object.camera_add(location=(width*.62,-width*.88,width*.90+height));cam=bpy.context.object;cam.name='Capture / '+view
    cam.data.type='ORTHO';cam.data.ortho_scale=width*1.36;cam.data.clip_start=.0001;sc.camera=cam
    target=(0,0,height)
    if view=='detail':cam.location=(width*.27,-width*.37,width*.45+height);cam.data.ortho_scale=width*.46;target=(0,-width*.13,height)
    aim(cam,target)
    for name,pos,size,power,color in [('key',(-.6,-.55,1.1),.9,24,(1,.94,.88)),('rim',(.8,.3,.8),.65,20,(.84,.93,1)),('fill',(-.2,.8,.75),.75,4,(1,1,1))]:
        ld=bpy.data.lights.new('Capture / '+name,'AREA');ld.energy=power*width*width;ld.shape='RECTANGLE';ld.size=size*width;ld.size_y=.35*width;ld.color=color
        ob=bpy.data.objects.new(ld.name,ld);sc.collection.objects.link(ob);ob.location=(pos[0]*width,pos[1]*width,pos[2]*width+height);aim(ob,target)
    sc['recipe_id']=id;sc['review_status']='Unrendered candidate; pixel review required';sc['capture_view']=view
    return sc

def chart(i,j,nx,ny):
    """One RGB value per whole pixel; stripe partition happens only in geometry."""
    u=(i+.5)/nx;v=(j+.5)/ny
    bars=[(.82,.82,.82),(.8,.68,.025),(.015,.62,.70),(.03,.62,.055),(.58,.025,.55),(.72,.035,.02),(.025,.09,.7)]
    rgb=bars[min(6,int(u*7))] if v>.64 else ((.012,.018,.029) if v>.2 else (u*.78,)*3)
    # Fine white border and central circle show legibility beneath actual glass.
    if .22<v<.61:
        rr=((u-.5)*1.25)**2+((v-.415)*2.05)**2
        if .111<rr<.139:rgb=(.72,.79,.84)
        if abs(u-.5)<.006 or abs(v-.415)<.008:rgb=(.3,.45,.55)
    if i in (1,nx-2) or j in (1,ny-2):rgb=(.38,.48,.56)
    return rgb

def display(id):
    s=get(id);crt=id=='19_crt';pitch=s['pixel_pitch_m'];nx,ny=((76,56) if crt else (152,96));w=nx*pitch;h=ny*pitch
    z0=.004;bow=.0016 if crt else 0
    def z(x,y):return z0+bow*(1-(x/(w*.69))**2-(y/(h*.85))**2)
    backing=material('Display / opaque rear enclosure',(.011,.015,.02),.47)
    box('Display / rear enclosure',(0,0,.001),(w+.005,h+.005,.004),backing,.0012)
    # Four frame sides leave screen aperture genuinely open.
    bezel=material('Display / satin bezel',(.018,.023,.028),.32)
    thickness=s['glass_thickness_m'];edge=z(w/2,h/2);ztop=edge+thickness+.00025
    for axis,sign in [(0,-1),(0,1),(1,-1),(1,1)]:
        size=(.0019,h+.0038,ztop) if axis==0 else (w,.0019,ztop)
        pos=(sign*(w/2+.00095),0,ztop/2) if axis==0 else (0,sign*(h/2+.00095),ztop/2)
        box('Display / aperture bezel',pos,size,bezel,.00035)
    # Curved opaque black face behind physical islands, with no transmitted illumination.
    gv=[];gf=[];n=32
    for j in range(n+1):
        for i in range(n+1):
            x=(i/n-.5)*w;y=(j/n-.5)*h;gv.append((x,y,z(x,y)-.00008))
    for j in range(n):
        for i in range(n):q=j*(n+1)+i;gf.append((q,q+1,q+n+2,q+n+1))
    mesh('Display / black pixel mask',gv,gf,backing)
    # Analytic RGB islands: no image texture, resampling, baked highlight or duplicate stripe mask.
    verts=[];faces=[];indices=[];mats=[];lookup={}
    for j in range(ny):
        for i in range(nx):
            rgb=chart(i,j,nx,ny)
            for c in range(3):
                intensity=round(rgb[c]*32)/32
                key=(c,intensity)
                if key not in lookup:
                    color=[0.,0.,0.];color[c]=max(.0005,intensity)
                    ma=material('Display / '+('phosphor' if crt else 'LCD RGB')+' '+str(key),(.001,.001,.001),.7)
                    p=ma.node_tree.nodes['Principled BSDF'];p.inputs['Emission Color'].default_value=(*color,1);p.inputs['Emission Strength'].default_value=5.5 if crt else 4.5
                    lookup[key]=len(mats);mats.append(ma)
                cx=(i+.5-nx/2)*pitch+(c-1)*pitch/3;cy=(j+.5-ny/2)*pitch
                sx=pitch*.133;sy=pitch*(.43 if crt else .445);q=len(verts)
                for dx,dy in [(-sx,-sy),(sx,-sy),(sx,sy),(-sx,sy)]:verts.append((cx+dx,cy+dy,z(cx+dx,cy+dy)))
                faces.append((q,q+1,q+2,q+3));indices.append(lookup[key])
    pixels=mesh('Display / independently colored real RGB islands',verts,faces)
    for m in mats:pixels.data.materials.append(m)
    for p,idx in zip(pixels.data.polygons,indices):p.material_index=idx
    pixels['pixel_pitch_m']=pitch;pixels['pixel_grid']=[nx,ny];pixels['subpixel_count']=len(faces)
    # Separate watertight optical solid, thickness unaltered from source specification.
    verts=[];faces=[];n=48
    for side in (0,1):
        for j in range(n+1):
            for i in range(n+1):
                x=(i/n-.5)*w;y=(j/n-.5)*h;verts.append((x,y,z(x,y)+.00018+side*thickness))
    stride=(n+1)**2
    for side in (0,1):
        for j in range(n):
            for i in range(n):
                q=side*stride+j*(n+1)+i;f=(q,q+1,q+n+2,q+n+1);faces.append(f if side else tuple(reversed(f)))
    boundary=list(range(n+1))+[j*(n+1)+n for j in range(1,n+1)]+[n*(n+1)+i for i in range(n-1,-1,-1)]+[j*(n+1) for j in range(n-1,0,-1)]
    for a,b in zip(boundary,boundary[1:]+boundary[:1]):faces.append((a,b,b+stride,a+stride))
    glass=material('Display / separate coverglass IOR 1.52',(.99,.997,1),.023,transmission=1,ior=1.52)
    cover=mesh('Display / curved faceplate' if crt else 'Display / planar coverglass',verts,faces,glass)
    for p in cover.data.polygons:p.use_smooth=True
    cover['glass_thickness_m']=thickness;cover['center_bow_m']=bow;cover['separate_optical_solid']=True
    if crt:
        grille=material('CRT / metal aperture grille',(.06,.07,.08),.42,metal=.9);t=Tubes()
        for i in range(nx+1):
            x=(i-nx/2)*pitch
            t.path([(x,(j/24-.5)*h,z(x,(j/24-.5)*h)+.000025) for j in range(25)],pitch*.018,4)
        t.object('CRT / curved physical aperture grille',grille)
    return s,w*1.16,.004

def plastic(id):
    s=get(id);textured=id=='21_plastic_texture';m=material(s['name'],s['color'],.49 if textured else .25)
    p=m.node_tree.nodes['Principled BSDF'];p.inputs['Coat Weight'].default_value=.07 if textured else .17
    nd=m.node_tree.nodes;lk=m.node_tree.links
    tex=nd.new('ShaderNodeTexNoise');tex.inputs['Scale'].default_value=1850 if textured else 8500;tex.inputs['Detail'].default_value=2;tex.inputs['Roughness'].default_value=.62
    tc=nd.new('ShaderNodeTexCoord');lk.new(tc.outputs['Position'] if 'Position' in tc.outputs else tc.outputs['Object'],tex.inputs['Vector'])
    bump=nd.new('ShaderNodeBump');bump.inputs['Distance'].default_value=.000065 if textured else .000008;bump.inputs['Strength'].default_value=.5 if textured else .28;lk.new(tex.outputs['Fac'],bump.inputs['Height']);lk.new(bump.outputs[0],p.inputs['Normal'])
    ramp=nd.new('ShaderNodeMapRange');ramp.inputs['From Min'].default_value=0;ramp.inputs['From Max'].default_value=1;ramp.inputs['To Min'].default_value=.40 if textured else .225;ramp.inputs['To Max'].default_value=.56 if textured else .285
    lk.new(tex.outputs['Fac'],ramp.inputs[0]);lk.new(ramp.outputs[0],p.inputs['Roughness'])
    body=box('Plastic / curved molded shell',(0,0,.020),(.12,.08,.034),m,.013)
    body.modifiers['Physical eased edges'].segments=12
    base=box('Plastic / lower mold half',(0,0,.007),(.1204,.0804,.012),m,.0055);base.modifiers['Physical eased edges'].segments=10
    seam=box('Plastic / 0.35mm mold parting line',(0,0,.0128),(.1205,.0805,.00035),material('Mold line',tuple(v*.67 for v in s['color']),.52),.00015)
    for o in (body,base):o['texture_scale']='Object coordinates in meters, scale applied';o['micro_relief_m']=.000065 if textured else .000008
    return s,.145,.018

def petg(id):
    s=get(id);mat=mapped(s);r=.014;wall=s['wall_thickness_m'];height=.026;pitch=s['layer_height_m'];base=.0012;zbase=.0012
    nt=160;nz=round(height/pitch)*8;ni=round((height-base)/pitch)*8;verts=[];faces=[];uv=[]
    def radius(z,a,inner=False):
        phase=(z/pitch)%1;bead=math.sqrt(max(0,1-((phase-.5)/.52)**2))
        theta=math.atan2(math.sin(a+1.72),math.cos(a+1.72))
        seam=0 if inner else .000055*math.exp(-.5*(theta/.035)**2)*(.65+.35*bead)
        return r+(bead-.65)*s['bead_width_m']*.155-(wall if inner else 0)+seam
    for inner,n,lo in [(False,nz,0),(True,ni,base)]:
        for j in range(n+1):
            z=lo+(height-lo)*j/n
            for i in range(nt):
                a=math.tau*i/nt;rr=radius(z,a,inner);verts.append((rr*math.cos(a),rr*math.sin(a),z+zbase));uv.append((a*r/s['tile_m'],z/s['tile_m']))
    stride=(nz+1)*nt
    for offset,n,inner in [(0,nz,False),(stride,ni,True)]:
        for j in range(n):
            for i in range(nt):
                a=offset+j*nt+i;b=offset+j*nt+(i+1)%nt;f=(a,b,b+nt,a+nt);faces.append(tuple(reversed(f)) if inner else f)
    # The wall, floor and rim form one connected watertight optical boundary.
    for i in range(nt):
        a=nz*nt+i;b=nz*nt+(i+1)%nt;c=stride+ni*nt+i;d=stride+ni*nt+(i+1)%nt;faces.append((a,b,d,c))
    for inner,offset,z in [(False,0,zbase),(True,stride,zbase+base)]:
        center=len(verts);verts.append((0,0,z));uv.append((.5,.5))
        for i in range(nt):
            a=offset+i;b=offset+(i+1)%nt;faces.append((center,a,b) if inner else (center,b,a))
    ob=mesh('PETG / unified printed wall rim and solid floor',verts,faces,mat,uv)
    for p in ob.data.polygons:p.use_smooth=True
    ob['layer_height_m']=pitch;ob['wall_thickness_m']=wall;ob['base_thickness_m']=base;ob['layer_count']=round(height/pitch);ob['seam_peak_m']=.000055
    ob['approximation']='Continuous fused bead shell, explicit start-stop seam, internal voids excluded; nominal radial wall 1.2 mm'
    if s.get('transmission'):
        for i,c in enumerate([(.055,.065,.075),(.76,.74,.64),(.31,.047,.025)]):
            box('Optics / external background reference '+str(i),((i-1)*.007,.021,.011),(.0065,.001,.018),material('Reference stripe '+str(i),c,.65),.0001)
    return s,.045,.012

def ceramic(id):
    s=get(id);pitch=s['tile_pitch_m'];width=pitch*3
    grout=material('Ceramic / matte mineral grout',(.30,.29,.25),.88)
    clay=material('Ceramic / unglazed fired body',(.48,.39,.29),.74)
    glaze=material('Ceramic / opaque pale celadon glaze',(.58,.67,.61),.145)
    p=glaze.node_tree.nodes['Principled BSDF'];p.inputs['Coat Weight'].default_value=.55;p.inputs['Coat Roughness'].default_value=.07;p.inputs['IOR'].default_value=1.5
    tex=glaze.node_tree.nodes.new('ShaderNodeTexNoise');tex.inputs['Scale'].default_value=900;tex.inputs['Detail'].default_value=2
    coord=glaze.node_tree.nodes.new('ShaderNodeTexCoord');glaze.node_tree.links.new(coord.outputs['Object'],tex.inputs['Vector'])
    bump=glaze.node_tree.nodes.new('ShaderNodeBump');bump.inputs['Distance'].default_value=.000035;bump.inputs['Strength'].default_value=.25
    glaze.node_tree.links.new(tex.outputs['Fac'],bump.inputs['Height']);glaze.node_tree.links.new(bump.outputs[0],p.inputs['Normal'])
    box('Ceramic / recessed grout bed',(0,0,.0028),(width,width,.0056),grout,.0001)
    for j in range(3):
        for i in range(3):
            o=box('Ceramic / fired tile '+str((i,j)),((i-1)*pitch,(j-1)*pitch,.008),(pitch-.003,pitch-.003,.008),clay,.00065)
            o.data.materials.append(glaze)
            for f in o.data.polygons:f.material_index=1 if f.normal.z>.5 else 0
            o.modifiers['Physical eased edges'].segments=5;o['glaze_face_only']=True;o['slab_thickness_m']=.008
    return s,width*1.1,.006

def engineered(id):
    s=get(id);back=material('Engineered / exposed alloy substrate',(.065,.085,.105),.34,metal=.92)
    box('Engineered / separate structural backing',(0,0,.002),(.112,.112,.004),back,.0004)
    if id=='15_future_ceramic':
        m=material('Engineered / satin ceramic plates',(.15,.27,.29),.21);p=m.node_tree.nodes['Principled BSDF'];p.inputs['Coat Weight'].default_value=.3;p.inputs['Coat Roughness'].default_value=.11
        pitch=.016;r=pitch/1.8
        for j in range(-3,4):
            for i in range(-3,4):
                x=i*pitch+(j%2)*pitch/2;y=j*pitch*math.sqrt(3)/2
                if abs(x)+r>.056:continue
                bpy.ops.mesh.primitive_cylinder_add(vertices=6,radius=r,depth=.003,location=(x,y,.006),rotation=(0,0,math.pi/6));o=bpy.context.object;o.name='Engineered / discrete ceramic hex';o.data.materials.append(m)
                mod=o.modifiers.new('Rounded sintered ceramic edges','BEVEL');mod.width=.00035;mod.segments=4;o.modifiers.new('Weighted corner normals','WEIGHTED_NORMAL')
    else:
        m=material('Engineered / anodized blue aluminum',(.045,.13,.19),.27,metal=1);p=m.node_tree.nodes['Principled BSDF'];p.inputs['Anisotropic'].default_value=.36
        for i in range(36):
            box('Engineered / individually extruded cooling fin',((i-17.5)*.003,0,.0075),(.0009,.11,.007),m,.00012)
    return s,.136,.004

def chain(id):
    s=get(id);m=mapped(s);t=Tubes();pitch=s['pitch_m'];rad=s['ring_radius_m'];wire=s['wire_radius_m'];n=9;R=.14
    # Each link remains a rigid circle. A shallow cylindrical drape rotates whole links.
    def rigid_ring(cx,cy,r,axis):
        points=[];angle=cx/R;center=Vector((R*math.sin(angle),cy,.010+R*(1-math.cos(angle))))
        ex=Vector((math.cos(angle),0,math.sin(angle)));ey=Vector((0,1,0));ez=Vector((-math.sin(angle),0,math.cos(angle)))
        for k in range(64):
            a=math.tau*k/64;u=r*math.cos(a);v=r*math.sin(a)
            p=center+(ex*u+ey*v if axis=='z' else ex*u+ez*v if axis=='y' else ey*u+ez*v);points.append(tuple(p))
        t.path(points,wire,10,True)
    for j in range(n):
        for i in range(n):
            x=(i-(n-1)/2)*pitch;y=(j-(n-1)/2)*pitch;rigid_ring(x,y,rad,'z')
            if i<n-1:rigid_ring(x+pitch/2,y,pitch*.35,'y')
            if j<n-1:rigid_ring(x,y+pitch/2,pitch*.35,'x')
    o=t.object('Chainmail / rigid closed rings shallow drape',m);o['construction']='Japanese four-in-one, rigid rings tilted along cylinder';o['drape_radius_m']=R;o['link_count']=t.components
    return s,.109,.01

def main():
    a=argparse.ArgumentParser();a.add_argument('--only',nargs='+',required=True);a.add_argument('--view',default='hero',choices=['hero','detail']);args=a.parse_args(sys.argv[sys.argv.index('--')+1:])
    for id in args.only:
        bpy.ops.wm.read_factory_settings(use_empty=True)
        if id in ('18_lcd','19_crt'):s,w,h=display(id)
        elif id in ('20_plastic_smooth','21_plastic_texture'):s,w,h=plastic(id)
        elif id in ('22_petg','23_petg_transparent'):s,w,h=petg(id)
        elif id=='25_tile_ceramic':s,w,h=ceramic(id)
        elif id in ('15_future_ceramic','16_future_ribbed'):s,w,h=engineered(id)
        elif id=='14_chainmail':s,w,h=chain(id)
        elif id=='12_lattice_wire':s,objects,w=build(id);h=.004
        else:raise ValueError(id)
        sc=configure(id,w,h,args.view)
        for o in sc.objects:
            if o.type=='MESH':o['recipe_id']=id
        sc['physical_model']='Authored representative material, not measured optics';sc['optical_validation']='Cycles required for transmission; EEVEE alone not authoritative'
        path=ROOT/'scenes'/f'{id}_r4_{args.view}.blend';sc.render.filepath=str(ROOT/'evidence'/f'{id}_r4_{args.view}.png')
        for image in bpy.data.images:
            if image.source=='FILE':image.pack()
        bpy.ops.wm.save_as_mainfile(filepath=str(path))
        print('READY_SCENE',str(path),flush=True)

if __name__=='__main__':main()
