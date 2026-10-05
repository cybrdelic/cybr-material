"""Deterministic, meter-scale interiors assembled from real mesh geometry.

Parquet cuts preserve the finite timber atlas scale. No picture is used as a
scene, backdrop, floor, furniture asset or lighting input.
"""
import math,random
import bpy
from mathutils import Vector
from fabrication import cushion,rug,slab_uv

ARCHITECTURES=[('14_walnut_salon','Walnut Parquet Salon'),
               ('15_oak_library','Oak Reading Room'),
               ('16_stone_kitchen','Stone Kitchen Atelier'),
               ('17_parquet_detail','Parquet / Joints and Surface Wear')]


class Interior:
    def __init__(self,api,key,title,camera,target,lens=30):
        self.a=api;self.rng=random.Random(7100+int(key[:2]));self.cache={};self.ids=set()
        self.scene=api.studio('CYCLES / Architecture / '+title,api.args.size,int(api.args.size*.75),target,camera,lens)
        self.scene.view_settings.exposure=-.45
        world=self.scene.world.node_tree;sky=world.nodes.new('ShaderNodeTexSky');sky.sky_type='NISHITA'
        sky.sun_disc=False;sky.sun_elevation=math.radians(30);sky.sun_rotation=1.3
        world.links.new(sky.outputs['Color'],world.nodes['Background'].inputs['Color'])
        world.nodes['Background'].inputs['Strength'].default_value=.08
        self.scene['scene_scale']='Meter-scale architecture; individual boards and joinery'
        self.scene['source_image_inputs']=0

    def mat(self,id):
        self.ids.add(id)
        if id not in self.cache:self.cache[id]=self.a.surface(id)
        return self.cache[id]

    def box(self,name,pos,size,id,bevel=.006,angle=0,uv_offset=(0,0)):
        mat,spec=self.mat(id)
        ob=self.a.atelier.bevel_cube(name,pos,size,mat,min(bevel,min(size)*.40))
        for face in ob.data.polygons:face.use_smooth=True
        edge=next(m for m in ob.modifiers if m.type=='BEVEL')
        edge.harden_normals=True;edge.segments=12 if bevel>=.025 else 5
        self.a.atelier.cube_uv_meters(ob,spec['tile_m'],spec.get('tile_y_m'))
        for loop in ob.data.uv_layers.active.data:loop.uv+=Vector(uv_offset)
        ob.rotation_euler.z=angle
        if spec['family']=='Wood' and any(word in name.lower() for word in ('front','face','panel','door','drawer')):
            ob.data.materials[0]=mat.copy()
            group=next(n for n in ob.data.materials[0].node_tree.nodes if n.type=='GROUP')
            group.inputs['Timber Axis Rotation'].default_value=(math.pi/2,0,0)
            ob['grain_axis']='Local Z / vertically fabricated panel'
        return ob

    def solid(self,name,pos,size,color,rough=.65,bevel=.005,metal=0):
        mat=self.a.atelier.plain(name+' / finish',color,rough,metal)
        ob=self.a.atelier.bevel_cube(name,pos,size,mat,min(bevel,min(size)*.4))
        for face in ob.data.polygons:face.use_smooth=True
        edge=next(m for m in ob.modifiers if m.type=='BEVEL')
        edge.harden_normals=True;edge.segments=12 if bevel>=.025 else 5
        return ob

    def tube(self,name,points,id,radius):
        return self.a.curve_line(name,points,self.mat(id)[0],radius)

    def vase(self,name,pos,id,height=.32,radius=.13):
        mat,spec=self.mat(id)
        profile=[(.006,0),(radius*.68,0),(radius*.85,.02*height),
                 (radius,.36*height),(radius*.92,.63*height),
                 (radius*.47,.88*height),(radius*.43,height),
                 (radius*.38,height),(radius*.37,.92*height),
                 (radius*.80,.60*height),(radius*.87,.35*height),(.006,.06*height)]
        return self.a.lathe(name,profile,mat,spec,pos,segments=64)

    def furniture_group(self,start,pos,angle):
        ca,sa=math.cos(angle),math.sin(angle)
        for ob in set(self.scene.objects)-start:
            x,y,z=ob.location
            ob.location=(pos[0]+ca*x-sa*y,pos[1]+sa*x+ca*y,pos[2]+z)
            ob.rotation_euler.z+=angle

    def chair(self,pos,angle=0):
        start=set(self.scene.objects)
        self.box('Leather lounge / steel seat cradle',(0,0,.335),(.70,.65,.032),'06_blackened_steel',.005)
        cushion(self,'Leather lounge / seat',(0,-.025,.425),(.72,.69,.15),'08_saddle_leather',.055,.014)
        cushion(self,'Leather lounge / padded back',(0,.29,.78),(.72,.16,.60),'08_saddle_leather',.048,.009)
        for x in (-.28,.28):
            self.tube('Leather lounge / back support',[(x,.26,.32),(x,.37,1.01)],'06_blackened_steel',.013)
        for x in (-.43,.43):
            self.box('Oak arm / finite cut',(x,.0,.66),(.085,.42,.07),'04_fumed_oak',.02,uv_offset=(x*.27,0))
            self.tube('Brass chair / structural bent side frame',[(x,-.30,.020),(x,-.29,.66),(x,.34,.66),(x,.34,.020)],'05_champagne_brass',.015)
        from interior_finishes import finish,cylinder
        pad=finish('Hardware / black elastomer pads',(.016,.018,.017),.83)
        for x in (-.43,.43):
            for y in (-.30,.34):cylinder(self,'Leather lounge / frame floor pad',(x,y,0),.018,.026,pad)
        for y in (-.29,.34):self.tube('Brass chair / structural under-seat tie',[(-.43,y,.335),(.43,y,.335)],'05_champagne_brass',.012)
        nickel=finish('Fittings / satin nickel',(.48,.49,.475),.24,1,True)
        for x in (-.43,.43):
            for y in (-.13,.13):cylinder(self,'Leather lounge / arm fixing screw',(x,y,.696),.004,.0015,nickel)
        self.furniture_group(start,pos,angle)

    def sofa(self,pos,angle=0):
        start=set(self.scene.objects)
        self.box('Linen sofa / upholstered base',(0,0,.36),(2.35,.94,.27),'10_natural_linen',.095)
        for x in (-.73,0,.73):
            cushion(self,'Linen sofa / seat cushion',(x,-.07,.56),(.70,.80,.19),'10_natural_linen',.066,.008)
            self.solid('Seat cushion / foam beneath woven holes',(x,-.07,.56),(.694,.794,.184),(.60,.59,.55),.9,.075)
            back=cushion(self,'Linen sofa / back cushion',(x,.30,.89),(.70,.24,.61),'10_natural_linen',.066,.003)
            core=self.solid('Back cushion / foam beneath woven holes',(x,.30,.89),(.694,.234,.604),(.60,.59,.55),.9,.087)
            # Actual welt cord follows the cushion edge rather than being painted.
            pts=[(x-.28,-.40,.625),(x+.28,-.40,.625),(x+.32,-.35,.625),
                 (x+.32,.20,.625),(x+.27,.25,.625),(x-.27,.25,.625),
                 (x-.32,.20,.625),(x-.32,-.35,.625),(x-.28,-.40,.625)]
            self.tube('Linen sofa / sewn welt',pts,'10_natural_linen',.0035)
        for x in (-1.17,1.17):self.box('Linen sofa / arm',(x,0,.67),(.22,.95,.44),'10_natural_linen',.075)
        for x in (-.98,.98):
            for y in (-.31,.31):self.box('Sofa / steel leg',(x,y,.1125),(.055,.055,.225),'06_blackened_steel',.01)
        self.furniture_group(start,pos,angle)

    def room(self,width,depth,height=3.1,window_side='left',back_opening=None):
        self.width=width;self.depth=depth;self.height=height
        ps='09_lime_plaster'
        if back_opening:
            bx,bw,bottom,top=back_opening
            for lo,hi in ((-width/2,bx-bw/2),(bx+bw/2,width/2)):
                self.box('Architecture / fireplace side wall',((lo+hi)/2,depth/2+.10,height/2),(hi-lo,.20,height),ps,.006)
            for lo,hi in ((0,bottom),(top,height)):
                self.box('Architecture / fireplace head or sill',(bx,depth/2+.10,(lo+hi)/2),(bw,.20,hi-lo),ps,.006)
        else:
            self.box('Architecture / back plaster wall',(0,depth/2+.10,height/2),(width,.20,height),ps,.006)
        self.box('Architecture / ceiling',(0,0,height+.07),(width,.0+depth,.14),ps,.006)
        # Front opening places the camera inside an actual doorway.
        self.box('Architecture / front left',( -width*.17,-depth/2-.10,height/2),(width*.66,.20,height),ps,.006)
        self.box('Architecture / door lintel',(width*.33,-depth/2-.10,height-.22),(width*.34,.20,.44),ps,.006)
        side=-1 if window_side=='left' else 1
        wall_x=side*(width/2+.10);win_y=.50;win_width=2.60;bottom=.62;top=height-.32
        self.box('Architecture / window sill wall',(wall_x,0,bottom/2),(.20,depth,bottom),ps,.006)
        self.box('Architecture / window head wall',(wall_x,0,(top+height)/2),(.20,depth,height-top),ps,.006)
        for lo,hi in ((-depth/2,win_y-win_width/2),(win_y+win_width/2,depth/2)):
            if hi>lo:self.box('Architecture / window side wall',(wall_x,(lo+hi)/2,(bottom+top)/2),(.20,hi-lo,top-bottom),ps,.006)
        self.box('Architecture / opposite plaster wall',(-wall_x,0,height/2),(.20,depth,height),ps,.006)
        for y in (win_y-win_width/2,win_y,win_y+win_width/2):
            self.box('Steel window / mullion',(side*(width/2-.005),y,(bottom+top)/2),(.085,.065,top-bottom),'06_blackened_steel',.008)
        for z in (bottom,top):self.box('Steel window / rail',(side*(width/2-.005),win_y,z),(.095,win_width+.07,.065),'06_blackened_steel',.008)
        # Open casement: daylight passes through the modeled opening.
        self.box('Travertine window / sill',(side*(width/2-.08),win_y,bottom),(.28,win_width+.14,.055),'02_roman_travertine',.012)
        for x in (-width/2+.015,width/2-.015):self.box('Architecture / lime skirting',(x,0,.075),(.028,depth,.15),ps,.003)
        self.box('Architecture / back skirting',(0,depth/2-.015,.075),(width,.028,.15),ps,.003)
        a=self.a.atelier
        a.area(self.scene,'Daylight / open casement',(side*(width/2+.15),win_y,1.86),(0,.25,.35),650,2.20,1.85,color=(1,.96,.90))
        a.area(self.scene,'Daylight / doorway',(width*.30,-depth/2+.06,2.0),(0,.60,.75),120,1.25,2.15,color=(.84,.91,1))
        data=bpy.data.lights.new('Sun / late afternoon','SUN');data.energy=1.20;data.angle=.035;data.color=(1,.93,.82)
        ob=bpy.data.objects.new(data.name,data);self.scene.collection.objects.link(ob)
        ob.location=(side*10,-4,7);a.look_at(ob,(0,.2,0))
        self.box('Courtyard / garden retaining wall',(side*(width/2+5.8),.8,.48),(.20,10,.96),ps,.008)
        self.box('Courtyard / stone paving',(side*(width/2+1.9),.6,-.02),(3.6,8.0,.04),'02_roman_travertine',.005)
        self.solid('Courtyard / continuous garden ground',(side*(width/2+12),0,-.065),(24,30,.08),(.20,.22,.14),.95,.002)
        self.plant((side*(width/2+1.45),win_y+.45,.0),ground=True)
        for y in (-1.5,2.4,4.0):self.plant((side*(width/2+4.6),y,0),ground=True)
        for y in (-1.65,2.65):
            self.box('Courtyard / pergola post',(side*(width/2+3.15),y,1.45),(.085,.085,2.9),'06_blackened_steel',.002)
        self.box('Courtyard / pergola header',(side*(width/2+3.15),.50,2.85),(.095,4.5,.16),'04_fumed_oak',.003)
        for j in range(9):self.box('Courtyard / open pergola rafter',(side*(width/2+1.65),-1.55+j*.51,2.96),(3.15,.07,.14),'04_fumed_oak',.003)
        for j in range(12):
            self.box('Courtyard / recessed paving joint',(side*(width/2+1.9),-3+j*.65,.002),(3.6,.003,.003),'06_blackened_steel',.0002)
        self.scene['room_dimensions_m']=[width,depth,height]

    def clip(self,poly):
        for axis,limit,keep_greater in ((0,-self.width/2,True),(0,self.width/2,False),(1,-self.depth/2,True),(1,self.depth/2,False)):
            if not poly:return []
            result=[];p=poly[-1];pin=(p[axis]>=limit) if keep_greater else (p[axis]<=limit)
            for q in poly:
                qin=(q[axis]>=limit) if keep_greater else (q[axis]<=limit)
                if pin!=qin:
                    t=(limit-p[axis])/(q[axis]-p[axis]);result.append((p[0]+t*(q[0]-p[0]),p[1]+t*(q[1]-p[1])))
                if qin:result.append(q)
                p,pin=q,qin
            poly=result
        return poly

    def parquet(self,id,pattern='herringbone'):
        mat,spec=self.mat(id);length=spec['tile_m']*.985;width=length/5;gap=.0008
        variants=[]
        for i in range(18):
            v=mat.copy();v.name=mat.name+' / board finish %02d'%i
            group=next(n for n in v.node_tree.nodes if n.type=='GROUP')
            tint=.88+.24*self.rng.random();group.inputs['Tint'].default_value=(tint,tint,tint,1)
            group.inputs['Roughness Offset'].default_value=self.rng.uniform(-.035,.025)
            # Contact polish follows the room's circulation strip. It is
            # evaluated from world position, independently of the atlas.
            nodes=v.node_tree.nodes;links=v.node_tree.links
            geo=nodes.new('ShaderNodeNewGeometry');sep=nodes.new('ShaderNodeSeparateXYZ')
            links.new(geo.outputs['Position'],sep.inputs[0])
            off=nodes.new('ShaderNodeMath');off.operation='SUBTRACT';off.inputs[1].default_value=.8
            links.new(sep.outputs['X'],off.inputs[0])
            sq=nodes.new('ShaderNodeMath');sq.operation='MULTIPLY';links.new(off.outputs[0],sq.inputs[0]);links.new(off.outputs[0],sq.inputs[1])
            zone=nodes.new('ShaderNodeMapRange');zone.clamp=True
            zone.inputs['From Min'].default_value=.10;zone.inputs['From Max'].default_value=2.10
            zone.inputs['To Min'].default_value=-.075;zone.inputs['To Max'].default_value=0
            links.new(sq.outputs[0],zone.inputs['Value'])
            offset=nodes.new('ShaderNodeMath');offset.operation='ADD';offset.inputs[1].default_value=group.inputs['Roughness Offset'].default_value
            links.new(zone.outputs[0],offset.inputs[0]);links.new(offset.outputs[0],group.inputs['Roughness Offset'])
            variants.append(v)
        def board(cx,cy,angle):
            ca,sa=math.cos(angle),math.sin(angle);w=width-gap;l=length-gap
            corners=[(-w/2,-l/2),(w/2,-l/2),(w/2,l/2),(-w/2,l/2)]
            poly=self.clip([(cx+ca*x-sa*y,cy+sa*x+ca*y) for x,y in corners])
            if len(poly)<3:return
            area=abs(sum(poly[i][0]*poly[(i+1)%len(poly)][1]-poly[(i+1)%len(poly)][0]*poly[i][1] for i in range(len(poly))))*.5
            if area<.0003:return
            local=[(ca*(x-cx)+sa*(y-cy),-sa*(x-cx)+ca*(y-cy)) for x,y in poly]
            h=.016+self.rng.uniform(-.00017,.00017);n=len(poly)
            verts=[(x,y,z) for z in (0,h) for x,y in local]
            faces=[tuple(reversed(range(n))),tuple(range(n,2*n))]+[(i,(i+1)%n,(i+1)%n+n,i+n) for i in range(n)]
            center=self.rng.uniform(.12,.88);flip=-1 if self.rng.random()<.5 else 1
            ob=self.a.assign_mesh('Parquet / individually cut board',verts,faces,self.rng.choice(variants),
                 lambda k,p,l:(center+verts[k][0]/spec['tile_m'],.5+flip*verts[k][1]/spec.get('tile_y_m',spec['tile_m'])))
            for f in ob.data.polygons:f.use_smooth=True
            ob.location=(cx,cy,.004);ob.rotation_euler.z=angle
            mod=ob.modifiers.new('Parquet / eased edge, 0.45 mm','BEVEL');mod.width=.00045;mod.segments=3;mod.harden_normals=True
            ob.modifiers.new('Parquet / weighted normals','WEIGHTED_NORMAL')
            ob['physical_cut_m']=[width,length];ob['joint_gap_m']=gap;ob['material_id']=id
        if pattern=='herringbone':
            theta=math.pi/4;ca,sa=math.cos(theta),math.sin(theta)
            bound=math.ceil((self.width+self.depth)/width)
            for i in range(-bound//5,bound//5+1):
                for j in range(-bound,bound+1):
                    for x,y,a in ((i*5+j+.5,i*5-j+2.5,theta),(i*5+j+3.5,i*5-j+4.5,theta+math.pi/2)):
                        px,py=x*width,y*width;cx=ca*px-sa*py;cy=sa*px+ca*py
                        if abs(cx)<self.width/2+length and abs(cy)<self.depth/2+length:board(cx,cy,a)
        else:
            for i in range(-math.ceil(self.width/length),math.ceil(self.width/length)+1):
                for j in range(-math.ceil(self.depth/length),math.ceil(self.depth/length)+1):
                    for k in range(5):
                        if (i+j)%2:board((i+.5)*length,(j+.1+k*.2)*length,math.pi/2)
                        else:board((i+.1+k*.2)*length,(j+.5)*length,0)
        self.solid('Parquet / dark joint substrate',(0,0,-.014),(self.width,self.depth,.035),(.09,.065,.04),.82)
        self.scene['floor_material']=id;self.scene['floor_pattern']=pattern
        self.scene['floor_board_length_m']=length;self.scene['floor_board_width_m']=width

    def plant(self,pos,ground=False):
        if not ground:self.vase('Porcelain / olive planter',pos,'07_bone_porcelain',.37,.21)
        bark=self.a.atelier.plain('Olive / bark',(.18,.13,.085),.86)
        leaf=self.a.atelier.plain('Olive / leaf',(.105,.16,.072),.72)
        self.a.curve_line('Olive / trunk',[(pos[0],pos[1],pos[2]+(0 if ground else .3)),(pos[0]+.025,pos[1],pos[2]+.8),(pos[0]-.04,pos[1]+.02,pos[2]+1.65)],bark,.019)
        for i in range(160 if ground else 26):
            a=self.rng.uniform(0,math.tau);z=pos[2]+self.rng.uniform(.85,1.75);r=self.rng.uniform(.16,.46)
            end=(pos[0]+math.cos(a)*r,pos[1]+math.sin(a)*r,z)
            self.a.curve_line('Olive / twig',[(pos[0],pos[1],z-.12),end],bark,.0025)
            for k in range(5):
                t=.35+k*.13;cx=pos[0]+t*(end[0]-pos[0]);cy=pos[1]+t*(end[1]-pos[1]);cz=z-.1+t*.1
                angle=a+self.rng.uniform(-1.1,1.1);length=self.rng.uniform(.045,.080);w=length*.20
                verts=[(0,0,0),(length*.4,-w,.006),(length,0,.003),(length*.4,w,.006),(length*.4,0,.012)]
                ob=self.a.assign_mesh('Olive / curved leaf',verts,[(0,1,4),(1,2,4),(2,3,4),(3,0,4)],leaf,lambda k,p,l:(0,0))
                ob.location=(cx,cy,cz);ob.rotation_euler=(self.rng.uniform(-.4,.4),self.rng.uniform(-.5,.5),angle)

    def finish(self,key):
        self.scene['material_ids']=','.join(sorted(self.ids))
        self.scene['floor_board_count']=sum(ob.name.startswith('Parquet / individually') for ob in self.scene.objects)
        self.scene['architecture_object_count']=len(self.scene.objects)
        return self.scene,key


def salon(api):
    r=Interior(api,'14_walnut_salon','Walnut Parquet Salon',(2.95,-2.55,1.62),(-.6,.70,1.02),28)
    fx=1.85
    r.room(7,6,3.1,back_opening=(fx,.96,.08,1.52));r.parquet('03_american_walnut')
    r.sofa((-.75,1.35,.02));r.chair((1.18,.45,.02),-.32)
    table=r.box('Calacatta / low coffee table',(-.52,.0,.46),(1.25,.72,.045),'01_calacatta_oro',.002)
    slab_uv(table,2.8)
    for x in (-.95,-.12):r.box('Travertine / table pier',(x,0,.23),(.22,.49,.43),'02_roman_travertine',.012)
    for x in (fx-.55,fx+.55):r.box('Travertine / fireplace jamb',(x,2.79,.82),(.18,.30,1.44),'02_roman_travertine',.018)
    r.box('Travertine / fireplace lintel',(fx,2.79,1.58),(1.28,.30,.18),'02_roman_travertine',.016)
    r.box('Steel / recessed firebox back',(fx,3.06,.81),(.96,.025,1.38),'06_blackened_steel',.003)
    for x in (fx-.47,fx+.47):r.box('Steel / recessed firebox lining',(x,2.83,.81),(.018,.47,1.38),'06_blackened_steel',.002)
    r.box('Steel / recessed firebox ceiling',(fx,2.83,1.50),(.96,.47,.018),'06_blackened_steel',.002)
    r.box('Travertine / hearth',(fx,2.55,.12),(1.52,.84,.19),'02_roman_travertine',.015)
    r.box('Brass / firebox trim',(fx,2.655,1.49),(.96,.025,.035),'05_champagne_brass',.003)
    r.vase('Porcelain / table vessel',(-.72,.10,.483),'07_bone_porcelain',.23,.08)
    r.plant((-2.78,1.10,.02))
    r.box('Steel / sculptural wall panel',(.50,2.865,1.89),(1.05,.04,.90),'06_blackened_steel',.009)
    # A woven rug leaves most of the parquet exposed; holes are real opacity.
    rug(r)
    from complete_rooms import apply
    apply(r,'14_walnut_salon')
    return r.finish('14_walnut_salon')


def library(api):
    r=Interior(api,'15_oak_library','Oak Reading Room',(2.55,-2.33,1.52),(-.60,.95,.91),29)
    r.room(6.4,5.5,3.0,'left');r.parquet('04_fumed_oak','basket')
    r.chair((-.80,.50,.02),.12);r.chair((1.18,.80,.02),-.52)
    for x in (-2.65,-1.32,0,1.32,2.65):r.box('Library / steel upright',(x,2.58,1.40),(.035,.045,2.115),'06_blackened_steel',.004)
    for x in (-2.65,2.65):
        for z in (.36,1.40,2.44):r.box('Library / recessed wall bracket',(x,2.69,z),(.03,.25,.03),'06_blackened_steel',.003)
    for z in (.36,.88,1.40,1.92,2.44):r.box('Library / steel shelf',(0,2.51,z),(5.38,.35,.035),'06_blackened_steel',.003)
    # Shelf plates are 35 mm thick. Covers sit 0.3 mm above their actual tops.
    for z in (.3778,.8978,1.4178,1.9378):
        for section in range(4):
            x=-2.49+section*1.32
            for j in range(15):
                w=r.rng.uniform(.022,.055);h=r.rng.uniform(.19,.36)
                id='08_saddle_leather' if j%4==0 else '10_natural_linen'
                r.solid('Library / paper fore-edge',(x+w/2,2.48,z+h/2),(w-.005,.234,h-.009),(.61,.57,.48),.90,.001)
                tint=.22+r.rng.random()*.68
                palette=[(.42,.48,.32),(.28,.36,.49),(.58,.24,.20),(.59,.45,.27),(.83,.75,.59),(.32,.26,.20)]
                binding_tint=tuple(c*(.72+tint*.45) for c in palette[r.rng.randrange(len(palette))])+(1,)
                for xx in (x+.0015,x+w-.0015):
                    ob=r.box('Library / cloth book cover',(xx,2.48,z+h/2),(.003,.247,h),id,.001)
                    ob.data.materials[0]=ob.data.materials[0].copy();g=next(n for n in ob.data.materials[0].node_tree.nodes if n.type=='GROUP')
                    g.inputs['Tint'].default_value=binding_tint
                ob=r.box('Library / bound spine',(x+w/2,2.352,z+h/2),(w,.014,h),id,.005)
                ob.data.materials[0]=ob.data.materials[0].copy();g=next(n for n in ob.data.materials[0].node_tree.nodes if n.type=='GROUP')
                g.inputs['Tint'].default_value=binding_tint
                # Raised spine bands and a small stamped title label give
                # bindings thickness, scale and varied plausible identities.
                for dz in (.07,h-.055):
                    r.box('Library / raised binding band',(x+w/2,2.341,z+dz),(w*.94,.004,.003),id,.0005)
                label_colors=[(.25,.12,.045),(.11,.15,.13),(.27,.20,.12),(.38,.30,.18)]
                r.solid('Library / foil spine title',(x+w/2,2.342,z+h*.67),(w*.65,.0015,.019),label_colors[j%4],.5,.0003)
                for line in range(3):r.solid('Library / title rule',(x+w/2,2.3405,z+h*.67+.005-line*.005),(w*.46,.0005,.0007),(.57,.45,.22),.48,.0001)
                x+=w+.008
            # A horizontal stack occupies the breathing space after the run.
            for stack in range(2):
                bx=x+.14;bottom=z+stack*.030
                r.solid('Library / horizontal page block',(bx,2.47,bottom+.014),(.20,.22,.025),(.64,.60,.51),.88,.001)
                for zz in (bottom+.001,bottom+.028):r.box('Library / horizontal cloth cover',(bx,2.47,zz),(.208,.232,.003),'10_natural_linen',.001)
    r.box('Steel / low cabinet carcass',(-2.42,1.72,.35),(.42,.43,.65),'06_blackened_steel',.014)
    for z in (.20,.51):r.box('Oak / finite cabinet face',(-2.42,1.495,z),(.39,.02,.285),'04_fumed_oak',.004)
    for z in (.20,.51):r.box('Brass / cabinet pull',(-2.42,1.485,z),(.12,.017,.017),'05_champagne_brass',.004)
    r.box('Calacatta / reading table',(.15,-.16,.55),(.58,.50,.035),'01_calacatta_oro',.015)
    r.box('Travertine / table support',(.15,-.16,.28),(.20,.27,.53),'02_roman_travertine',.012)
    r.vase('Porcelain / reading cup',(.18,-.13,.568),'07_bone_porcelain',.105,.045)
    r.plant((-2.42,.18,.02))
    from complete_rooms import apply
    apply(r,'15_oak_library')
    return r.finish('15_oak_library')


def kitchen(api):
    from kitchen_fabrication import build as fabricated_kitchen
    return fabricated_kitchen(api,Interior)


def floor_detail(api):
    r=Interior(api,'17_parquet_detail','Parquet / Joints and Surface Wear',(.72,-1.05,.67),(0,0,.024),60)
    r.scene.cycles.use_denoising=False
    if not api.args.draft:
        r.scene.cycles.samples=768;r.scene.cycles.adaptive_min_samples=128;r.scene.cycles.adaptive_threshold=.005
    r.scene['postprocessing']='AgX display transform only; raw path-traced floor detail'
    r.width=3.7;r.depth=3.3;r.parquet('03_american_walnut')
    # One real-height inspection board among the floor cuts exposes fine relief.
    # Routed 18 mm recess through the flooring, top flush at 20 mm.
    for ob in list(r.scene.objects):
        if not ob.name.startswith('Parquet / individually'):continue
        # A real boolean kerf leaves no overlapping wood beneath the inlay.
        cutter=r.solid('Temporary inlay router',(-.59,0,.020),(.018,4,.050),(.1,.1,.1),.8,0)
        mod=ob.modifiers.new('Parquet / routed brass recess','BOOLEAN');mod.operation='DIFFERENCE';mod.object=cutter
        bpy.context.view_layer.objects.active=ob
        bpy.ops.object.modifier_apply(modifier=mod.name)
        bpy.data.objects.remove(cutter,do_unlink=True)
    r.box('Brass / flush recessed threshold',(-.59,.0,.0185),(.0175,3.3,.003),'05_champagne_brass',.00035)
    r.box('Lime plaster / skirting',(0,.69,.10),(20,.045,.20),'09_lime_plaster',.004)
    r.solid('Floor detail / continuous neutral surround',(0,0,-.041),(200,200,.02),(.16,.15,.13),.9,.001)
    a=api.atelier
    a.area(r.scene,'Grazing window / board grain',(-1.1,-.45,.33),(0,0,.025),65,.90,1.20,color=(1,.955,.89))
    a.area(r.scene,'Soft sky / grain fill',(.60,.30,1.3),(0,0,.025),28,1.1,color=(.85,.92,1))
    a.area(r.scene,'Window / bevel reflection',(.10,.8,.60),(0,0,.025),22,.13,1.2)
    return r.finish('17_parquet_detail')


def build(api,key):
    return {'14_walnut_salon':salon,'15_oak_library':library,
            '16_stone_kitchen':kitchen,'17_parquet_detail':floor_detail}[key](api)
