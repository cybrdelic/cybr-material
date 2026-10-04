"""Deterministic, meter-scale interiors assembled from real mesh geometry.

Parquet cuts preserve the finite timber atlas scale. No picture is used as a
scene, backdrop, floor, furniture asset or lighting input.
"""
import math,random
import bpy
from mathutils import Vector

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
        self.a.atelier.cube_uv_meters(ob,spec['tile_m'],spec.get('tile_y_m'))
        for loop in ob.data.uv_layers.active.data:loop.uv+=Vector(uv_offset)
        ob.rotation_euler.z=angle
        return ob

    def solid(self,name,pos,size,color,rough=.65,bevel=.005,metal=0):
        mat=self.a.atelier.plain(name+' / finish',color,rough,metal)
        return self.a.atelier.bevel_cube(name,pos,size,mat,min(bevel,min(size)*.4))

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
        self.box('Leather lounge / seat',(0,0,.43),(.75,.75,.20),'08_saddle_leather',.085)
        back=self.box('Leather lounge / tilted back',(0,.29,.80),(.76,.20,.66),'08_saddle_leather',.08)
        back.rotation_euler.x=math.radians(-12)
        for x in (-.43,.43):
            self.box('Oak arm / finite cut',(x,.0,.66),(.085,.42,.07),'04_fumed_oak',.02,uv_offset=(x*.27,0))
            self.tube('Brass chair / bent rail',[(x,-.30,.11),(x,-.29,.66),(x,.34,.66),(x,.34,.10)],'05_champagne_brass',.015)
        for x in (-.29,.29):
            for y in (-.27,.27):self.box('Steel chair / foot',(x,y,.1653),(.038,.038,.330),'06_blackened_steel',.008)
        for y in (-.29,.34):self.tube('Brass chair / under-seat tie',[(-.43,y,.37),(.43,y,.37)],'05_champagne_brass',.012)
        self.furniture_group(start,pos,angle)

    def sofa(self,pos,angle=0):
        start=set(self.scene.objects)
        self.box('Linen sofa / upholstered base',(0,0,.36),(2.35,.94,.27),'10_natural_linen',.095)
        for x in (-.73,0,.73):
            self.box('Linen sofa / seat cushion',(x,-.07,.56),(.70,.80,.19),'10_natural_linen',.078)
            self.solid('Seat cushion / foam beneath woven holes',(x,-.07,.56),(.694,.794,.184),(.60,.59,.55),.9,.075)
            back=self.box('Linen sofa / back cushion',(x,.30,.89),(.70,.24,.61),'10_natural_linen',.09)
            back.rotation_euler.x=math.radians(-9)
            core=self.solid('Back cushion / foam beneath woven holes',(x,.30,.89),(.694,.234,.604),(.60,.59,.55),.9,.087)
            core.rotation_euler.x=math.radians(-9)
            # Actual welt cord follows the cushion edge rather than being painted.
            pts=[(x-.28,-.40,.625),(x+.28,-.40,.625),(x+.32,-.35,.625),
                 (x+.32,.20,.625),(x+.27,.25,.625),(x-.27,.25,.625),
                 (x-.32,.20,.625),(x-.32,-.35,.625),(x-.28,-.40,.625)]
            self.tube('Linen sofa / sewn welt',pts,'10_natural_linen',.0035)
        for x in (-1.17,1.17):self.box('Linen sofa / arm',(x,0,.67),(.22,.95,.44),'10_natural_linen',.075)
        for x in (-.98,.98):
            for y in (-.31,.31):self.box('Sofa / steel leg',(x,y,.1125),(.055,.055,.225),'06_blackened_steel',.01)
        self.furniture_group(start,pos,angle)

    def room(self,width,depth,height=3.1,window_side='left'):
        self.width=width;self.depth=depth;self.height=height
        ps='09_lime_plaster'
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
        self.box('Courtyard / lime boundary',(side*(width/2+3.6),.8,1.45),(.20,8.4,2.9),ps,.008)
        self.box('Courtyard / stone paving',(side*(width/2+1.9),.6,-.02),(3.6,8.0,.04),'02_roman_travertine',.005)
        self.plant((side*(width/2+1.45),win_y+.45,.0))
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
            for f in ob.data.polygons:f.use_smooth=False
            ob.location=(cx,cy,.004);ob.rotation_euler.z=angle
            mod=ob.modifiers.new('Parquet / eased edge, 0.45 mm','BEVEL');mod.width=.00045;mod.segments=2
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

    def plant(self,pos):
        self.vase('Porcelain / olive planter',pos,'07_bone_porcelain',.37,.21)
        bark=self.a.atelier.plain('Olive / bark',(.18,.13,.085),.86)
        leaf=self.a.atelier.plain('Olive / leaf',(.105,.16,.072),.72)
        self.a.curve_line('Olive / trunk',[(pos[0],pos[1],pos[2]+.3),(pos[0]+.025,pos[1],pos[2]+.8),(pos[0]-.04,pos[1]+.02,pos[2]+1.65)],bark,.019)
        for i in range(26):
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
    r.room(7,6,3.1);r.parquet('03_american_walnut')
    r.sofa((-.75,1.35,.02));r.chair((1.18,.45,.02),-.32)
    r.box('Calacatta / low coffee table',(-.52,.0,.46),(1.25,.72,.045),'01_calacatta_oro',.020)
    for x in (-.95,-.12):r.box('Travertine / table pier',(x,0,.23),(.22,.49,.43),'02_roman_travertine',.012)
    for x in (-2.30,-1.20):r.box('Travertine / fireplace jamb',(x,2.79,.82),(.18,.30,1.44),'02_roman_travertine',.018)
    r.box('Travertine / fireplace lintel',(-1.75,2.79,1.58),(1.28,.30,.18),'02_roman_travertine',.016)
    r.box('Steel / recessed firebox back',(-1.75,2.91,.68),(.90,.025,1.02),'06_blackened_steel',.003)
    for x in (-2.20,-1.30):r.box('Steel / recessed firebox lining',(x,2.79,.68),(.018,.23,1.02),'06_blackened_steel',.002)
    r.box('Travertine / hearth',(-1.75,2.55,.12),(1.52,.84,.19),'02_roman_travertine',.015)
    r.box('Brass / firebox trim',(-1.75,2.655,1.09),(.96,.025,.035),'05_champagne_brass',.003)
    r.vase('Porcelain / table vessel',(-.72,.10,.483),'07_bone_porcelain',.23,.08)
    r.plant((-2.78,1.10,.02))
    r.box('Steel / sculptural wall panel',(.50,2.865,1.89),(1.05,.04,.90),'06_blackened_steel',.009)
    # A woven rug leaves most of the parquet exposed; holes are real opacity.
    mat,spec=r.mat('10_natural_linen');rug=api.patch('Linen / low woven rug',mat,spec,(2.35,1.55),(-.72,-.40,.024),True,n=129)
    # This room-scale sheet uses full normal; remove the inspection displacement.
    for mod in list(rug.modifiers):rug.modifiers.remove(mod)
    return r.finish('14_walnut_salon')


def library(api):
    r=Interior(api,'15_oak_library','Oak Reading Room',(2.55,-2.33,1.52),(-.60,.95,.91),29)
    r.room(6.4,5.5,3.0,'left');r.parquet('04_fumed_oak','basket')
    r.chair((-.80,.50,.02),.12);r.chair((1.18,.80,.02),-.52)
    for x in (-2.65,-1.32,0,1.32,2.65):r.box('Library / steel upright',(x,2.51,1.46),(.035,.34,2.73),'06_blackened_steel',.004)
    for z in (.36,.88,1.40,1.92,2.44):r.box('Library / steel shelf',(0,2.51,z),(5.38,.35,.035),'06_blackened_steel',.003)
    for z in (.395,.915,1.435,1.955):
        for section in range(4):
            x=-2.49+section*1.32
            for j in range(15):
                w=r.rng.uniform(.022,.055);h=r.rng.uniform(.19,.36)
                id='08_saddle_leather' if j%4==0 else '10_natural_linen'
                r.solid('Library / paper fore-edge',(x+w/2,2.48,z+h/2),(w-.005,.234,h-.009),(.61,.57,.48),.90,.001)
                tint=.36+r.rng.random()*.77
                for xx in (x+.0015,x+w-.0015):
                    ob=r.box('Library / cloth book cover',(xx,2.48,z+h/2),(.003,.247,h),id,.001)
                    ob.data.materials[0]=ob.data.materials[0].copy();g=next(n for n in ob.data.materials[0].node_tree.nodes if n.type=='GROUP')
                    g.inputs['Tint'].default_value=(tint,tint*.88,tint*.75,1)
                ob=r.box('Library / bound spine',(x+w/2,2.352,z+h/2),(w,.014,h),id,.005)
                ob.data.materials[0]=ob.data.materials[0].copy();g=next(n for n in ob.data.materials[0].node_tree.nodes if n.type=='GROUP')
                g.inputs['Tint'].default_value=(tint,tint*.88,tint*.75,1)
                x+=w+.008
    r.box('Steel / low cabinet carcass',(-2.42,1.72,.35),(.42,.43,.65),'06_blackened_steel',.014)
    for z in (.20,.51):r.box('Oak / finite cabinet face',(-2.42,1.495,z),(.39,.02,.285),'04_fumed_oak',.004)
    for z in (.20,.51):r.box('Brass / cabinet pull',(-2.42,1.485,z),(.12,.017,.017),'05_champagne_brass',.004)
    r.box('Calacatta / reading table',(.15,-.16,.55),(.58,.50,.035),'01_calacatta_oro',.015)
    r.box('Travertine / table support',(.15,-.16,.28),(.20,.27,.53),'02_roman_travertine',.012)
    r.vase('Porcelain / reading cup',(.18,-.13,.568),'07_bone_porcelain',.105,.045)
    r.plant((-2.42,.18,.02))
    return r.finish('15_oak_library')


def kitchen(api):
    r=Interior(api,'16_stone_kitchen','Stone Kitchen Atelier',(2.78,-2.45,1.59),(-.46,.77,1.00),28)
    r.room(6.8,5.8,3.1,'left')
    tile=.60
    for i in range(-6,6):
        for j in range(-5,5):
            x=(i+.5)*tile;y=(j+.5)*tile
            if abs(x)<3.4 and abs(y)<2.9:r.box('Travertine / floor tile',(x,y,.01),(tile-.002,tile-.002,.022),'02_roman_travertine',.0012,uv_offset=(i*.31,j*.19))
    r.solid('Stone floor / joint bed',(0,0,-.01),(6.8,5.8,.03),(.29,.25,.20),.90)
    for x in (-2.4,-1.6,-.8,0,.8,1.6,2.4):
        r.box('Steel kitchen / cabinet carcass',(x,2.51,.46),(.78,.62,.88),'06_blackened_steel',.009)
        for dx in (-.196,.196):
            for z in (.27,.67):r.box('Oak kitchen / finite door field',(x+dx,2.185,z),(.37,.02,.38),'04_fumed_oak',.004,uv_offset=(r.rng.uniform(-.045,.045),0))
        r.box('Brass kitchen / pull',(x,2.161,.73),(.18,.018,.018),'05_champagne_brass',.003)
    r.box('Calacatta / rear counter',(0,2.45,.936),(5.70,.80,.055),'01_calacatta_oro',.012)
    for i in range(8):
        panel=r.box('Calacatta / bookmatched backsplash panel',((i-3.5)*.72,2.865,1.36),(.718,.025,.72),'01_calacatta_oro',.0015)
        for loop in panel.data.uv_layers.active.data:
            if i%2:loop.uv.x=1-loop.uv.x
            loop.uv+=Vector(((i//2)*.137,(i//2)*.193))
    r.box('Steel kitchen / range hood',(-1.48,2.54,2.28),(1.12,.68,.62),'06_blackened_steel',.016)
    r.box('Steel kitchen / cooktop',(-1.48,2.37,.972),(1.06,.52,.022),'06_blackened_steel',.010)
    for x in (-1.73,-1.23):
        for y in (2.23,2.53):r.tube('Cooktop / actual burner ring',[(x+.10*math.cos(k*math.tau/48),y+.10*math.sin(k*math.tau/48),.999) for k in range(49)],'06_blackened_steel',.006)
    # Island is assembled from stone panels with real reveals and a waterfall edge.
    r.box('Calacatta / island top',(.20,.28,.97),(2.55,1.10,.065),'01_calacatta_oro',.018)
    for x in (-1.07,1.47):r.box('Calacatta / waterfall panel',(x,.28,.48),(.045,1.10,.95),'01_calacatta_oro',.006)
    r.box('Steel / island recessed plinth',(.20,.28,.36),(2.30,.91,.64),'06_blackened_steel',.010)
    for x in (-.63,.23,1.09):
        r.box('Walnut / island door',(x,-.282,.53),(.51,.028,.52),'03_american_walnut',.005,uv_offset=(r.rng.uniform(-.015,.015),0))
        r.box('Brass / island handle',(x,-.310,.67),(.16,.017,.017),'05_champagne_brass',.004)
    for x in (-.35,.74):
        r.box('Leather / counter stool',(x,-1.05,.70),(.44,.42,.09),'08_saddle_leather',.04)
        for dx in (-.16,.16):
            for dy in (-.14,.14):r.box('Steel / stool leg',(x+dx,-1.05+dy,.35),(.023,.023,.69),'06_blackened_steel',.004)
        r.tube('Brass / stool footrest',[(x-.16,-1.19,.25),(x+.16,-1.19,.25)],'05_champagne_brass',.009)
    r.vase('Porcelain / kitchen vessel',(.67,.34,1.003),'07_bone_porcelain',.29,.105)
    r.box('Linen / folded tea cloth',(-.31,.38,1.007),(.38,.44,.009),'10_natural_linen',.003,angle=.10)
    r.plant((-2.89,1.20,.03))
    return r.finish('16_stone_kitchen')


def floor_detail(api):
    r=Interior(api,'17_parquet_detail','Parquet / Joints and Surface Wear',(.72,-1.05,.67),(0,0,.024),60)
    r.scene.cycles.use_denoising=False
    if not api.args.draft:
        r.scene.cycles.samples=768;r.scene.cycles.adaptive_min_samples=128;r.scene.cycles.adaptive_threshold=.005
    r.scene['postprocessing']='AgX display transform only; raw path-traced floor detail'
    r.width=1.85;r.depth=1.65;r.parquet('03_american_walnut')
    # One real-height inspection board among the floor cuts exposes fine relief.
    mat,spec=r.mat('04_fumed_oak');displaced,_=api.surface('04_fumed_oak',True)
    board=api.patch('Oak / displaced finish comparison',displaced,spec,(.09,.44),(.53,.12,.0355),n=385)
    body=board.modifiers.new('Oak / actual solid batten thickness','SOLIDIFY');body.thickness=.0155;body.offset=-1
    r.box('Brass / threshold strip',(-.59,.0,.025),(.018,1.6,.018),'05_champagne_brass',.0015)
    r.box('Lime plaster / skirting',(0,.69,.10),(1.85,.045,.20),'09_lime_plaster',.004)
    a=api.atelier
    a.area(r.scene,'Grazing window / board grain',(-1.1,-.45,.33),(0,0,.025),65,.90,1.20,color=(1,.955,.89))
    a.area(r.scene,'Soft sky / grain fill',(.60,.30,1.3),(0,0,.025),28,1.1,color=(.85,.92,1))
    a.area(r.scene,'Window / bevel reflection',(.10,.8,.60),(0,0,.025),22,.13,1.2)
    return r.finish('17_parquet_detail')


def build(api,key):
    return {'14_walnut_salon':salon,'15_oak_library':library,
            '16_stone_kitchen':kitchen,'17_parquet_detail':floor_detail}[key](api)
