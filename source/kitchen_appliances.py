"""Integrated refrigerator and sink with explicitly specified metal finishes."""
import math
import bpy
from interior_finishes import finish,box,cylinder
from fabrication import slab_uv

def refrigerator(r):
    steel=finish('Appliance / vertical brushed stainless',(.42,.455,.47),.28,1,True)
    nickel=finish('Fittings / satin nickel',(.48,.49,.475),.24,1,True)
    gasket=finish('Appliance / black door gasket',(.014,.016,.017),.82)
    x=2.47;y=2.48
    box(r,'Refrigerator / insulated enclosure',(x,y,1.085),(.85,.71,2.13),steel,.012)
    for z,h in ((1.415,1.39),(.365,.67)):
        box(r,'Refrigerator / continuous door seal',(x,2.111,z),(.841,.018,h+.008),gasket,.005)
        box(r,'Refrigerator / stainless door',(x,2.092,z),(.832,.029,h),steel,.009)
    # Mounted pulls have stand-offs, not bars touching the door skin.
    r.a.curve_line('Refrigerator / vertical nickel pull',[(x-.29,2.045,.98),(x-.29,2.045,1.68)],nickel,.012)
    for z in (1.00,1.66):box(r,'Refrigerator / handle stand-off',(x-.29,2.065,z),(.026,.038,.026),nickel,.004)
    r.a.curve_line('Refrigerator / freezer drawer pull',[(x-.26,2.045,.555),(x+.26,2.045,.555)],nickel,.012)
    for dx in (-.24,.24):box(r,'Refrigerator / freezer pull mount',(x+dx,2.065,.555),(.025,.037,.025),nickel,.004)
    box(r,'Refrigerator / recessed ventilation grille',(x,2.112,.062),(.72,.036,.048),gasket,.002)
    for i in range(12):box(r,'Refrigerator / ventilation louvre',(x-.33+i*.06,2.092,.062),(.013,.015,.036),steel,.001)
    for z in (.74,2.10):box(r,'Refrigerator / upper door hinge',(x+.393,2.107,z),(.025,.022,.045),nickel,.002)
    r.box('Kitchen / overhead refrigerator cabinet',(x,2.52,2.38),(.89,.64,.40),'04_fumed_oak',.002)
    r.box('Kitchen / overhead fitted door',(x,2.188,2.38),(.875,.020,.386),'04_fumed_oak',.001)
    r.box('Kitchen / overhead brass pull',(x,2.164,2.225),(.21,.018,.012),'05_champagne_brass',.003)
    r.scene['appliance_integration']='850 mm refrigerator, mounted nickel pulls, door seals, vent grille and fitted overhead cabinet'


def sink(r,counter):
    stainless=finish('Sink / fine brushed stainless',(.38,.415,.43),.32,1,True)
    nickel=finish('Fittings / satin nickel',(.48,.49,.475),.24,1,True)
    black=finish('Sink / recessed drain shadow',(.015,.018,.020),.65)
    x,y=.94,2.42
    cutter=box(r,'Temporary sink aperture',(x,y,.92),(.570,.424,.18),black,.055)
    mod=counter.modifiers.new('Worktop / fabricated sink cutout','BOOLEAN');mod.operation='DIFFERENCE';mod.object=cutter
    bpy.context.view_layer.objects.active=counter;bpy.ops.object.modifier_apply(modifier=mod.name)
    bpy.data.objects.remove(cutter,do_unlink=True)
    slab_uv(counter,2.8)
    profile=[(.302,.234,.951),(.287,.218,.951),(.276,.206,.936),(.245,.174,.801),(.221,.15,.785),(.002,.002,.785)]
    segments=128;verts=[]
    for width,depth,z in profile:
        for i in range(segments):
            a=i*math.tau/segments;c=math.cos(a);s=math.sin(a)
            verts.append((x+width*math.copysign(abs(c)**.32,c),y+depth*math.copysign(abs(s)**.32,s),z))
    faces=[(j*segments+i,j*segments+(i+1)%segments,(j+1)*segments+(i+1)%segments,(j+1)*segments+i) for j in range(len(profile)-1) for i in range(segments)]
    bowl=r.a.assign_mesh('Kitchen sink / formed steel basin',verts,faces,stainless,lambda k,p,l:(verts[k][0],verts[k][1]))
    solid=bowl.modifiers.new('Sink / 1 mm steel gauge','SOLIDIFY');solid.thickness=.001
    cylinder(r,'Sink / recessed drain',(x,y,.783),.034,.003,black)
    cylinder(r,'Sink / drain flange',(x,y,.786),.040,.002,stainless)
    for i in range(7):
        a=i*math.tau/7;box(r,'Sink / strainer aperture',(x+.021*math.cos(a),y+.021*math.sin(a),.7885),(.006,.006,.0005),black,.001)
    cylinder(r,'Faucet / mounted base',(x,2.716,.951),.033,.018,nickel)
    points=[(x,2.716,.960),(x,2.716,1.24)]
    for k in range(33):
        a=math.pi-k*math.pi/32
        points.append((x,2.596-.12*math.cos(a),1.24+.12*math.sin(a)))
    points.append((x,2.476,1.21))
    r.a.curve_line('Faucet / nickel gooseneck',points,nickel,.012)
    cylinder(r,'Faucet / outlet aerator',(x,2.476,1.199),.014,.014,nickel)
    box(r,'Faucet / mixer lever',(x+.044,2.716,1.035),(.055,.013,.012),nickel,.003)
    r.scene['sink_fabrication']='570 mm countertop aperture, formed 1 mm stainless bowl, drain and mounted gooseneck mixer'
