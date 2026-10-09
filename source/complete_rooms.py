"""Inhabited room use, restrained fixtures and fabrication details."""
import math
from fabrication import cushion
from interior_finishes import finish,paint,box,cylinder

def wall_finish(r,key):
    base=(.52,.495,.45) if key=='14_walnut_salon' else (.13,.175,.155) if key=='15_oak_library' else (.56,.545,.505)
    back=paint('Walls / '+('mineral sage' if key=='15_oak_library' else 'warm mineral white'),base)
    for ob in r.scene.objects:
        if ob.type=='MESH' and ob.name.startswith('Architecture / back plaster wall'):
            ob.data.materials[0]=back
    r.scene['wall_finish']='Intact mineral paint with plaster reveals and skirting'

def floor_lamp(r,pos):
    bronze=finish('Hardware / brushed warm bronze',(.36,.245,.10),.29,1,True)
    rubber=finish('Hardware / black elastomer pads',(.016,.018,.017),.83)
    paper=finish('Lighting / warm paper shade',(.62,.565,.435),.81)
    spec={'tile_m':1,'height_scale_m':0}
    cylinder(r,'Reading lamp / protected floor base',pos,.15,.018,rubber)
    cylinder(r,'Reading lamp / bronze base',(pos[0],pos[1],pos[2]+.018),.145,.012,bronze)
    cylinder(r,'Reading lamp / stem',(pos[0],pos[1],pos[2]+.030),.010,1.36,bronze)
    r.a.lathe('Reading lamp / open paper shade',[(.19,0),(.125,.26),(.122,.26),(.187,0)],paper,spec,(pos[0],pos[1],pos[2]+1.36),True,96)
    bulb=finish('Lighting / warm LED diffuser',(.88,.64,.32),.42,emission=2.8)
    cylinder(r,'Reading lamp / LED diffuser',(pos[0],pos[1],pos[2]+1.40),.040,.06,bulb)

def sideboard(r,pos):
    x,y,z=pos
    r.box('Salon sideboard / walnut enclosure',(x,y,z+.45),(1.10,.39,.64),'03_american_walnut',.002)
    for dx in (-.274,.274):
        r.box('Salon sideboard / fitted door',(x+dx,y-.202,z+.45),(.543,.018,.624),'03_american_walnut',.001)
        r.box('Salon sideboard / brass mounted pull',(x+dx,y-.224,z+.48),(.13,.016,.012),'05_champagne_brass',.003)
    for dx in (-.44,.44):
        for dy in (-.12,.12):r.box('Salon sideboard / brass foot',(x+dx,y+dy,z+.075),(.023,.023,.15),'05_champagne_brass',.002)
    for dx in (-.52,.52):
        for dz in (.27,.64):r.box('Salon sideboard / concealed hinge plate',(x+dx,y-.19,z+dz),(.013,.013,.033),'06_blackened_steel',.001)
    r.vase('Salon / table lamp ceramic base',(x-.30,y,z+.776),'07_bone_porcelain',.22,.072)
    paper=finish('Lighting / warm paper shade',(.62,.565,.435),.81)
    r.a.lathe('Salon / open table lamp shade',[(.145,0),(.105,.19),(.101,.19),(.141,0)],paper,{'tile_m':1},(x-.30,y,z+.98),True,64)

def wall_art(r):
    bronze=finish('Artwork / dark bronze frame',(.15,.095,.039),.38,1)
    canvas=paint('Artwork / chalk canvas',(.59,.54,.46))
    ochre=paint('Artwork / clay relief',(.35,.205,.105))
    x,y,z=.50,2.837,1.89
    box(r,'Salon / framed relief canvas',(x,y,z),(1.035,.024,.88),canvas,.001)
    for dx in (-.533,.533):box(r,'Salon / bronze picture frame',(x+dx,y-.006,z),(.014,.034,.91),bronze,.001)
    for dz in (-.455,.455):box(r,'Salon / bronze picture frame',(x,y-.006,z+dz),(1.08,.034,.014),bronze,.001)
    for i in range(3):
        cx=x-.22+i*.21;cz=z-.19+i*.105
        pts=[(cx+.18*math.cos(k*math.pi/48),y-.016,cz+.18*math.sin(k*math.pi/48)) for k in range(49)]
        r.a.curve_line('Salon / geometric clay relief',pts,ochre,.014)

def apply(r,key):
    wall_finish(r,key)
    if key=='14_walnut_salon':
        for ob in list(r.scene.objects):
            if ob.name.startswith('Steel / sculptural wall panel'):
                import bpy
                bpy.data.objects.remove(ob,do_unlink=True)
        wall_art(r);sideboard(r,(-2.48,2.47,.02));floor_lamp(r,(.71,1.58,.02))
        # Separate pillows have actual sewn perimeters and rest on the sofa.
        cushion(r,'Salon / loose linen cushion',(-1.47,1.46,.844),(.43,.18,.38),'10_natural_linen',.06,.004)
        r.scene['room_use']='Furnished sitting room with storage, lighting, framed artwork and sewn upholstery'
    elif key=='15_oak_library':
        floor_lamp(r,(1.83,1.12,.02))
        # A working writing console fits the window wall and has a drawer,
        # joinery rails, a reading folio and a proper task light.
        r.box('Library / writing console top',(-2.28,.73,.76),(.60,1.08,.034),'04_fumed_oak',.0015)
        for dx in (-.23,.23):
            for dy in (-.43,.43):r.box('Library / console joined leg',(-2.28+dx,.73+dy,.387),(.04,.04,.72),'04_fumed_oak',.002)
        r.box('Library / writing drawer front',(-1.968,.73,.671),(.018,.89,.11),'04_fumed_oak',.001)
        r.box('Library / drawer brass pull',(-1.95,.73,.674),(.022,.18,.012),'05_champagne_brass',.003)
        folio=finish('Books / oxblood cloth binding',(.095,.032,.025),.72)
        box(r,'Library / writing folio',(-2.24,.62,.792),(.32,.22,.027),folio,.002)
        r.vase('Library / pen cup',(-2.29,1.08,.779),'07_bone_porcelain',.11,.047)
        ink=finish('Stationery / graphite lacquer',(.021,.025,.021),.4)
        for i in range(3):cylinder(r,'Library / writing pencil',(-2.305+i*.014,1.08,.794),.0023,.17,ink)
        for x in (-2.60,2.60):
            for z in (.3778,.8978,1.4178,1.9378):
                r.box('Library / steel bookend',(x,2.48,z+.10),(.012,.22,.20),'06_blackened_steel',.001)
        r.scene['room_use']='Reading and writing room with bound books, writing console, task lighting and fitted hardware'
