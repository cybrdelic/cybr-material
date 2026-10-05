"""Cabinet modules and stone fabrication in meters, with useful clearances."""
import math
from fabrication import cushion,slab_uv
from kitchen_appliances import refrigerator,sink
from complete_rooms import wall_finish
from interior_finishes import finish,box

def build(api,Interior):
    r=Interior(api,'16_stone_kitchen','Stone Kitchen Atelier',(2.78,-2.45,1.59),(-.46,.77,1.00),28)
    # Retain the original angle for controlled comparison; the main camera
    # opens the framing enough to show the refrigerator and full work run.
    original=r.scene.camera.copy();original.data=r.scene.camera.data.copy()
    original.name='Camera / v3.1.1 comparison';r.scene.collection.objects.link(original)
    r.scene['comparison_camera_name']=original.name
    r.scene.camera.location=(2.88,-2.55,1.68);r.scene.camera.data.lens=24
    api.atelier.look_at(r.scene.camera,(.30,1.0,1.05))
    r.room(6.8,5.8,3.1,'left')
    for i in range(-6,6):
        for j in range(-5,5):
            x=(i+.5)*.60;y=(j+.5)*.60
            if abs(x)<3.4 and abs(y)<2.9:r.box('Travertine / honed floor tile',(x,y,.010),(.598,.598,.022),'02_roman_travertine',.0007,uv_offset=(i*.31,j*.19))
    r.solid('Stone floor / joint bed',(0,0,-.012),(6.8,5.8,.026),(.31,.28,.23),.9)
    # Cabinet fronts sit against real 18 mm carcass sides and recessed plinth.
    for i,x in enumerate((-2.4,-1.6,-.8,0,.8,1.6)):
        r.solid('Rear cabinet / plywood enclosure',(x,2.51,.51),(.792,.60,.78),(.17,.12,.065),.64,.002)
        r.box('Rear cabinet / recessed toe kick',(x,2.56,.065),(.786,.46,.13),'06_blackened_steel',.001)
        for z,h in ((.755,.285),(.36,.49)):
            front=r.box('Oak kitchen / inset drawer front',(x,2.194,z),(.794,.022,h),'04_fumed_oak',.0015)
            r.box('Drawer / brass pull bar',(x,2.161,z+h*.28),(.25,.016,.012),'05_champagne_brass',.003)
            for dx in (-.105,.105):r.box('Drawer / handle standoff',(x+dx,2.177,z+h*.28),(.013,.026,.013),'05_champagne_brass',.002)
    counter=r.box('Calacatta / continuous rear worktop',(-.41,2.45,.93),(4.88,.78,.04),'01_calacatta_oro',.0015);slab_uv(counter,2.8)
    sink(r,counter);refrigerator(r)
    # Three slabs with narrow real butt joints; no eightfold motif stamping.
    for x in (-2.038,-.41,1.218):
        panel=r.box('Calacatta / fabricated backsplash slab',(x,2.868,1.30),(1.626,.022,.70),'01_calacatta_oro',.0005)
        slab_uv(panel,2.8,'wall')
    bronze=finish('Hardware / brushed warm bronze',(.36,.245,.10),.29,1,True)
    box(r,'Backsplash / fitted bronze cap',(-.41,2.850,1.657),(4.88,.020,.014),bronze,.001)
    for x in (-2.52,1.70):
        plate=finish('Electrical / warm satin outlet plate',(.40,.385,.35),.47,.4)
        box(r,'Kitchen / backsplash outlet plate',(x,2.847,1.145),(.073,.006,.095),plate,.003)
        for dx in (-.014,.014):box(r,'Kitchen / recessed socket',(x+dx,2.842,1.145),(.007,.003,.025),finish('Electrical / dark socket insert',(.014,.014,.012),.75),.001)
    hood=r.box('Steel kitchen / folded hood casing',(-1.48,2.56,2.26),(1.12,.64,.62),'06_blackened_steel',.002)
    r.solid('Range hood / dark recessed intake',(-1.48,2.53,1.947),(.99,.48,.014),(.014,.018,.02),.6,.001)
    for i in range(15):r.box('Range hood / intake baffle',(-1.94+i*.066,2.53,1.938),(.024,.45,.012),'06_blackened_steel',.0015)
    for x in (-1.84,-1.12):r.solid('Range hood / task light',(x,2.26,1.94),(.045,.015,.007),(.79,.76,.65),.3,.001)
    r.solid('Cooktop / black glass insert',(-1.48,2.36,.957),(1.03,.49,.012),(.009,.012,.014),.14,.003)
    for x in (-1.73,-1.23):
        for y in (2.24,2.49):
            r.box('Cooktop / burner cap',(x,y,.969),(.11,.11,.011),'06_blackened_steel',.004)
            r.tube('Cooktop / burner ring',[(x+.065*math.cos(k*math.tau/48),y+.065*math.sin(k*math.tau/48),.977) for k in range(49)],'06_blackened_steel',.005)
            for dx in (-.09,.09):r.box('Cooktop / pot support',(x+dx,y,.986),(.017,.17,.02),'06_blackened_steel',.002)
    for x in (-1.71,-1.56,-1.41,-1.26):r.box('Cooktop / rotary control',(x,2.126,.971),(.033,.029,.022),'06_blackened_steel',.007)
    # Island: coherent 18 mm carcass, continuous panel backs, working drawers,
    # 450 mm seating overhang, 145 mm recessed toe kick and 40 mm slab edges.
    r.solid('Island / plywood cabinet body',(.20,.45,.54),(2.39,.72,.78),(.12,.075,.036),.7,.001)
    r.box('Island / recessed plinth',(.20,.50,.082),(2.25,.60,.125),'06_blackened_steel',.002)
    for x in (-.596,.20,.996):
        r.box('Walnut island / fitted seating-side panel',(x,.078,.535),(.790,.024,.776),'03_american_walnut',.0013)
        for z,h in ((.76,.30),(.372,.468)):
            r.box('Walnut island / working drawer',(x,.826,z),(.790,.022,h),'03_american_walnut',.0013)
            r.box('Brass island / working-side pull',(x,.851,z+h*.27),(.26,.018,.012),'05_champagne_brass',.003)
    top=r.box('Calacatta / island fabricated slab',(.20,.24,.97),(2.55,1.24,.040),'01_calacatta_oro',.0012)
    slab_uv(top,2.8)
    for x in (-1.055,1.455):
        panel=r.box('Calacatta / mitered waterfall leg',(x,.24,.484),(.040,1.24,.930),'01_calacatta_oro',.0006)
        slab_uv(panel,2.8,'waterfall',x,.949)
    for x in (-.35,.74):
        cushion(r,'Leather / sewn counter stool',(x,-.83,.695),(.43,.39,.072),'08_saddle_leather',.025,.004)
        r.box('Counter stool / seat support',(x,-.83,.65),(.36,.33,.022),'06_blackened_steel',.002)
        for dx in (-.155,.155):
            for dy in (-.135,.135):r.box('Steel / stool leg',(x+dx,-.83+dy,.34),(.025,.025,.64),'06_blackened_steel',.002)
        for y in (-.965,-.695):r.tube('Counter stool / footrest rail',[(x-.155,y,.24),(x+.155,y,.24)],'05_champagne_brass',.007)
    r.vase('Porcelain / kitchen vessel',(.67,.38,.991),'07_bone_porcelain',.29,.105)
    linen,ls=r.mat('10_natural_linen')
    cloth=api.patch('Linen / resting tea cloth',linen,ls,(.36,.39),(-.31,.34,.992),True,n=97)
    # Fold geometry is explicit; this room-scale cloth uses the full normal
    # rather than sampling yarn displacement on a coarse fold grid.
    for mod in list(cloth.modifiers):
        if mod.type=='DISPLACE':cloth.modifiers.remove(mod)
    body=cloth.modifiers.new('Tea cloth / sewn thickness','SOLIDIFY');body.thickness=.0006
    r.plant((-2.89,1.20,.03))
    r.scene['island_knee_clearance_m']=.45
    r.scene['cabinet_reveal_m']=.004
    r.scene['stone_thickness_m']=.04
    wall_finish(r,'16_stone_kitchen')
    r.scene['room_use']='Complete kitchen with refrigerator, sink, faucet, cabinetry, cooking appliances, fitted hardware and seating'
    return r.finish('16_stone_kitchen')
