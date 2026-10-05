"""Ten individually art-directed physical specimens under shared daylight.

The view includes edges, thickness and curvature. The old flat macros remain
in the comparison archive; this suite demonstrates fabrication and response.
"""
import math
import bpy

def build(a,spec):
    key=spec['id'];family=spec['family']
    width={'Wood':.30,'Stone':.36,'Metal':.18,'Ceramic':.18,'Leather':.20,'Plaster':.32,'Textile':.22}[family]
    scene=a.studio('CYCLES / Detail / '+spec['name'],a.args.size,a.args.size,
                   (0,0,width*.105),(width*.73,-width*1.08,width*.97),55)
    scene.view_settings.exposure=-.12
    backing=a.atelier.plain('Study / neutral warm paper',(.19,.18,.16),.86)
    a.atelier.bevel_cube('Study / seamless ground',(0,0,-.011),(width*12,width*12,.02),backing,.002)
    mat,s=a.surface(key)
    def box(name,loc,size,bevel=.001):
        ob=a.atelier.bevel_cube(name,loc,size,mat,bevel)
        a.atelier.cube_uv_meters(ob,s['tile_m'],s.get('tile_y_m'));return ob
    if family=='Wood':
        box('Timber / longitudinal face and coherent end grain',(-width*.15,0,.022),(width*.47,width*.85,.044),.0014)
        ob=box('Timber / second individual saw cut',(width*.22,.006,.038),(width*.20,width*.67,.076),.0014)
        ob.rotation_euler.z=.08
    elif family=='Stone':
        stone=box('Stone / honed top and eased fabricated edge',(-.018,.014,.025),(.265,.25,.05),.0012)
        if key.startswith('02'):
            raw,rs=a.surface(key,True)
            patch=a.patch('Travertine / exposed cut face',raw,rs,(.18,.17),(-.035,.011,.0501),n=321)
        else:
            for loop in stone.data.uv_layers.active.data:loop.uv.x+=.09
        ob=box('Stone / companion seam specimen',(.118,-.065,.013),(.028,.16,.026),.0008)
    elif family=='Metal':
        # Bent 2 mm sheet, two faces, rounded bend, sheared edges.
        profile=[(-.065,0),(.038,0),(.049,.003),(.056,.011),(.059,.025),(.059,.055)]
        verts=[(x,y,z+.003) for y in (-.052,.052) for x,z in profile]
        count=len(profile);faces=[(i,i+1,i+1+count,i+count) for i in range(count-1)]
        ob=a.assign_mesh('Metal / brake-formed sheet',verts,faces,mat,lambda k,p,l:(.5+verts[k][0]/s['tile_m'],.5+verts[k][1]/s['tile_m']))
        solid=ob.modifiers.new('Fabricated 2 mm sheet','SOLIDIFY');solid.thickness=.002
        edge=a.atelier.plain('Metal / freshly sheared edge',(.38,.40,.41) if key.startswith('06') else (.63,.46,.23),.26,1)
        ob.data.materials.append(edge);solid.material_offset_rim=1
        bevel=ob.modifiers.new('Sheet / eased edge','BEVEL');bevel.width=.0003;bevel.segments=3
        a.lathe('Metal / turned collar',[(.008,0),(.024,0),(.024,.014),(.022,.016),(.009,.016),(.008,.014),(.008,0)],mat,s,(-.02,-.016,.005),False,96)
    elif family=='Ceramic':
        profile=[(.001,0),(.028,0),(.030,.003),(.030,.008),(.049,.025),(.066,.060),(.067,.075),(.065,.078),(.062,.077),(.061,.066),(.056,.047),(.042,.025),(.025,.011),(.001,.011)]
        a.lathe('Porcelain / glazed open bowl with visible wall',profile,mat,s,(0,0,0),True,128)
        bisque=a.atelier.plain('Porcelain / unglazed foot ring',(.66,.62,.52),.69)
        a.lathe('Porcelain / ground foot ring',[(.021,0),(.030,0),(.030,.003),(.021,.003),(.021,0)],bisque,s,(0,0,.0001),False,96)
    elif family=='Leather':
        raw,rs=a.surface(key,True)
        ob=a.patch('Leather / flexible 1.8 mm hide',raw,rs,(.155,.157),(0,0,.002),n=241)
        for v in ob.data.vertices:
            x,y,z=v.co
            v.co.z+=.025*math.exp(-((x-.07)/.022)**2)*((y/.157+.5)**3)
        solid=ob.modifiers.new('Hide / visible cut thickness','SOLIDIFY');solid.thickness=.0018;solid.offset=1
        thread=a.atelier.plain('Leather / flax saddle stitch',(.51,.39,.22),.8)
        for i in range(36):
            y=-.064+i*.0035
            a.curve_line('Leather / saddle stitch',[( -.066,y,.0039),(-.0657,y+.001,.0044),(-.066,y+.0024,.0039)],thread,.00022)
    elif family=='Textile':
        raw,rs=a.surface(key,True)
        ob=a.patch('Linen / resting folded cloth',raw,rs,(.185,.20),(0,0,.001),True,n=321)
        solid=ob.modifiers.new('Linen / physical thickness','SOLIDIFY');solid.thickness=.0005
        # A double-fold hem follows the actual edge coordinates.
        points=[tuple(v.co+ob.location) for v in ob.data.vertices[:321]]
        a.curve_line('Linen / sewn double hem',[(x,y,z+.0004) for x,y,z in points],mat,.00065)
    else:
        box('Plaster / troweled sample panel',(0,0,.012),(.25,.24,.024),.0008)
        raw,rs=a.surface(key,True)
        a.patch('Plaster / overlapping trowel skim',raw,rs,(.25,.24),(0,0,.024),n=321)
    a.atelier.area(scene,'Study / grazing daylight',(-width*.9,-width*.55,width*.9),(0,0,.02),width*width*115,width*.8,width*1.4,color=(1,.96,.89))
    a.atelier.area(scene,'Study / cool fill',(width,.12,width*1.4),(0,0,.02),width*width*30,width*1.3,color=(.85,.92,1))
    a.atelier.area(scene,'Study / reflection strip',(0,width,width*.9),(0,0,.02),width*width*38,width*.14,width*1.1)
    scene['inspection_width_m']=width;scene['material_ids']=key;scene['specimen']='Fabricated edge, physical thickness and surface response'
    scene['geometry_height']='Height_Macro.png';scene['normal']='Normal_Micro_OpenGL.png'
    return scene
