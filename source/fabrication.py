"""Sewn upholstery, constructed rug and coherent stone slab mapping."""
import math
import bpy

def cushion(r,name,pos,size,id,bevel=.045,compression=.006):
    ob=r.box(name,pos,size,id,bevel)
    # Apply rounded construction before localized sewn-edge gathers.
    bpy.context.view_layer.objects.active=ob;ob.select_set(True)
    edge=next(m for m in ob.modifiers if m.type=='BEVEL')
    bpy.ops.object.modifier_apply(modifier=edge.name)
    sub=ob.modifiers.new('Upholstery / tension grid','SUBSURF');sub.subdivision_type='SIMPLE';sub.levels=2
    bpy.ops.object.modifier_apply(modifier=sub.name)
    w,d,h=size
    back_panel=h>d*1.5
    for v in ob.data.vertices:
        x,y,z=v.co
        if back_panel and y<-d*.1:
            face=max(0,min(1,(-y/d-.1)/.4))
            center=math.exp(-((x/(w*.36))**4+(z/(h*.36))**4))
            gathers=math.exp(-((abs(x)-w*.38)/(w*.08))**2)
            v.co.y+=face*(compression*center+.0018*math.sin(z/.020)*gathers)
        elif not back_panel and z>h*.1:
            top=max(0,min(1,(z/h-.1)/.4))
            center=math.exp(-((x/(w*.36))**4+(y/(d*.36))**4))
            near_seam=math.exp(-((abs(x)-w*.39)/(w*.085))**2)
            gather=.0013*math.sin(y/.018)*near_seam*(1-center)
            v.co.z+=top*(-compression*center+gather)
    # Upholstered weave has a backing beneath it; pores do not see daylight.
    mat=ob.data.materials[0].copy();mat.node_tree.nodes.get('');ob.data.materials[0]=mat
    group=next(n for n in mat.node_tree.nodes if n.type=='GROUP');group.node_tree=group.node_tree.copy()
    p=next(n for n in group.node_tree.nodes if n.type=='BSDF_PRINCIPLED')
    for link in list(p.inputs['Alpha'].links):group.node_tree.links.remove(link)
    p.inputs['Alpha'].default_value=1
    points=[]
    # Welt follows the rounded panel edge: horizontal on seats, vertical on
    # backs. A dark waxed seam is readable without painting a wear blotch.
    b=min(bevel,min(size)*.40);inset=min(.020,b*.60)
    radius=math.sqrt(max(0,b*b-(b-inset)**2))+.0006
    other=h if back_panel else d
    for corner in range(4):
        sx=1 if corner in (0,3) else -1;sy=1 if corner in (0,1) else -1
        cx=sx*(w/2-b);cy=sy*(other/2-b)
        for j in range(33):
            a=(corner*90+j*90/32)*math.pi/180
            # Sweep clockwise in X/other coordinates around each corner.
            angles=(0,90,180,270)
            a=(angles[corner]+j*90/32)*math.pi/180
            u=cx+radius*math.cos(a);v=cy+radius*math.sin(a)
            if back_panel:points.append((pos[0]+u,pos[1]-d/2+inset,pos[2]+v))
            else:points.append((pos[0]+u,pos[1]+v,pos[2]+h/2-inset))
    points.append(points[0])
    if id.startswith('08'):
        from interior_finishes import finish
        seam=finish('Upholstery / dark waxed leather welt',(.041,.015,.006),.65)
        r.a.curve_line(name+' / sewn perimeter',points,seam,.0021)
    else:r.tube(name+' / sewn perimeter',points,id,.0024)
    ob['construction']='Closed upholstered cushion, tension grid, sewn perimeter, contact compression'
    return ob


def rug(r):
    mat,spec=r.mat('10_natural_linen');mat=mat.copy();g=next(n for n in mat.node_tree.nodes if n.type=='GROUP')
    g.inputs['UV Scale'].default_value=(.26,.26,1)
    width,depth=2.35,1.55;nx,ny=161,111
    verts=[]
    for j in range(ny):
        y=(j/(ny-1)-.5)*depth
        for i in range(nx):
            x=(i/(nx-1)-.5)*width
            # Flat heavy construction with one local lifted corner and a
            # compressed edge beside the chair, not global sine corrugation.
            corner=.013*math.exp(-((x-width*.47)/.12)**2-((y+depth*.46)/.16)**2)
            z=.0025+corner+.0015*math.exp(-((x-.79)/.11)**2-((y-.49)/.18)**2)
            verts.append((x,y,z))
    faces=[(j*nx+i,j*nx+i+1,(j+1)*nx+i+1,(j+1)*nx+i) for j in range(ny-1) for i in range(nx-1)]
    ob=r.a.assign_mesh('Rug / heavy woven body',verts,faces,mat,lambda k,p,l:(verts[k][0]/spec['tile_m'],verts[k][1]/spec['tile_m']))
    ob.location=(-.72,-.40,.024)
    mod=ob.modifiers.new('Rug / 5 mm woven body','SOLIDIFY');mod.thickness=.005;mod.offset=-1
    bound=r.a.atelier.plain('Rug / bound flax selvedge',(.38,.32,.24),.9)
    perimeter=list(range(nx))+[j*nx+nx-1 for j in range(1,ny)]+list(range((ny-1)*nx+nx-2,(ny-1)*nx-1,-1))+[j*nx for j in range(ny-2,0,-1)]+[0]
    r.a.curve_line('Rug / continuous edge binding',[(verts[k][0]-.72,verts[k][1]-.40,verts[k][2]+.024) for k in perimeter],bound,.003)
    ob['construction']='5 mm woven body; continuous bound edge; localized gravity-resting corner'
    return ob


def slab_uv(ob,tile,kind='top',edge_x=0,top_z=1):
    """Unfold a slab: waterfall continues the top vein across the miter."""
    uv=ob.data.uv_layers.active
    for p in ob.data.polygons:
        for k in p.loop_indices:
            co=ob.data.vertices[ob.data.loops[k].vertex_index].co+ob.location
            if kind=='waterfall':
                sign=-1 if edge_x<0 else 1
                x=edge_x+sign*(top_z-co.z);y=co.y
            elif kind=='wall':x=co.x;y=co.z
            else:x=co.x;y=co.y
            uv.data[k].uv=(.50+x/tile,.50+y/tile)
    ob['slab_mapping']='Unfolded physical slab coordinates / meters'
