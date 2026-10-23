"""Editable three-dimensional timber finish, in object-space meters.

The axis is local Y. Each object's stable random value selects a different
saw plane through the same growth volume. Face and end grain therefore meet
at the edge rather than re-projecting a longitudinal image on every face.
Native finite 4K anatomical maps remain available for other DCCs.
"""
def install(group,spec):
    n=group.nodes;l=group.links
    p=next(v for v in n if v.type=='BSDF_PRINCIPLED')
    gi=next(v for v in n if v.type=='GROUP_INPUT')
    go=next(v for v in n if v.type=='GROUP_OUTPUT')
    oak=spec['id'].startswith('04')
    def node(t,label):
        v=n.new(t);v.label=label;return v
    def math(op,a,b=0):
        v=node('ShaderNodeMath','Timber / '+op);v.operation=op
        for i,x in enumerate((a,b)):
            if isinstance(x,(int,float)):v.inputs[i].default_value=x
            else:l.new(x,v.inputs[i])
        return v.outputs[0]
    def ramp(source,stops,label):
        v=node('ShaderNodeValToRGB',label);r=v.color_ramp
        while len(r.elements)>2:r.elements.remove(r.elements[-1])
        for i,(position,col) in enumerate(stops):
            e=r.elements[i] if i<2 else r.elements.new(position);e.position=position;e.color=(*col,1)
        l.new(source,v.inputs[0]);return v.outputs[0]
    def noise(vector,scale,detail=2):
        v=node('ShaderNodeTexNoise','Timber / finite cell bundles');v.noise_dimensions='3D'
        v.inputs['Scale'].default_value=scale;v.inputs['Detail'].default_value=detail
        v.inputs['Roughness'].default_value=.65;l.new(vector,v.inputs['Vector']);return v.outputs['Fac']
    def vector_scale(vector,scale):
        v=node('ShaderNodeVectorMath','Timber / physical cell proportions');v.operation='MULTIPLY'
        l.new(vector,v.inputs[0]);v.inputs[1].default_value=scale;return v.outputs['Vector']
    geom=node('ShaderNodeNewGeometry','Timber / meter-space position')
    trans=node('ShaderNodeVectorTransform','Timber / coherent local volume');trans.vector_type='POINT';trans.convert_from='WORLD';trans.convert_to='OBJECT'
    l.new(geom.outputs['Position'],trans.inputs[0])
    mapping=node('ShaderNodeMapping','Timber / fabrication grain orientation')
    l.new(trans.outputs[0],mapping.inputs['Vector']);l.new(gi.outputs['Timber Axis Rotation'],mapping.inputs['Rotation'])
    pos=mapping.outputs[0]
    info=node('ShaderNodeObjectInfo','Timber / individual saw cut')
    sep=node('ShaderNodeSeparateXYZ','Timber / Y is fiber axis');l.new(pos,sep.inputs[0])
    # A different pith offset and cut depth per board, never arbitrary UV tint.
    rx=math('ADD',sep.outputs['X'],math('ADD',math('MULTIPLY',info.outputs['Random'],.20),-.03))
    rz=math('ADD',sep.outputs['Z'],math('ADD',math('MULTIPLY',info.outputs['Random'],.040),.024))
    rz=math('ADD',rz,math('ADD',math('MULTIPLY',sep.outputs['Y'],.025),math('MULTIPLY',math('MULTIPLY',sep.outputs['Y'],sep.outputs['Y']),.10)))
    radial=math('SQRT',math('ADD',math('MULTIPLY',rx,rx),math('MULTIPLY',rz,rz)))
    broad=noise(vector_scale(pos,(1,.055,1)),75,3)
    wobble=math('MULTIPLY',math('SUBTRACT',broad,.5),.00075)
    rings=math('MULTIPLY',math('ADD',radial,wobble),1250 if oak else 1570)
    phase=math('MULTIPLY',math('ADD',math('SINE',rings),1),.5)
    # Cell bundles follow radial growth coordinates, so fibers bend through
    # a flat-sawn cathedral instead of crossing growth as straight streaks.
    theta=math('ARCTAN2',rx,rz)
    cells=node('ShaderNodeCombineXYZ','Timber / radial-longitudinal anatomy')
    l.new(math('ADD',radial,wobble),cells.inputs['X']);l.new(sep.outputs['Y'],cells.inputs['Y']);l.new(math('MULTIPLY',theta,.018),cells.inputs['Z'])
    bundles=noise(vector_scale(cells.outputs[0],(1,.022,1)),760,3)
    fibers=noise(vector_scale(cells.outputs[0],(1,.009,1)),4700,2)
    early=math('POWER',phase,9 if oak else 3)
    pore=math('MULTIPLY',math('POWER',fibers,7),math('ADD',math('MULTIPLY',early,.8),.2))
    pigment=math('ADD',math('MULTIPLY',bundles,.48),math('ADD',math('MULTIPLY',broad,.28),math('MULTIPLY',phase,.15)))
    pigment=math('SUBTRACT',pigment,math('MULTIPLY',pore,.7 if oak else .3))
    colors=[(.060,.038,.018),(.13,.081,.039),(.225,.148,.076)] if oak else [(.027,.010,.004),(.065,.025,.010),(.13,.059,.028)]
    base=ramp(pigment,[(.18,colors[0]),(.48,colors[1]),(.82,colors[2])],'Timber / intrinsic heartwood')
    if oak:
        rays=noise(vector_scale(cells.outputs[0],(160,1600,80)),1,1)
        mask=math('MULTIPLY',math('MAXIMUM',math('SUBTRACT',rays,.63),0),1.4)
        mix=node('ShaderNodeMixRGB','Oak / radial medullary plates');l.new(mask,mix.inputs[0]);l.new(base,mix.inputs[1]);mix.inputs[2].default_value=(.30,.225,.125,1);base=mix.outputs[0]
    tint=node('ShaderNodeMixRGB','Timber / selected board finish');tint.blend_type='MULTIPLY';tint.inputs[0].default_value=1
    l.new(base,tint.inputs[1]);l.new(gi.outputs['Tint'],tint.inputs[2]);l.new(tint.outputs[0],p.inputs['Base Color']);l.new(tint.outputs[0],go.inputs['Base Color'])
    rough=math('ADD',.39 if not oak else .46,math('ADD',math('MULTIPLY',pore,.13),math('MULTIPLY',bundles,.034)))
    l.new(math('ADD',rough,gi.outputs['Roughness Offset']),p.inputs['Roughness'])
    height=math('SUBTRACT',math('MULTIPLY',fibers,.08),math('MULTIPLY',pore,.35))
    bump=node('ShaderNodeBump','Timber / lumen-scale finish relief');bump.inputs['Distance'].default_value=.000085;bump.inputs['Strength'].default_value=.6
    l.new(height,bump.inputs['Height']);l.new(bump.outputs['Normal'],p.inputs['Normal']);l.new(bump.outputs['Normal'],p.inputs['Coat Normal'])
    p.inputs['Coat Weight'].default_value=.08;p.inputs['Coat Roughness'].default_value=.32
    group['timber_axis']='Local Y / meters';group['end_grain']='Same radial growth volume as longitudinal faces'
