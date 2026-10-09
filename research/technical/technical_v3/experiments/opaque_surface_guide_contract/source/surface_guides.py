"""Declared OIDN first-hit features for opaque single-Principled materials.

This does not edit the BSDF graph or radiance. It adds shader AOV outputs.
OIDN permits normal-incidence reflectivity for metals and diffuse color for
matte surfaces: https://www.openimagedenoise.org/documentation.html#albedos
Unsupported transmission, emission and multi-closure materials fail closed.
The guide is approximate surface reflectivity, not a solved optical albedo.
"""
import bpy

ALBEDO='OIDN_SurfaceReflectivity'
NORMAL='OIDN_SurfaceNormal'

def configure(scene):
    mats={m for o in scene.objects if not o.hide_render for m in getattr(o.data,'materials',[]) if m}
    records=[]
    for m in mats:
        if not m.use_nodes:raise ValueError('Node material required: '+m.name)
        nd=m.node_tree.nodes;lk=m.node_tree.links
        outputs=[n for n in nd if n.type=='OUTPUT_MATERIAL' and n.is_active_output]
        if len(outputs)!=1:raise ValueError('One active material output required')
        surf=outputs[0].inputs['Surface']
        if not surf.is_linked or surf.links[0].from_node.type!='BSDF_PRINCIPLED':
            raise ValueError('Only direct Principled output supported: '+m.name)
        p=surf.links[0].from_node
        alpha=p.inputs['Alpha']
        if alpha.is_linked or alpha.default_value != 1:
            raise ValueError('Opaque Alpha=1 required; linked Alpha is unsupported: '+m.name)
        for name in ['Transmission Weight','Subsurface Weight']:
            a=p.inputs[name]
            if a.is_linked or a.default_value!=0:raise ValueError('Unsupported '+name+': '+m.name)
        em=p.inputs['Emission Color'];strength=p.inputs['Emission Strength']
        if em.is_linked or strength.is_linked or strength.default_value*max(em.default_value[:3])!=0:
            raise ValueError('Emissive source requires a separate feature contract')
        if outputs[0].inputs['Volume'].is_linked:raise ValueError('Volume is unsupported')
        for name in [ALBEDO,NORMAL]:
            if any(n.type=='OUTPUT_AOV' and n.name==name for n in nd):raise ValueError('AOV already configured')
        alb=nd.new('ShaderNodeOutputAOV');alb.name=ALBEDO
        src=p.inputs['Base Color']
        if src.is_linked:lk.new(src.links[0].from_socket,alb.inputs['Color'])
        else:
            if min(src.default_value[:3])<0 or max(src.default_value[:3])>1:raise ValueError('Unbounded reflectivity')
            alb.inputs['Color'].default_value=src.default_value
        normal=nd.new('ShaderNodeOutputAOV');normal.name=NORMAL
        if p.inputs['Normal'].is_linked:
            normal_source=p.inputs['Normal'].links[0].from_socket
            if normal_source.node.type not in {'NORMAL_MAP','BUMP'} or normal_source.name!='Normal':
                raise ValueError('Unsupported linked shading-normal semantics: '+m.name)
            lk.new(normal_source,normal.inputs['Color'])
        else:
            geom=nd.new('ShaderNodeNewGeometry');geom.name='OIDN feature normal only'
            lk.new(geom.outputs['Normal'],normal.inputs['Color'])
        records.append({'material':m.name,'reflectivity':'linked Base Color or identical constant','normal':'same linked shading normal or Geometry Normal','BSDF_output_unchanged':True})
    for vl in scene.view_layers:
        for name in [ALBEDO,NORMAL]:
            if any(a.name==name for a in vl.aovs):raise ValueError('View layer AOV already configured')
            a=vl.aovs.add();a.name=name;a.type='COLOR'
    return records
