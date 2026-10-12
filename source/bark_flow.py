"""Original opt-in bark relief layered on the retained wood palette.

Spacing follows prototype-local metre coordinates and the authored instance
scale. Height is a world-metre Bump input. This is a visual bark recipe, not a
measured botanical material or scattering law.
"""
import json,math
from pathlib import Path

def load_profile(path,identifier="victorville_dv3_dry_bark"):
 p=json.loads(Path(path).read_text())["profiles"][identifier];assert p["schema_version"]==1
 for band in p["bands"]:
  assert all(math.isfinite(v) and v>0 for v in band.values())
  assert band["height_m"]<band["cross_m"]<band["along_m"]
 return p

def install(material,profile):
 n=material.node_tree.nodes;l=material.node_tree.links;b=n.get('Principled BSDF');assert b is not None
 def node(t,name):a=n.new(t);a.name='DV bark '+name;return a
 def wire(v,s):
  if hasattr(v,'is_output'):l.new(v,s)
  else:s.default_value=v
 def math(op,*v):
  a=node('ShaderNodeMath',op);a.operation=op
  for i,x in enumerate(v):wire(x,a.inputs[i])
  return a.outputs[0]
 attr=node('ShaderNodeAttribute','axis flow in prototype metres');attr.attribute_name='dv_wood_flow_m'
 valid=node('ShaderNodeAttribute','fused wood mask');valid.attribute_name='dv_wood_flow_valid'
 heights=[];variations=[]
 for i,band in enumerate(profile['bands']):
  cross,along,height=band['cross_m'],band['along_m'],band['height_m'];noise_profile=profile['noise']
  mapping=node('ShaderNodeVectorMath','metric anisotropy');mapping.operation='MULTIPLY';l.new(attr.outputs['Vector'],mapping.inputs[0]);mapping.inputs[1].default_value=(1/cross,1/cross,1/along)
  noise=node('ShaderNodeTexNoise','irregular longitudinal weathering');noise.noise_dimensions=noise_profile['dimensions'];noise.inputs['Scale'].default_value=noise_profile['scale'];noise.inputs['Detail'].default_value=noise_profile['detail'];noise.inputs['Roughness'].default_value=noise_profile['roughness'];noise.inputs['W'].default_value=noise_profile['phase_origin']+noise_profile['phase_step']*i;l.new(mapping.outputs[0],noise.inputs['Vector'])
  value=math('SUBTRACT',noise.outputs['Fac'],.5);heights.append(math('MULTIPLY',value,height));variations.append(value)
 h=math('MULTIPLY',math('ADD',*heights),valid.outputs['Fac']);bump=node('ShaderNodeBump','metre fine bark relief');bump.inputs['Strength'].default_value=1.;bump.inputs['Distance'].default_value=1.;l.new(h,bump.inputs['Height'])
 if b.inputs['Normal'].is_linked:l.new(b.inputs['Normal'].links[0].from_socket,bump.inputs['Normal'])
 l.new(bump.outputs['Normal'],b.inputs['Normal'])
 source=b.inputs['Base Color'].links[0].from_socket if b.inputs['Base Color'].is_linked else tuple(b.inputs['Base Color'].default_value)
 factor=math('ADD',1.,math('MULTIPLY',math('MULTIPLY',variations[0],profile['pigment_scalar_contrast']),valid.outputs['Fac']));mul=node('ShaderNodeMixRGB','retained chromatic endpoints');mul.blend_type='MULTIPLY';mul.inputs[0].default_value=1.;wire(source,mul.inputs[1]);wire(factor,mul.inputs[2]);l.new(mul.outputs[0],b.inputs['Base Color'])
 return {'profile':profile['profile'],'world_height_coefficient_m':[b['height_m'] for b in profile['bands']],'prototype_spacing_cross_along_m':[[b['cross_m'],b['along_m']] for b in profile['bands']],'chromatic_endpoints_retained':True,'material_owner':'cybrdelic/cybr-material'}
