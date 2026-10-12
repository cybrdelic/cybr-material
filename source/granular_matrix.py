"""Reusable CYBR granular matrix fields in world-space metres.

Opt-in material building blocks only. Existing selected recipes and catalogue
entries are untouched. Scene deposition, plant placement and finite stock do
not belong here. The profile's unresolved roughness term is a documented
heuristic, not a claim of measured BRDF or exact normal-variance integration.
"""
import json,math
from pathlib import Path

def load_profile(path,identifier='victorville_du5_optical_matrix'):
 profile=json.loads(Path(path).read_text())['profiles'][identifier]
 assert profile['schema_version']==1 and profile['units']=='metres'
 assert all(0<b['wavelength_m']<=.1 and 0<b['height_coefficient_m']<.01 for b in profile['bands'])
 assert all(math.isfinite(v)for b in profile['bands']for v in b.values())
 assert profile['filter']['full_response_cell_pixels']>profile['filter']['zero_response_cell_pixels']>0
 return profile

class Nodes:
 def __init__(self,tree,prefix='DU5 '):self.n=tree.nodes;self.l=tree.links;self.prefix=prefix
 def node(self,kind,name):
  p=self.n.new(kind);p.name=self.prefix+name;return p
 def wire(self,value,socket):
  if hasattr(value,'is_output'):self.l.new(value,socket)
  else:socket.default_value=value
 def calc(self,op,*values):
  p=self.node('ShaderNodeMath',op);p.operation=op
  for i,v in enumerate(values):self.wire(v,p.inputs[i])
  return p.outputs[0]
 def remap(self,value,a,b,c,d):
  p=self.node('ShaderNodeMapRange','smooth metric reconstruction');p.interpolation_type='SMOOTHERSTEP';p.clamp=True
  for key,v in [('Value',value),('From Min',a),('From Max',b),('To Min',c),('To Max',d)]:self.wire(v,p.inputs[key])
  return p.outputs['Result']

def surface_footprint(builder,normal,incoming,view_distance,angular_pitch):
 """Conservative scalar grazing footprint; metric inputs, approximate filter."""
 n,l=builder.n,builder.l;node,calc=builder.node,builder.calc
 dot=node('ShaderNodeVectorMath','grazing surface footprint');dot.operation='DOT_PRODUCT';l.new(normal,dot.inputs[0]);l.new(incoming,dot.inputs[1]);cosine=calc('MAXIMUM',.12,calc('ABSOLUTE',dot.outputs['Value']))
 return calc('DIVIDE',calc('MULTIPLY',view_distance,angular_pitch),cosine)

def aggregate_fields(builder,position_m,footprint_m,relief_gain,profile):
 """Return coupled height [m], achromatic pigment delta and roughness proxy."""
 n,l=builder.n,builder.l;node,calc,remap=builder.node,builder.calc,builder.remap
 warp=profile['warp'];pos=position_m;footprint=footprint_m
 noise=node('ShaderNodeTexNoise','sub-centimetre aggregate warp');noise.inputs['Scale'].default_value=1/warp['wavelength_m'];noise.inputs['Detail'].default_value=warp['detail'];noise.inputs['Roughness'].default_value=warp['roughness'];l.new(pos,noise.inputs['Vector'])
 vsub=node('ShaderNodeVectorMath','zero mean warp');vsub.operation='SUBTRACT';l.new(noise.outputs['Color'],vsub.inputs[0]);vsub.inputs[1].default_value=(.5,.5,.5)
 vs=node('ShaderNodeVectorMath','metric warp 0.8mm');vs.operation='SCALE';l.new(vsub.outputs['Vector'],vs.inputs[0]);vs.inputs['Scale'].default_value=warp['amplitude_m']
 va=node('ShaderNodeVectorMath','irregular aggregate coordinates');va.operation='ADD';l.new(pos,va.inputs[0]);l.new(vs.outputs['Vector'],va.inputs[1]);coord=va.outputs['Vector']
 pigment=0.;height=0.;unresolved=0.
 for i,band in enumerate(profile['bands']):
  size=band['wavelength_m'];amplitude=band['height_coefficient_m'];pigment_contrast=band['scalar_pigment_contrast'];noise_profile=profile['noise']
  grain=node('ShaderNodeTexNoise',f'{size*1000:g}mm irregular aggregate band');grain.noise_dimensions=noise_profile['dimensions'];grain.inputs['Scale'].default_value=1/size;grain.inputs['Detail'].default_value=noise_profile['detail'];grain.inputs['Roughness'].default_value=noise_profile['roughness'];grain.inputs['W'].default_value=noise_profile['phase_origin']+i*noise_profile['phase_step'];l.new(coord,grain.inputs['Vector'])
  value=grain.outputs['Fac'];filterweight=remap(footprint,size/profile['filter']['full_response_cell_pixels'],size/profile['filter']['zero_response_cell_pixels'],1.,0.)
  h=calc('MULTIPLY',calc('SUBTRACT',value,.5),calc('MULTIPLY',filterweight,amplitude));height=calc('ADD',height,h)
  mono=node('ShaderNodeRGBToBW','coherent mineral and aggregate field');l.new(grain.outputs['Color'],mono.inputs['Color'])
  bodytone=calc('ADD',calc('MULTIPLY',calc('SUBTRACT',mono.outputs[0],.5),profile['pigment']['mineral_weight']),calc('MULTIPLY',calc('SUBTRACT',value,.5),profile['pigment']['height_weight']))
  pigment=calc('ADD',pigment,calc('MULTIPLY',bodytone,calc('MULTIPLY',filterweight,pigment_contrast)))
  unresolved=calc('ADD',unresolved,calc('MULTIPLY',calc('SUBTRACT',1.,filterweight),band['unresolved_roughness_proxy']))
 return {'height_m':calc('MULTIPLY',height,relief_gain),'pigment_delta':pigment,'roughness_proxy':unresolved}

def attach_closure(builder,principled,height_m,reflectance_factor,roughness_offset,profile):
 """Layer fields over an existing approved material without replacing hue."""
 n,l=builder.n,builder.l;node,wire,calc=builder.node,builder.wire,builder.calc;b=principled
 source=b.inputs['Base Color'].links[0].from_socket if b.inputs['Base Color'].is_linked else tuple(b.inputs['Base Color'].default_value)
 mul=node('ShaderNodeMixRGB','coupled exposed grain pigment');mul.blend_type='MULTIPLY';mul.inputs[0].default_value=1.;wire(source,mul.inputs[1]);wire(reflectance_factor,mul.inputs[2]);l.new(mul.outputs[0],b.inputs['Base Color'])
 bb=node('ShaderNodeBump','metre aggregate relief over source micro-normal');bb.inputs['Strength'].default_value=1.;bb.inputs['Distance'].default_value=1.;l.new(height_m,bb.inputs['Height'])
 if b.inputs['Normal'].is_linked:l.new(b.inputs['Normal'].links[0].from_socket,bb.inputs['Normal'])
 l.new(bb.outputs['Normal'],b.inputs['Normal'])
 r=b.inputs['Roughness'].links[0].from_socket if b.inputs['Roughness'].is_linked else b.inputs['Roughness'].default_value
 l.new(calc('MINIMUM',profile['roughness_limit'],calc('ADD',r,roughness_offset)),b.inputs['Roughness'])
