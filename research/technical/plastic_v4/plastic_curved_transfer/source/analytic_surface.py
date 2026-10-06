"""Exact original geometric construction, without a new molding-physics claim."""
import numpy as np
DRAFT=np.tan(np.deg2rad(.8));HALF=np.array([.037,.037,.001])
def base_and_direction(q):
 q=np.asarray(q,dtype='f8');c=np.clip(q,-HALF,HALF);dv=q-c;length=np.linalg.norm(dv,axis=-1,keepdims=True);n=np.divide(dv,length,out=np.zeros_like(dv),where=length>0);s=c+.003*n;s=s+np.array([0,0,.004]);p=s.copy();p[...,:2]*=(1-DRAFT*np.abs(s[...,2]-.004)/.04)[...,None];seam=.000022*np.exp(-.5*((s[...,2]-.004)/.000027)**2)*np.maximum(0,1-n[...,2]**2);p[...,:2]+=n[...,:2]*seam[...,None]
 underside=n[...,2]<-.999
 for cx in [-.023,.023]:
  for cy in [-.023,.023]:
   dist=np.hypot(p[...,0]-cx,p[...,1]-cy);ring=.000006*np.exp(-.5*((dist-.0022)/.000095)**2);core=.000002/(1+np.exp(np.minimum(50,(dist-.0022)/.00006)));p[...,2]+=underside*(ring+core)
 return p,n*np.maximum(n[...,2],0)[...,None]**6

def make_shader_displacement(mat,height_image=None,constant_height=None):
 """Object transform must be identity. All metric expressions match base_and_direction."""
 nodes=mat.node_tree.nodes;links=mat.node_tree.links
 def m(op,*args):
  nd=nodes.new('ShaderNodeMath');nd.operation=op
  for i,a in enumerate(args):
   if isinstance(a,(int,float,np.floating)):nd.inputs[i].default_value=float(a)
   else:links.new(a,nd.inputs[i])
  return nd.outputs[0]
 def v(op,*args):
  nd=nodes.new('ShaderNodeVectorMath');nd.operation=op
  for i,a in enumerate(args):
   ix=3 if op=='SCALE' and i==1 else i
   if isinstance(a,(tuple,list,np.ndarray,int,float,np.floating)):nd.inputs[ix].default_value=a
   else:links.new(a,nd.inputs[ix])
  return nd.outputs[0]
 def xyz(socket):
  nd=nodes.new('ShaderNodeSeparateXYZ');links.new(socket,nd.inputs[0]);return list(nd.outputs)
 def combine(x,y,z):
  nd=nodes.new('ShaderNodeCombineXYZ')
  for i,a in enumerate([x,y,z]):
   if isinstance(a,(int,float)):nd.inputs[i].default_value=a
   else:links.new(a,nd.inputs[i])
  return nd.outputs[0]
 attr=nodes.new('ShaderNodeAttribute');attr.attribute_name='original_rounding_chart_q';q=attr.outputs['Vector'];c=v('MINIMUM',v('MAXIMUM',q,(-.037,-.037,-.001)),(.037,.037,.001));n=v('NORMALIZE',v('SUBTRACT',q,c));s=v('ADD',v('ADD',c,v('SCALE',n,.003)),(0,0,.004));sx,sy,sz=xyz(s);nx,ny,nz=xyz(n);draft=m('SUBTRACT',1,m('MULTIPLY',DRAFT/.04,m('ABSOLUTE',m('SUBTRACT',sz,.004))));px=m('MULTIPLY',sx,draft);py=m('MULTIPLY',sy,draft)
 seam=m('MULTIPLY',m('MULTIPLY',.000022,m('EXPONENT',m('MULTIPLY',-.5,m('POWER',m('DIVIDE',m('SUBTRACT',sz,.004),.000027),2)))),m('MAXIMUM',0,m('SUBTRACT',1,m('MULTIPLY',nz,nz))))
 px=m('ADD',px,m('MULTIPLY',nx,seam));py=m('ADD',py,m('MULTIPLY',ny,seam));pz=sz;mask=m('LESS_THAN',nz,-.999)
 for cx in [-.023,.023]:
  for cy in [-.023,.023]:
   dist=m('SQRT',m('ADD',m('POWER',m('SUBTRACT',px,cx),2),m('POWER',m('SUBTRACT',py,cy),2)));ring=m('MULTIPLY',.000006,m('EXPONENT',m('MULTIPLY',-.5,m('POWER',m('DIVIDE',m('SUBTRACT',dist,.0022),.000095),2))));core=m('DIVIDE',.000002,m('ADD',1,m('EXPONENT',m('MINIMUM',50,m('DIVIDE',m('SUBTRACT',dist,.0022),.00006)))));pz=m('ADD',pz,m('MULTIPLY',mask,m('ADD',ring,core)))
 base=combine(px,py,pz);direction=v('SCALE',n,m('POWER',m('MAXIMUM',nz,0),6));geometry=nodes.new('ShaderNodeNewGeometry');geometry.name='Actual evaluation position and normals'
 if constant_height is not None:h=float(constant_height)
 else:
  assert height_image;tex=nodes.new('ShaderNodeTexImage');tex.image=height_image;tex.interpolation='Linear';tex.extension='EXTEND';uv=v('ADD',v('SCALE',base,12.5),(.5,.5,0));links.new(uv,tex.inputs['Vector']);h=m('MULTIPLY',m('SUBTRACT',tex.outputs['Color'],.5),80e-6)
 delta=v('ADD',v('SUBTRACT',base,geometry.outputs['Position']),v('SCALE',direction,h));disp=nodes.new('ShaderNodeVectorDisplacement');disp.space='OBJECT';disp.inputs['Midlevel'].default_value=0;disp.inputs['Scale'].default_value=1;links.new(delta,disp.inputs['Vector']);out=next(n for n in nodes if n.type=='OUTPUT_MATERIAL');links.new(disp.outputs[0],out.inputs['Displacement']);mat.displacement_method='DISPLACEMENT'
 return {'chart':q,'base':base,'direction':direction,'height':h,'geometry':geometry,'math':m,'vector':v}
