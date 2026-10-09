"""Exact retained expressions and native bilinear slopes, shared by CPU/node backends.
No finite-difference bump, filtered derivative, or fitted material coefficient.
"""
from pathlib import Path
import numpy as np

class Dual:
    def __init__(self,op,v,u=0.,w=0.):self.op=op;self.v=v;self.u=u;self.w=w
    def cast(self,x):return x if isinstance(x,Dual) else Dual(self.op,x)
    def __add__(self,x):
        x=self.cast(x);o=self.op;return Dual(o,o('ADD',self.v,x.v),o('ADD',self.u,x.u),o('ADD',self.w,x.w))
    __radd__=__add__
    def __neg__(self):return self*-1
    def __sub__(self,x):return self+-self.cast(x)
    def __rsub__(self,x):return self.cast(x)+-self
    def __mul__(self,x):
        x=self.cast(x);o=self.op;return Dual(o,o('MULTIPLY',self.v,x.v),o('ADD',o('MULTIPLY',self.u,x.v),o('MULTIPLY',self.v,x.u)),o('ADD',o('MULTIPLY',self.w,x.v),o('MULTIPLY',self.v,x.w)))
    __rmul__=__mul__
    def __truediv__(self,x):return self*self.cast(x).power(-1)
    def __rtruediv__(self,x):return self.cast(x)*self.power(-1)
    def power(self,p):
        o=self.op;v=o('POWER',self.v,p);scale=o('MULTIPLY',p,o('POWER',self.v,p-1))
        return Dual(o,v,o('MULTIPLY',scale,self.u),o('MULTIPLY',scale,self.w))
    def exp(self):
        o=self.op;v=o('EXPONENT',self.v);return Dual(o,v,o('MULTIPLY',v,self.u),o('MULTIPLY',v,self.w))
    def absolute(self):
        o=self.op;sign=o('SUBTRACT',o('GREATER_THAN',self.v,0),o('LESS_THAN',self.v,0));return Dual(o,o('ABSOLUTE',self.v),o('MULTIPLY',sign,self.u),o('MULTIPLY',sign,self.w))
    def clip(self,lo,hi):
        o=self.op;inside=o('MULTIPLY',o('GREATER_THAN',self.v,lo),o('LESS_THAN',self.v,hi));return Dual(o,o('MINIMUM',o('MAXIMUM',self.v,lo),hi),o('MULTIPLY',inside,self.u),o('MULTIPLY',inside,self.w))
    def minimum(self,hi):
        o=self.op;inside=o('LESS_THAN',self.v,hi);return Dual(o,o('MINIMUM',self.v,hi),o('MULTIPLY',inside,self.u),o('MULTIPLY',inside,self.w))

def field_expressions(o,q,tu,tv,native_cell_override=None):
    """q coordinates and derivatives are in the object's retained meter chart."""
    Q=[Dual(o,q[i],tu[i],tv[i]) for i in range(3)]
    half=[.037,.037,.001];c=[Q[i].clip(-half[i],half[i]) for i in range(3)];w=[Q[i]-c[i] for i in range(3)];r=sum(v*v for v in w).power(.5);n=[v/r for v in w]
    s=[c[i]+.003*n[i]+(0,0,.004)[i] for i in range(3)];z=s[2]-.004;draft=1-np.tan(np.deg2rad(.8))/.04*z.absolute();seam=22e-6*(-.5*(z/27e-6).power(2)).exp()*(1-n[2]*n[2]).clip(0,1)
    b=[s[0]*draft+n[0]*seam,s[1]*draft+n[1]*seam,s[2]]
    underside=o('LESS_THAN',n[2].v,-.999)
    for cx in [-.023,.023]:
        for cy in [-.023,.023]:
            dist=((b[0]-cx).power(2)+(b[1]-cy).power(2)).power(.5)
            ring=6e-6*(-.5*((dist-.0022)/95e-6).power(2)).exp();core=2e-6/(1+((dist-.0022)/60e-6).minimum(50).exp());b[2]=b[2]+(ring+core)*underside
    d=[v*n[2].clip(0,1).power(6) for v in n]
    # Pixel-center EXTEND, bottom-up native map. Four closest samples own one cell.
    raw=[o('SUBTRACT',o('MULTIPLY',o('ADD',o('DIVIDE',b[i].v,.08),.5),4096),.5) for i in range(2)]
    tex=[o('MINIMUM',o('MAXIMUM',r,0),4095) for r in raw];idx=[o('FLOOR',x) for x in tex]
    if native_cell_override is not None:idx=list(np.asarray(native_cell_override).T)
    f=[o('SUBTRACT',tex[i],idx[i]) for i in range(2)]
    x0,y0=idx;x1=o('MINIMUM',o('ADD',x0,1),4095);y1=o('MINIMUM',o('ADD',y0,1),4095)
    h00,h10,h01,h11=[o.sample_height(x,y) for x,y in [(x0,y0),(x1,y0),(x0,y1),(x1,y1)]]
    x,y=f;omx=o('SUBTRACT',1,x);omy=o('SUBTRACT',1,y)
    add=lambda a,b:o('ADD',a,b);mul=lambda a,b:o('MULTIPLY',a,b);sub=lambda a,b:o('SUBTRACT',a,b)
    H=add(add(mul(mul(omx,omy),h00),mul(mul(x,omy),h10)),add(mul(mul(omx,y),h01),mul(mul(x,y),h11)))
    hx=o('DIVIDE',add(mul(omy,sub(h10,h00)),mul(y,sub(h11,h01))),.08/4096)
    hy=o('DIVIDE',add(mul(omx,sub(h01,h00)),mul(x,sub(h11,h10))),.08/4096)
    hx=mul(hx,mul(o('SUBTRACT',1,o('LESS_THAN',raw[0],0)),o('LESS_THAN',raw[0],4095)));hy=mul(hy,mul(o('SUBTRACT',1,o('LESS_THAN',raw[1],0)),o('LESS_THAN',raw[1],4095)))
    Hd=Dual(o,H,add(mul(hx,b[0].u),mul(hy,b[1].u)),add(mul(hx,b[0].w),mul(hy,b[1].w)))
    F=[b[i]+d[i]*Hd for i in range(3)];u=[v.u for v in F];v=[v.w for v in F];cross=[sub(mul(u[1],v[2]),mul(u[2],v[1])),sub(mul(u[2],v[0]),mul(u[0],v[2])),sub(mul(u[0],v[1]),mul(u[1],v[0]))];length=o('SQRT',add(add(mul(cross[0],cross[0]),mul(cross[1],cross[1])),mul(cross[2],cross[2])));normal=[o('DIVIDE',c,length) for c in cross]
    return {'position':[v.v for v in F],'normal':normal,'cell':idx,'height':H,'height_gradient':[hx,hy],'base':[v.v for v in b]}

class NumericOps:
    def __init__(self,height,dtype='f8'):self.height=height;self.dtype=np.dtype(dtype)
    def __call__(self,op,*args):
        a=[np.asarray(x,dtype=self.dtype) for x in args]
        with np.errstate(divide='ignore',invalid='ignore',over='ignore'):
            if op=='ADD':r=a[0]+a[1]
            elif op=='SUBTRACT':r=a[0]-a[1]
            elif op=='MULTIPLY':r=a[0]*a[1]
            elif op=='DIVIDE':r=np.divide(a[0],a[1],out=np.zeros(np.broadcast_shapes(a[0].shape,a[1].shape),self.dtype),where=a[1]!=0)
            elif op=='POWER':r=np.where(a[0]==0,0,np.power(a[0],a[1]))
            elif op=='EXPONENT':r=np.exp(a[0])
            elif op=='ABSOLUTE':r=np.abs(a[0])
            elif op=='MINIMUM':r=np.minimum(a[0],a[1])
            elif op=='MAXIMUM':r=np.maximum(a[0],a[1])
            elif op=='GREATER_THAN':r=a[0]>a[1]
            elif op=='LESS_THAN':r=a[0]<a[1]
            elif op=='FLOOR':r=np.floor(a[0])
            elif op=='SQRT':r=np.sqrt(a[0])
            else:raise ValueError(op)
        return np.asarray(r,dtype=self.dtype)
    def sample_height(self,x,y):
        # Match stored uint16 normalized samples and float32 shader decode order.
        encoded=np.asarray(self.height[y.astype('i4'),x.astype('i4')]/80e-6+.5,dtype=self.dtype)
        return self('MULTIPLY',self('SUBTRACT',encoded,.5),80e-6)

class BlenderOps:
    def __init__(self,nodes,links,image):self.nodes=nodes;self.links=links;self.image=image
    def __call__(self,op,*args):
        if all(isinstance(a,(int,float,np.number)) for a in args):
            return float(NumericOps(None)(op,*args))
        n=self.nodes.new('ShaderNodeMath');n.operation=op
        for i,a in enumerate(args):
            if isinstance(a,(int,float,np.number)):n.inputs[i].default_value=float(a)
            else:self.links.new(a,n.inputs[i])
        return n.outputs[0]
    def combine(self,xyz):
        n=self.nodes.new('ShaderNodeCombineXYZ')
        for i,a in enumerate(xyz):
            if isinstance(a,(int,float,np.number)):n.inputs[i].default_value=float(a)
            else:self.links.new(a,n.inputs[i])
        return n.outputs[0]
    def separate(self,v):
        n=self.nodes.new('ShaderNodeSeparateXYZ');self.links.new(v,n.inputs[0]);return list(n.outputs)
    def sample_height(self,x,y):
        xyz=[self('DIVIDE',self('ADD',p,.5),4096) for p in [x,y]]+[0.]
        n=self.nodes.new('ShaderNodeTexImage');n.image=self.image;n.interpolation='Closest';n.extension='EXTEND';self.links.new(self.combine(xyz),n.inputs['Vector']);sep=self.nodes.new('ShaderNodeSeparateColor');sep.mode='RGB';self.links.new(n.outputs['Color'],sep.inputs[0]);return self('MULTIPLY',self('SUBTRACT',sep.outputs[0],.5),80e-6)


def build_group(image):
    import bpy
    g=bpy.data.node_groups.new('Retained physical plastic / analytic field','ShaderNodeTree')
    for name in ['Original q','Tangent U','Tangent V']:g.interface.new_socket(name=name,in_out='INPUT',socket_type='NodeSocketVector')
    for name in ['Normal object','Position object','Native cell','Base object']:g.interface.new_socket(name=name,in_out='OUTPUT',socket_type='NodeSocketVector')
    g.interface.new_socket(name='Native height',in_out='OUTPUT',socket_type='NodeSocketFloat')
    n=g.nodes;l=g.links;i=n.new('NodeGroupInput');out=n.new('NodeGroupOutput');o=BlenderOps(n,l,image);r=field_expressions(o,o.separate(i.outputs['Original q']),o.separate(i.outputs['Tangent U']),o.separate(i.outputs['Tangent V']))
    for name,value in [('Normal object',r['normal']),('Position object',r['position']),('Native cell',r['cell']+[0]),('Base object',r['base'])]:l.new(o.combine(value),out.inputs[name])
    l.new(r['height'],out.inputs['Native height'])
    return g
