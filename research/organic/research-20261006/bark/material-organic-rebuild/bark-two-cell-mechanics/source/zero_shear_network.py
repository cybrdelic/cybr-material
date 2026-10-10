"""Same stock and energy, with an exact zero mean top-minus-bottom x constraint.

The right-top x DOF is dependent, not pinned: x5 = x0+x1+x2-x3-x4.
All five unconstrained junction x coordinates can move through this elimination.
The generalized constraint coordinate is mean(top x)-mean(bottom x).
"""
import numpy as np
from two_cell_network import TwoCellNetwork


class ZeroAverageShearNetwork(TwoCellNetwork):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.dependent=10
        self.free=self.free[self.free!=self.dependent]
        self.coeff=np.zeros(len(self.free))
        for dof,value in [(0,1),(2,1),(4,1),(6,-1),(8,-1)]:
            self.coeff[self.free==dof]=value
        self.shear_gradient=np.zeros(self.ndof)
        self.shear_gradient[[0,2,4]]=-1/3
        self.shear_gradient[[6,8,10]]=1/3

    def control(self,q,strain):
        out=super().control(q,strain)
        out[10]=out[0]+out[2]+out[4]-out[6]-out[8]
        return out

    def reduced_gradient(self,g):
        return g[self.free]+self.coeff*g[self.dependent]

    def reduced_hessian(self,H):
        f=self.free;c=self.coeff;d=self.dependent
        return (H[np.ix_(f,f)]+np.outer(H[f,d],c)+np.outer(c,H[d,f])
                +H[d,d]*np.outer(c,c))

    def increment(self,q,dq,scale=1.):
        out=q.copy();out[self.free]+=scale*dq
        out[self.dependent]+=scale*(self.coeff@dq)
        return out

    def support_gradient(self,g):
        # g = C^T lambda + prescribed-coordinate reactions at equilibrium.
        multiplier=3*g[self.dependent]
        support=multiplier*self.shear_gradient
        support[self.fixed]=g[self.fixed]
        return support

    def solve(self,*args,**kwargs):
        q,receipt=super().solve(*args,**kwargs)
        _,g,_=self.evaluate(q)
        receipt['mean_shear_displacement_m']=float(self.shear_gradient@q*self.L)
        receipt['generalized_shear_reaction_N']=float(3*g[self.dependent]*self.unit/self.L)
        receipt['support_force_and_moment_include_shear_constraint']=True
        receipt['passed']=bool(receipt['passed'] and abs(self.shear_gradient@q)<1e-12)
        return q,receipt

    def topology(self):
        result=super().topology()
        result['boundary_conditions']='Six prescribed radial junction coordinates; one x translation anchor; exact zero mean(top x)-mean(bottom x); all other compatible junction motions, interior positions, and six junction rotations free'
        result['shear_constraint']='mean(top junction x)-mean(bottom junction x)=reference value=0, imposed by exact linear elimination'
        result['dependent_dof']=10
        return result
