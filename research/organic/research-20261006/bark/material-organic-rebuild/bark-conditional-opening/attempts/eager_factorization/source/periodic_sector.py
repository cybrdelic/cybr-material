"""Rotation-aware periodic condensation of the fixed-reference prism coupon."""
import numpy as np
from scipy.sparse import coo_matrix,kron,eye
from scipy.sparse.linalg import splu
from volume_coupon import Coupon,configuration,ROOT


def rotation_y(angle):
    c=np.cos(angle);s=np.sin(angle)
    return np.array([[c,0.,s],[0.,1.,0.],[-s,0.,c]])


class PeriodicSector:
    def __init__(self,pitch=.006,sector_multiplier=1,radial_pitch=.001):
        config=configuration();config['width_m']=config['outer_radius_m']*2*np.pi/14*sector_multiplier
        config['height_m']=.012;config['radial_pitch_m']=radial_pitch
        self.base=Coupon(pitch,precondition=True,config_override=config)
        self.config=config;self.pitch=pitch;self.multiplier=sector_multiplier
        self.alpha=config['width_m']/config['outer_radius_m']
        b=self.base;X=b.X;r=np.linalg.norm(X[:,[0,2]],axis=1)
        sx=np.arctan2(X[:,0],X[:,2])*config['outer_radius_m'];y=X[:,1]
        right=np.isclose(sx,config['width_m']/2,atol=1e-11,rtol=0)
        upper=np.isclose(y,config['height_m']/2,atol=1e-11,rtol=0)
        keys=np.column_stack((r,np.where(right,-config['width_m']/2,sx),np.where(upper,-config['height_m']/2,y)))
        canonical_values=keys.copy();keys=np.round(keys,12)
        _,first,groups=np.unique(keys,axis=0,return_index=True,return_inverse=True)
        self.groups=groups;self.X=X[first].copy();self.representatives=first
        self.rotations=np.broadcast_to(np.eye(3),(len(X),3,3)).copy()
        self.rotations[right]=rotation_y(self.alpha)
        # Representatives selected by np.unique may be on the right/upper edge.
        # Define coordinates in the canonical left/lower chart explicitly.
        theta=canonical_values[first,1]/config['outer_radius_m']
        self.X=np.column_stack((canonical_values[first,0]*np.sin(theta),canonical_values[first,2],canonical_values[first,0]*np.cos(theta)))
        rows=np.broadcast_to(np.arange(len(X))[:,None,None]*3+np.arange(3)[None,:,None],self.rotations.shape)
        cols=np.broadcast_to(groups[:,None,None]*3+np.arange(3)[None,None,:],self.rotations.shape)
        self.T=coo_matrix((self.rotations.ravel(),(rows.ravel(),cols.ravel())),shape=(len(X)*3,len(self.X)*3)).tocsr()
        self.fixed_nodes=np.unique(groups[b.fixed_nodes]);self.free_nodes=np.setdiff1d(np.arange(len(self.X)),self.fixed_nodes)
        self.free=(self.free_nodes[:,None]*3+np.arange(3)).ravel()
        self.K=(self.T.T@b.K@self.T).tocsc()
        self.factor=splu(self.K[self.free,:][:,self.free])
        self.mass_matrix=(self.T.T@kron(b.mass_matrix,eye(3),format='csr')@self.T).tocsr()
        self.mass=b.mass;self.volume=b.volume;self.right_nodes=np.flatnonzero(right);self.upper_nodes=np.flatnonzero(upper)
        self.reference_reconstruction_error=float(np.max(np.abs(self.expand(self.X)+np.column_stack((np.zeros(len(X)),np.where(upper,config['height_m'],0.),np.zeros(len(X))))-X)))
        assert self.reference_reconstruction_error<2e-12

    @property
    def load(self):return self.base.load
    def set_load(self,s):self.base.set_load(s)
    def expand(self,v):return np.asarray(self.T@np.asarray(v).ravel()).reshape(self.base.X.shape)
    def evaluate(self,v,*args,**kwargs):
        E,g,r=self.base.evaluate(self.expand(v),*args,**kwargs)
        condensed=np.asarray(self.T.T@g.ravel()).reshape(self.X.shape)
        r['maximum_free_force_N']=float(np.max(np.abs(condensed[self.free_nodes])))
        r['sector_angle_radians']=self.alpha
        return E,condensed,r
    def equilibrate(self,v):return Coupon.equilibrate(self,v)
    def kinetic(self,v):return .5*float(v.ravel()@(self.mass_matrix@v.ravel()))


if __name__=='__main__':
    m=PeriodicSector();print(dict(full_nodes=len(m.base.X),periodic_nodes=len(m.X),reference_reconstruction_error=m.reference_reconstruction_error))
