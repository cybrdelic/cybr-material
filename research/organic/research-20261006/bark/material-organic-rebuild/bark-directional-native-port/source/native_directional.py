"""Opt-in native batch API for the unchanged, uncalibrated directional envelope.

H=F-I is an explicit primary field. Prism inputs are local displacement increments
relative to a supplied origin gradient H0; never pass world positions as du.
All geometry uses frozen reference gradients/weights, with axes as columns.
"""
from pathlib import Path
import ctypes
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'bark-constitutive-upgrade/source'))
from directional_candidate import synthetic_candidate

DOUBLE=ctypes.POINTER(ctypes.c_double)
INT64=ctypes.POINTER(ctypes.c_int64)
ERRORS={-1:'Empty or invalid batch dimensions',-2:'Nonfinite incremental input',
        -3:'Invalid branch parameters or stretch bounds',-4:'Material frame must be finite, orthonormal and right handed',
        -5:'Inverted, singular or nonfinite deformation',-6:'Degenerate material stretch',
        -7:'Outside declared experimental stretch domain',-8:'Degenerate normalized material frame',
        -9:'Nonfinite constitutive output',-10:'Invalid reference gradient or weight'}

def _ptr(x):return x.ctypes.data_as(DOUBLE)
def _field(x,shape,name,positive=False):
    try:out=np.array(np.broadcast_to(np.asarray(x,dtype=np.float64),shape),dtype=np.float64,order='C',copy=True)
    except (TypeError,ValueError) as exc:raise ValueError(f'{name} cannot broadcast to {shape}') from exc
    if not np.isfinite(out).all() or (positive and np.any(out<=0)):raise ValueError(f'Invalid {name}')
    out.setflags(write=False)
    return out

def _diagnostics(raw):
    return dict(J=raw[...,0],stretches=raw[...,1:4],axial_energy_J_m3=raw[...,4],angular_energy_J_m3=raw[...,5])

class NativeEnvelope:
    """One law, many explicit gradients/axes; no history or constitutive tangent."""
    def __init__(self,law=None):
        law=synthetic_candidate() if law is None else law
        self.axes=_field(law.axes,(3,3),'material axes')
        self.parameters=np.array([[b.E,b.H,b.y,b.D,b.d] for b in (*law.tension,*law.compression)],dtype=np.float64).ravel()
        self.parameters=np.r_[self.parameters,law.G,law.minimum_stretch,law.maximum_stretch]
        if self.parameters.shape!=(33,):raise ValueError('Expected exactly three tension and three compression branches')
        self.parameters.setflags(write=False)
        path=Path(__file__).with_name('libdirectional_native.so')
        if not path.is_file():raise RuntimeError('Run source/build_native.py before using the isolated native port')
        self.library=ctypes.CDLL(str(path))
        self.point_function=self.library.directional_points
        self.point_function.argtypes=[ctypes.c_int64]+[DOUBLE]*6+[INT64]
        self.point_function.restype=ctypes.c_int
        self.prism_function=self.library.directional_prisms
        self.prism_function.argtypes=[ctypes.c_int64]*2+[DOUBLE]*10+[INT64]
        self.prism_function.restype=ctypes.c_int
        # Eagerly validate the frozen parameters and default frame natively.
        self.evaluate_incremental(np.zeros((3,3)))
    def _check(self,rc,bad):
        if rc:raise ValueError(f'{ERRORS.get(rc,"Native evaluation failed")}; status {rc}, flat point {int(bad[0])}')
    def evaluate_incremental(self,H,axes=None):
        h=np.ascontiguousarray(H,dtype=np.float64)
        if h.ndim<2 or h.shape[-2:]!=(3,3) or not np.isfinite(h).all():raise ValueError('Expected finite (...,3,3) incremental gradients')
        count=h.size//9
        if not count:raise ValueError('Empty material-point batch')
        a=_field(self.axes if axes is None else axes,h.shape,'material axes')
        shape=h.shape[:-2];w=np.empty(shape);p=np.empty_like(h);d=np.empty((*shape,6));bad=np.array([-1],dtype=np.int64)
        rc=self.point_function(count,_ptr(h),_ptr(a),_ptr(self.parameters),_ptr(w),_ptr(p),_ptr(d),bad.ctypes.data_as(INT64))
        self._check(rc,bad)
        return w,p,_diagnostics(d)
    def evaluate(self,F,axes=None):
        """Compatibility convenience; tiny diagonal perturbations need the H API."""
        return self.evaluate_incremental(np.asarray(F,dtype=np.float64)-np.eye(3),axes)

class NativePrisms:
    """18-node contraction with one law and per-quadrature material frames.

    gradients: (C,Q,18,3), reference dN/dX in m^-1
    weights: (C,Q), strictly positive reference volume weights in m^3
    H0: broadcastable to (C,Q,3,3), F at the increment origin minus I
    axes: broadcastable to (C,Q,3,3), reference material directions as columns
    density: broadcastable to (C,), fixed reference density in kg/m^3
    evaluate(du): (C,18,3) local displacement increments, metres

    H = H0 + sum_i (du_i - mean_i du_i) tensor grad_i.
    H0 stays fixed during an increment solve; a caller changing origin must
    explicitly transfer the complete current H to a new H0. This API does not
    redefine reference volume, reference mass, axes or gradients.
    """
    def __init__(self,gradients,weights,H0,axes,density,law=None):
        g=np.asarray(gradients,dtype=np.float64)
        if g.ndim!=4 or g.shape[-2:]!=(18,3) or not all(g.shape[:2]):raise ValueError('Expected nonempty (cells,quadrature,18,3) reference gradients')
        self.count,self.quadrature=g.shape[:2]
        self.gradients=_field(g,g.shape,'reference gradients')
        self.weights=_field(weights,g.shape[:2],'reference weights',positive=True)
        self.H0=_field(H0,(*g.shape[:2],3,3),'origin gradient H0')
        self.axes=_field(axes,self.H0.shape,'material axes')
        self.density=_field(density,(self.count,),'reference density',positive=True)
        self.reference_volume=self.weights.sum(axis=1)
        self.mass=self.reference_volume*self.density
        if not np.isfinite(self.reference_volume).all() or not np.isfinite(self.mass).all() or np.any(self.mass<=0):raise ValueError('Nonfinite or degenerate reference volume/mass')
        self.reference_volume.setflags(write=False);self.mass.setflags(write=False)
        self.material=NativeEnvelope(law)
        # Validate origins and frames without changing or relaxing them.
        self.material.evaluate_incremental(self.H0,self.axes)
    def evaluate(self,du):
        u=np.ascontiguousarray(du,dtype=np.float64)
        if u.shape!=(self.count,18,3) or not np.isfinite(u).all():raise ValueError('Expected finite (cells,18,3) local displacement increments')
        gradient=np.empty_like(u);energies=np.empty(self.count);volumes=np.empty(self.count)
        raw=np.empty((self.count,self.quadrature,6));bad=np.array([-1],dtype=np.int64)
        rc=self.material.prism_function(self.count,self.quadrature,_ptr(u),_ptr(self.H0),_ptr(self.axes),_ptr(self.gradients),_ptr(self.weights),_ptr(self.material.parameters),_ptr(gradient),_ptr(energies),_ptr(volumes),_ptr(raw),bad.ctypes.data_as(INT64))
        self.material._check(rc,bad)
        diagnostics=_diagnostics(raw);j=diagnostics['J']
        current_density=self.density[:,None]/j
        mean_density=self.mass/volumes
        if not np.isfinite(current_density).all() or not np.isfinite(mean_density).all():raise FloatingPointError('Nonfinite current-density bookkeeping')
        diagnostics.update(cell_energy_J=energies,cell_volume_m3=volumes,minimum_det_F=float(j.min()),maximum_det_F=float(j.max()),
            cell_reference_volume_m3=self.reference_volume.copy(),cell_material_mass_kg=self.mass.copy(),
            cell_mean_density_kg_m3=mean_density,quadrature_current_density_kg_m3=current_density,
            cell_axial_energy_J=np.sum(self.weights*raw[...,4],axis=1),cell_angular_energy_J=np.sum(self.weights*raw[...,5],axis=1),dissipated_energy_J=0.)
        return float(energies.sum()),gradient,diagnostics
