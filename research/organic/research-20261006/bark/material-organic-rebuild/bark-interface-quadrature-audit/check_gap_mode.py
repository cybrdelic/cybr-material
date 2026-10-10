"""Analytic P2 face-gap mode invisible to two-point Gauss integration."""
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
K=1.e12 # Pa/m
amplitude=1.e-6 # m
area=.006*.006 # reference m^2
exact=2*K*amplitude**2*area/45
rows=[]
for order in (2,3,4,6):
    z,w=leggauss(order)
    gap=amplitude*(z*z-1/3)
    energy=.5*K*area/4*np.sum(w)*np.sum(w*gap*gap)
    rows.append(dict(order=order,energy_J=float(energy),relative_error=float(abs(energy/exact-1))))
assert rows[0]['energy_J']<exact*1.e-25
assert max(r['relative_error'] for r in rows[1:])<1.e-13
print(json.dumps(dict(exact_energy_J=exact,quadrature=rows,mode='gap(zeta)=a*(zeta^2-1/3)',claim='Undamaged normal spring on a planar P2 face. Not a complete nonlinear/contact quadrature qualification.'),indent=2))
