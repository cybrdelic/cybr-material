"""Exact rational algebra checks against the inspected Blender 4.3.2 basis."""
from fractions import Fraction as F
from pathlib import Path
import json

# Coefficients of u^0 through u^3, columns are the four Catmull control keys.
# Expanded directly from catmull_rom_basis_eval's n0,n1,n2,n3 times 1/2.
cycles=[
    [F(0),F(1),F(0),F(0)],
    [F(-1,2),F(0),F(1,2),F(0)],
    [F(1),F(-5,2),F(2),F(-1,2)],
    [F(-1,2),F(3,2),F(-3,2),F(1,2)],
]
bezier_controls=[
    [F(0),F(1),F(0),F(0)],
    [F(-1,6),F(1),F(1,6),F(0)],
    [F(0),F(1,6),F(1),F(-1,6)],
    [F(0),F(0),F(1),F(0)],
]
bernstein_power=[[1,0,0,0],[-3,3,0,0],[3,-6,3,0],[-1,3,-3,1]]
converted=[[sum(F(bernstein_power[k][j])*bezier_controls[j][i] for j in range(4))
            for i in range(4)] for k in range(4)]
assert converted==cycles
assert [sum(row) for row in cycles]==[1,0,0,0]
# Repeating P0=P1 yields first derivative (P2-P1)/2 at u=0.
assert [cycles[1][0]+cycles[1][1],cycles[1][2],cycles[1][3]]==[F(-1,2),F(1,2),0]
# Radius keys 0,1,1,0 overshoot their authored maximum at the midpoint.
r=[0,1,1,0]
mid=sum(sum(cycles[k][j]*r[j] for j in range(4))*F(1,2)**k for k in range(4))
assert mid==F(9,8)
receipt={'exact_rational_basis_match':True,
         'partition_of_unity_and_constant_radius_preserved':True,
         'repeated_endpoint_half_difference_tangent':True,
         'variable_radius_counterexample':{'keys':[0,1,1,0],'u':0.5,'interpolated':float(mid),'key_max':1},
         'scope':'Polynomial representation identity only; no scene execution or ray precision certification.'}
Path(__file__).with_name('cycles_basis_algebra_receipt.json').write_text(json.dumps(receipt,indent=2))
print(json.dumps(receipt,indent=2))
