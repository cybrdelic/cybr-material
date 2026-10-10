"""Exact rational audit of retained float32 owner0/native span36.

The uniform Catmull-Rom power coefficients are dyadic rationals. Define
S=gamma' dot gamma', T=|gamma' cross gamma''|^2, H=S^3-r^2*T.
S>0 and H>0 are exactly regularity and strict K*r<1 at each parameter.
Rational Bernstein subdivision proves polynomial signs continuously; sampled
analytic values are diagnostic only. No geometry, solver or gate is changed.
"""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import math
import sys
import time
import numpy as np
from physical_length_window import ROOT,UNIT
sys.path.insert(0,str(UNIT/'native-qualification'))
import native_curve_gate as gate

def add(a,b):
    n=max(len(a),len(b));out=[F(0)]*n
    for i,v in enumerate(a):out[i]+=v
    for i,v in enumerate(b):out[i]+=v
    return out

def mul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out

def evaluate(p,t):
    value=F(0)
    for coefficient in p[::-1]:value=value*t+coefficient
    return value

def bernstein(power):
    degree=len(power)-1
    return [sum((power[k]*F(math.comb(i,k),math.comb(degree,k)) for k in range(i+1)),F(0)) for i in range(degree+1)]

def split(controls):
    levels=[controls]
    while len(levels[-1])>1:levels.append([(a+b)/2 for a,b in zip(levels[-1][:-1],levels[-1][1:])])
    return [level[0] for level in levels],[level[-1] for level in levels[::-1]]

def restrict(controls,lo,hi):
    # Dyadic intervals produced by binary subdivision only.
    current=controls;left=F(0);right=F(1)
    while (left,right)!=(lo,hi):
        a,b=split(current);mid=(left+right)/2
        if hi<=mid:current=a;right=mid
        elif lo>=mid:current=b;left=mid
        else:raise ValueError('interval is not a binary subdivision leaf')
    return current

def prove_positive(controls,max_depth=18,max_work=8192):
    pending=[(F(0),F(1),controls,0)];leaves=[];unresolved=[];work=0
    while pending:
        lo,hi,c,depth=pending.pop();work+=1
        if min(c)>0:
            leaves.append((lo,hi,min(c),max(c)));continue
        if max(c)<0:return {'status':'negative_interval','lo':lo,'hi':hi,'controls':c,'depth':depth,'work':work}
        if work>=max_work:
            return {'status':'unresolved','lo':lo,'hi':hi,'controls':c,'depth':depth,'work':work}
        if depth>=max_depth:
            # A sign-changing polynomial has unresolved root leaves. Continue
            # the finite search to find a strict negative interval elsewhere.
            unresolved.append((lo,hi,c,depth));continue
        a,b=split(c);mid=(lo+hi)/2
        pending.extend([(mid,hi,b,depth+1),(lo,mid,a,depth+1)])
    if unresolved:
        lo,hi,c,depth=unresolved[0]
        return {'status':'unresolved','lo':lo,'hi':hi,'controls':c,'depth':depth,'work':work}
    return {'status':'positive_entire_interval','lower':min(x[2] for x in leaves),
            'upper':max(x[3] for x in leaves),'leaves':len(leaves),'work':work}

def rational_report(value):
    if isinstance(value,F):return str(value)
    if isinstance(value,dict):return {k:rational_report(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [rational_report(v) for v in value]
    return value

def trace_gate(body,radius):
    pending=gate.native_beziers(body);span=np.arange(len(pending));t0=np.zeros(len(pending));t1=np.ones(len(pending))
    pad=256*gate.EPS*max(1e-6,float(np.abs(body).max()))*19
    rows=[]
    for depth in range(17):
        speed,kappa=gate.continuous_bounds(pending,pad)
        good=np.isfinite(speed)&np.isfinite(kappa)&(speed>0)&(kappa*radius<1)
        if depth<3:good[:]=False
        rows.append({'depth':depth,'pending':len(pending),'accepted':int(good.sum()),'unresolved':int((~good).sum())})
        if np.all(good):return {'status':'pass','trace':rows}
        if depth==16:
            index=int(np.flatnonzero(~good)[0]);u=F(float(t0[index]));v=F(float(t1[index]))
            unpadded_speed,unpadded_kappa=gate.continuous_bounds(pending[index:index+1],0.)
            return {'status':'same_gate_failure_reproduced','native_span':int(span[index]),
                    'parameter_interval':[str(u),str(v)],'depth':depth,
                    'leaf_float64_bezier_controls':pending[index].tolist(),'position_pad_m':pad,
                    'leaf_parameter_speed_lower_m':float(speed[index]),
                    'native_parameter_speed_lower_m':float(speed[index]*2**depth),
                    'padded_kappa_radius_upper':float(kappa[index]*radius),
                    'unpadded_kappa_radius_upper':float(unpadded_kappa[0]*radius),
                    'trace':rows},u,v
        use=~good;mid=(t0[use]+t1[use])/2
        t0=np.stack((t0[use],mid),axis=1).ravel();t1=np.stack((mid,t1[use]),axis=1).ravel()
        span=np.repeat(span[use],2);pending=gate.split_beziers(pending[use])

def run():
    start=time.monotonic()
    # Analytic parabola: S=1+4u^2, T=4, so focal equality is exactly H(0)=0.
    analytic_S=[F(1),F(0),F(4)]
    analytic_cube=mul(mul(analytic_S,analytic_S),analytic_S)
    assert prove_positive(bernstein(add(analytic_cube,[F(-1)])))['status']=='unresolved'
    assert prove_positive(bernstein(add(analytic_cube,[F(-36,25)])))['status']=='negative_interval'
    assert prove_positive(bernstein(add(analytic_cube,[F(-1,4)])))['status']=='positive_entire_interval'
    frozen=UNIT/'revisions/r9_physical_resume/arrays/loop_968_partial.npz'
    candidate=ROOT/'arrays/trial_968_prefix_native.npz'
    state=np.load(frozen);new=np.load(candidate)
    old=state['world'].astype('f4')[:,0,:]
    keys=new['body'][0];points=old[35:39]
    assert np.array_equal(points,keys[35:39])
    assert np.array_equal(gate.native_beziers(old)[36],gate.native_beziers(keys)[36])
    radius=F(float(new['body_radius']))
    exact=[[F(float(v)) for v in row] for row in points]
    p0,p1,p2,p3=exact
    power=[p1,[(b-a)/2 for a,b in zip(p0,p2)],
           [a-F(5,2)*b+2*c-d/2 for a,b,c,d in zip(p0,p1,p2,p3)],
           [-a/2+F(3,2)*b-F(3,2)*c+d/2 for a,b,c,d in zip(p0,p1,p2,p3)]]
    D=[[power[1][axis],2*power[2][axis],3*power[3][axis]] for axis in range(3)]
    E=[[2*power[2][axis],6*power[3][axis]] for axis in range(3)]
    S=[F(0)]
    for axis in range(3):S=add(S,mul(D[axis],D[axis]))
    cross=[]
    for a,b in ((1,2),(2,0),(0,1)):
        cross.append(add(mul(D[a],E[b]),[-v for v in mul(D[b],E[a])]))
    T=[F(0)]
    for value in cross:T=add(T,mul(value,value))
    H=add(mul(mul(S,S),S),[-radius*radius*v for v in T])
    speed_proof=prove_positive(bernstein(S))
    nonfold_proof=prove_positive(bernstein(H))
    trace,lo,hi=trace_gate(keys,float(radius))
    assert trace['native_span']==36
    gate_speed_bounds=restrict(bernstein(S),lo,hi)
    gate_cross_bounds=restrict(bernstein(T),lo,hi)
    gate_nonfold_bounds=restrict(bernstein(H),lo,hi)
    witness_parameter=(lo+hi)/2
    witness_S=evaluate(S,witness_parameter);witness_T=evaluate(T,witness_parameter);witness_H=evaluate(H,witness_parameter)
    witness={'parameter':str(witness_parameter),'speed_squared_exact':str(witness_S),
             'nonfold_polynomial_exact':str(witness_H),'speed_m_per_native_parameter':math.sqrt(float(witness_S)),
             'kappa_radius':math.sqrt(float(radius*radius*witness_T/witness_S**3)),
             'actual_excess_curvature_proved_at_exact_parameter':witness_S>0 and witness_H<0}
    if nonfold_proof['status']=='negative_interval':
        a,b=nonfold_proof['lo'],nonfold_proof['hi']
        sa=restrict(bernstein(S),a,b);ta=restrict(bernstein(T),a,b)
        hc=nonfold_proof['controls']
        lower=max(math.sqrt(max(0,float(radius*radius*min(ta)/max(sa)**3))),
                  math.sqrt(float(1-max(hc)/max(sa)**3)))
        upper=math.sqrt(float(radius*radius*max(ta)/min(sa)**3)) if min(sa)>0 else None
        classification='actual_native_excess_curvature_interval_proved'
        correction='Add a native curvature constraint at the adaptive key transition around stations35..38 before accepting/exporting a continuation. Do not copy adaptive section stations directly into uniform Catmull-Rom keys without native qualification. Stage any tangent-compatible refit or redistribution separately; changing raw keys alone requires rechecking all contact/containment/root constraints.'
        bad_interval={'parameter_interval':[str(a),str(b)],'speed_lower_m':math.sqrt(float(min(sa))),
                      'speed_upper_m':math.sqrt(float(max(sa))),
                      'kappa_radius_lower':lower,'kappa_radius_upper':upper,
                      'nonfold_polynomial_upper_exact':str(max(hc)),
                      'strict_nonfold_polynomial_bernstein_all_negative':True}
    elif nonfold_proof['status']=='positive_entire_interval':
        classification='exact_native_span_regular_and_nonfolding_gate_bound_too_loose'
        correction='Replace or supplement the lossy absolute-position curvature bound with an audited exact-polynomial or translation-invariant relative-control bound for stored float32 cubics. Keep strictK*r<1 and regularity gates; no optimization or key mutation is justified by this span.'
        bad_interval=None
    else:
        classification='exact_polynomial_audit_unresolved'
        correction='Continue exact polynomial isolation on this one retained span; no geometry change or promotion is justified yet.'
        bad_interval=None
    grid=np.linspace(0,1,2049)
    pf=np.asarray([[float(v) for v in p] for p in power])
    derivative=pf[1]+2*grid[:,None]*pf[2]+3*grid[:,None]**2*pf[3]
    acceleration=2*pf[2]+6*grid[:,None]*pf[3]
    sampled_speed=np.linalg.norm(derivative,axis=1)
    sampled_kappa=np.linalg.norm(np.cross(derivative,acceleration),axis=1)/sampled_speed**3*float(radius)
    maximum_parameter=F(int(np.argmax(sampled_kappa)),2048)
    maximum_S=evaluate(S,maximum_parameter);maximum_T=evaluate(T,maximum_parameter);maximum_H=evaluate(H,maximum_parameter)
    exact_peak={'parameter':str(maximum_parameter),'speed_squared_exact':str(maximum_S),
                'nonfold_polynomial_exact':str(maximum_H),
                'speed_m_per_native_parameter':math.sqrt(float(maximum_S)),
                'kappa_radius':math.sqrt(float(radius*radius*maximum_T/maximum_S**3)),
                'strict_excess_curvature_point_proved':maximum_S>0 and maximum_H<0}
    assert speed_proof['status']=='positive_entire_interval'
    assert nonfold_proof['status']=='negative_interval' and exact_peak['strict_excess_curvature_point_proved']
    index=min(int(np.argmax(sampled_kappa)),2047)
    peak_lo,peak_hi=F(index,2048),F(index+1,2048)
    peak_h=restrict(bernstein(H),peak_lo,peak_hi)
    peak_s=restrict(bernstein(S),peak_lo,peak_hi)
    peak_t=restrict(bernstein(T),peak_lo,peak_hi)
    assert max(peak_h)<0 and min(peak_s)>0
    peak_interval={'parameter_interval':[str(peak_lo),str(peak_hi)],
                   'speed_lower_m':math.sqrt(float(min(peak_s))),
                   'speed_upper_m':math.sqrt(float(max(peak_s))),
                   'kappa_radius_lower':math.sqrt(float(1-max(peak_h)/max(peak_s)**3)),
                   'kappa_radius_upper':math.sqrt(float(radius*radius*max(peak_t)/min(peak_s)**3)),
                   'exact_speed_squared_bernstein':rational_report(peak_s),
                   'exact_nonfold_polynomial_bernstein':rational_report(peak_h),
                   'all_exact_nonfold_controls_negative':True}
    assert peak_interval['kappa_radius_lower']>1
    result={'classification':classification,'owner':0,'native_span':36,
        'frozen_native_controls_equal_rejected_candidate':True,'retained_history_outside_last_window_all_fibre_motion':True,
        'stored_float32_points':points.tolist(),'stored_float32_point_bits_hex':[[f'{int(v):08x}' for v in row] for row in points.view('u4')],
        'guide_station_parameters_35_through38_m':state['stations'][35:39].tolist(),
        'adjacent_parameter_steps_35_through38_m':np.diff(state['stations'][35:39]).tolist(),
        'adjacent_stored_key_chord_lengths_m':np.linalg.norm(np.diff(points.astype(float),axis=0),axis=1).tolist(),
        'exact_native_power_coefficients':rational_report(power),
        'exact_radius_m':str(radius),'exact_speed_squared_power':rational_report(S),
        'exact_cross_product_squared_power':rational_report(T),'exact_nonfold_polynomial_power':rational_report(H),
        'continuous_speed_proof':rational_report(speed_proof),'continuous_nonfold_proof':rational_report(nonfold_proof),
        'whole_native_span_speed_bounds_m':[math.sqrt(float(speed_proof['lower'])),math.sqrt(float(speed_proof['upper']))],
        'strict_excess_curvature_interval':bad_interval,'gate_reproduction':trace,'exact_point_witness':witness,
        'exact_excess_curvature_point':exact_peak,
        'strong_excess_curvature_interval':peak_interval,
        'gate_first_unresolved_leaf_nonfolding_exactly_proved':min(gate_speed_bounds)>0 and min(gate_nonfold_bounds)>0,
        'gate_first_leaf_exact_nonfold_bernstein_bounds':rational_report([min(gate_nonfold_bounds),max(gate_nonfold_bounds)]),
        'gate_leaf_exact_speed_squared_bernstein_bounds':rational_report([min(gate_speed_bounds),max(gate_speed_bounds)]),
        'gate_leaf_exact_cross_squared_bernstein_bounds':rational_report([min(gate_cross_bounds),max(gate_cross_bounds)]),
        'sampled_diagnostics_only':{'samples':len(grid),'minimum_speed_m':float(sampled_speed.min()),
                                    'maximum_kappa_radius':float(sampled_kappa.max()),
                                    'maximum_parameter':float(grid[np.argmax(sampled_kappa)]),'sampled_pass_is_not_a_proof':True},
        'minimal_upstream_correction':correction,'no_solver':True,'frozen_states_modified':False,
        'exact_polynomial_source_assertions':['analytic parabola strict focal equality fails closed',
            'analytic parabola radius0.6 has negative nonfold interval',
            'analytic parabola radius0.25 has strictly positive nonfold polynomial',
            'exact native span controls equal frozen history',
            'exact native entire-span speed positive',
            'exact native negative nonfold interval and exact point both prove excess curvature'],
        'selected_r5_promoted':False,'elapsed_seconds':time.monotonic()-start,
        'source_sha256':hashlib.file_digest(Path(__file__).open('rb'),'sha256').hexdigest()}
    path=ROOT/'receipts/retained_span36_diagnosis.json'
    path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    np.savez_compressed(ROOT/'arrays/retained_span36_coefficients.npz',stored_points=points,
                        power_float64=pf,bezier_float64=gate.native_beziers(old)[36],
                        body_radius=new['body_radius'])
    print(json.dumps({'classification':classification,'speed_proof':speed_proof['status'],
                      'nonfold_proof':nonfold_proof['status'],'witness':witness,
                      'bad_interval':bad_interval,'sampled_max_kappa_radius':float(sampled_kappa.max()),
                      'exact_excess_curvature_point':exact_peak,
                      'strong_interval_summary':{k:v for k,v in peak_interval.items() if not k.startswith('exact_')},
                      'gate_first_leaf_nonfolding_exactly_proved':result['gate_first_unresolved_leaf_nonfolding_exactly_proved'],
                      'seconds':result['elapsed_seconds']},indent=2))

if __name__=='__main__':run()
