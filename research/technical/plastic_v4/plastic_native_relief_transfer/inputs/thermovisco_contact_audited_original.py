"""Reduced thermo-viscoelastic surface imprint model, SI units.

This is a transient numerical model, not an injection-flow solver or a validated
polymer model. Every material number below is an assumed generic value.

At each tool depth H, a Maxwell foundation of effective thickness ell carries
prescribed pressure p. Its irreversible dashpot displacement s evolves under
sigma = E*(u-s)/ell, u=min(H,s+p*ell/E). On tool contact, sigma relaxes and the
nonnegative obstacle reaction p-sigma carries the remainder. The retained
replica after unloading is s. A depth-resolved FV heat equation is coupled via
Arrhenius viscosity and the actual gap H-u. Gap + contact + first half-cell
thermal resistances are in series. An adiabatic polymer back face closes the
thermal domain. Exact isothermal mechanical substeps conserve work/energy;
backward-Euler FV thermal steps conserve heat including viscous dissipation.
No random input, final replication coefficient, or invented appearance field.
"""
from dataclasses import dataclass, asdict, replace
import sys
sys.path.insert(0,"/workspace/shared/material-causal-rebuild/audit")
from fixed_process_math import elastic_pressure_transition
import numpy as np

@dataclass(frozen=True)
class Parameters:
    # Assumed generic polymer/air/contact properties, not measured for the panel.
    density_kg_m3: float = 1050.
    heat_capacity_J_kgK: float = 1800.
    conductivity_W_mK: float = .20
    air_conductivity_W_mK: float = .030
    contact_h_W_m2K: float = 6000.
    maxwell_modulus_Pa: float = 3e8
    viscosity_reference_Pas: float = 3e5
    viscosity_reference_K: float = 503.15
    activation_energy_J_mol: float = 160000.
    gas_constant_J_molK: float = 8.314462618
    max_viscosity_Pas: float = 1e20
    deformation_layer_m: float = .00018
    thermal_depth_m: float = .004
    melt_temperature_K: float = 503.15
    mold_temperature_K: float = 333.15
    skin_diagnostic_temperature_K: float = 383.15
    pressure_Pa: float = 20e6
    pressure_hold_s: float = .50
    end_time_s: float = .50
    dt_s: float = .0000625
    depth_cells: int = 96
    mesh_stretch: float = 7.
    include_dissipation: bool = True


def depth_mesh(p):
    v=np.linspace(0.,1.,p.depth_cells+1)
    edges=p.thermal_depth_m*np.expm1(p.mesh_stretch*v)/np.expm1(p.mesh_stretch)
    dz=np.diff(edges); centers=(edges[:-1]+edges[1:])*.5
    overlap=np.maximum(0.,np.minimum(edges[1:],p.deformation_layer_m)-edges[:-1])
    assert abs(overlap.sum()-p.deformation_layer_m)<1e-12
    return edges,dz,centers,overlap/p.deformation_layer_m


def viscosity(T,p):
    logeta=np.log(p.viscosity_reference_Pas)+(p.activation_energy_J_mol/p.gas_constant_J_molK)*(1./T-1./p.viscosity_reference_K)
    return np.exp(np.minimum(logeta,np.log(p.max_viscosity_Pas)))


def mechanical_step(s,H,pressure,eta,dt,p):
    """Exact piecewise Maxwell/rigid-obstacle step at fixed eta and pressure.

    Input/output s is irreversible displacement. Elastic displacement and
    contact reaction are instantaneous quasistatic states. Returns exact
    dissipated energy per area and pressure work for this fixed-load step.
    """
    ell=p.deformation_layer_m; E=p.maxwell_modulus_Pa
    u0=np.minimum(H,s+pressure*ell/E)
    if pressure<=0:
        return s.copy(),s.copy(),np.zeros_like(s),np.zeros_like(s),np.zeros_like(s)
    rate=pressure*ell/eta
    time_hit=np.maximum(0.,H-pressure*ell/E-s)/rate
    pre=np.minimum(dt,time_hit)
    sp=np.minimum(H,s+rate*pre)
    post=dt-pre
    decay=np.exp(-E*post/eta)
    sn=sp+(H-sp)*(-np.expm1(-E*post/eta))
    sn=np.maximum(s,np.minimum(H,sn))
    u1=np.minimum(H,sn+pressure*ell/E)
    diss=pressure*pressure*ell*pre/eta + .5*E/ell*(H-sp)**2*(-np.expm1(-2*E*post/eta))
    work=pressure*(u1-u0)
    energy0=.5*E/ell*(u0-s)**2
    energy1=.5*E/ell*(u1-sn)**2
    balance=work-(energy1-energy0)-diss
    return sn,u1,diss,work,balance


def thermal_step(T,gap,heat_J_m2,dt,p,mesh=None):
    """Conservative backward-Euler finite volumes, vectorized over columns.

    T is (depth_cells,columns), gap is columns. Boundary flux uses the new cell
    temperature. Source integrates exactly to the mechanical dissipated heat.
    """
    edges,dz,z,source_weights=depth_mesh(p) if mesh is None else mesh
    cap=p.density_kg_m3*p.heat_capacity_J_kgK*dz
    edge_g=p.conductivity_W_mK/np.diff(z)
    face_h=1./(1./p.contact_h_W_m2K+gap/p.air_conductivity_W_mK+dz[0]/(2*p.conductivity_W_mK))
    n,m=T.shape
    diagonal=np.broadcast_to((cap/dt)[:,None],(n,m)).copy()
    diagonal[:-1]+=edge_g[:,None];diagonal[1:]+=edge_g[:,None]
    diagonal[0]+=face_h
    rhs=T*(cap/dt)[:,None]+source_weights[:,None]*heat_J_m2[None,:]/dt
    rhs[0]+=face_h*p.mold_temperature_K
    # Thomas elimination; no columns communicate (explicit reduced assumption).
    for j in range(1,n):
        factor=-edge_g[j-1]/diagonal[j-1]
        diagonal[j]-=factor*(-edge_g[j-1])
        rhs[j]-=factor*rhs[j-1]
    out=np.empty_like(T);out[-1]=rhs[-1]/diagonal[-1]
    for j in range(n-2,-1,-1):out[j]=(rhs[j]+edge_g[j]*out[j+1])/diagonal[j]
    Q=dt*face_h*(out[0]-p.mold_temperature_K)
    deltaU=np.sum(cap[:,None]*(out-T),axis=0)
    residual=deltaU+Q-heat_J_m2
    return out,Q,residual,face_h


def simulate(tool_depth_m,p=Parameters(),keep_history=True):
    H=np.asarray(tool_depth_m,dtype=np.float64).reshape(-1)
    if np.any(H<0) or not np.isfinite(H).all():raise ValueError('Tool depth must be finite and nonnegative')
    if p.dt_s<=0 or p.end_time_s<0:raise ValueError('Invalid time grid')
    if p.pressure_Pa<0 or p.contact_h_W_m2K<=0 or p.mold_temperature_K<=0:
        raise ValueError("Pressure must be nonnegative and thermal parameters positive")
    mesh=depth_mesh(p);edges,dz,z,weights=mesh
    T=np.full((p.depth_cells,len(H)),p.melt_temperature_K,dtype=np.float64)
    s=np.zeros_like(H);u=s.copy();Q_total=np.zeros_like(H);D_total=np.zeros_like(H)
    mechanical_D_total=np.zeros_like(H); W_total=np.zeros_like(H); previous_pressure=0.
    startup_work=np.zeros_like(H); release_work=np.zeros_like(H)
    energy_max=0.;mech_max=0.;t=0.;steps=0;history=[]; trace=[]
    probes=np.unique(np.linspace(0,len(H)-1,min(len(H),9)).astype(int))
    def snapshot(pressure):
        frac=np.divide(s,H,out=np.ones_like(H),where=H>0)
        skin=np.sum(weights[:,None]*(T<p.skin_diagnostic_temperature_K),axis=0)
        trace.append(dict(time_s=float(t), tool_depth_m=H[probes].tolist(), retained_displacement_m=s[probes].tolist(), loaded_displacement_m=u[probes].tolist(), surface_temperature_K=T[0,probes].tolist(), effective_viscosity_Pas=np.sum(weights[:,None]*viscosity(T[:,probes],p),axis=0).tolist()))
        history.append(dict(time_s=float(t),mean_retained_fraction=float(frac.mean()),mean_gap_um=float(np.mean(H-u)*1e6),mean_surface_K=float(T[0].mean()),max_surface_K=float(T[0].max()),mean_skin_fraction=float(skin.mean()),pressure_Pa=float(pressure),mean_mold_heat_J_m2=float(Q_total.mean()),mean_dissipation_J_m2=float(D_total.mean())))
    next_record=0.;record_interval=max(p.end_time_s/50,p.dt_s)
    while t<p.end_time_s-1e-13:
        pressure=p.pressure_Pa if t<p.pressure_hold_s-1e-13 else 0.
        dt=min(p.dt_s,p.end_time_s-t)
        if t<p.pressure_hold_s<t+dt:dt=p.pressure_hold_s-t
        # Uniform strain reduced layer: arithmetic thickness-average viscosity.
        eta=np.sum(weights[:,None]*viscosity(T,p),axis=0)
        if pressure != previous_pressure:
            transition=elastic_pressure_transition(s,H,previous_pressure,pressure,p.deformation_layer_m,p.maxwell_modulus_Pa)
            W_total+=transition["pressure_work_J_m2"]
            if pressure>previous_pressure:startup_work+=transition["pressure_work_J_m2"]
            else:release_work+=transition["pressure_work_J_m2"]
            previous_pressure=pressure
        s,u,diss,work,mech_res=mechanical_step(s,H,pressure,eta,dt,p)
        T,Q,res,hface=thermal_step(T,np.maximum(0,H-u),diss if p.include_dissipation else np.zeros_like(diss),dt,p,mesh)
        Q_total+=Q;D_total+=diss if p.include_dissipation else 0
        W_total+=work;mechanical_D_total+=diss
        energy_max=max(energy_max,float(np.max(np.abs(res))))
        mech_max=max(mech_max,float(np.max(np.abs(mech_res))))
        t+=dt;steps+=1
        if keep_history and t>=next_record-1e-13:snapshot(pressure);next_record+=record_interval
    pressure=p.pressure_Pa if t<=p.pressure_hold_s+1e-13 else 0.
    # Final unloaded geometry is s, not the elastic displacement while packed.
    final_gap=H-s
    # Explicit quasistatic unloading closes the elastic reservoir; it does not
    # heat the polymer. Dynamic removal would need a different model.
    transition=elastic_pressure_transition(s,H,previous_pressure,0.,p.deformation_layer_m,p.maxwell_modulus_Pa)
    release_work+=transition["pressure_work_J_m2"];W_total+=transition["pressure_work_J_m2"]
    mechanical_global=W_total-mechanical_D_total
    skin=np.sum(weights[:,None]*(T<p.skin_diagnostic_temperature_K),axis=0)
    thermal_storage=np.sum((p.density_kg_m3*p.heat_capacity_J_kgK*dz)[:,None]*(T-p.melt_temperature_K),axis=0)
    energy_global=thermal_storage+Q_total-D_total
    frac=np.divide(s,H,out=np.ones_like(H),where=H>0)
    return dict(tool_depth_m=H,retained_displacement_m=s,loaded_displacement_m=u,
      unloaded_gap_m=final_gap,retained_fraction=frac,temperature_K=T,
      surface_temperature_K=T[0],skin_fraction=skin,depth_edges_m=edges,
      extracted_heat_J_m2=Q_total,dissipated_heat_J_m2=D_total,history=history,trace=trace,
      pressure_work_J_m2=W_total, mechanical_dissipation_J_m2=mechanical_D_total,
      startup_pressure_work_J_m2=startup_work, release_pressure_work_J_m2=release_work,
      pressure_transition_model="Quasistatic pressure ramp at fixed dashpot state; no dynamic impact or release heating",
      parameters=asdict(p),diagnostics=dict(max_full_cycle_mechanical_residual_J_m2=float(np.max(abs(mechanical_global))),
      max_full_cycle_combined_residual_J_m2=float(np.max(abs(W_total-thermal_storage-Q_total-(mechanical_D_total-D_total)))),steps=steps,max_step_thermal_residual_J_m2=energy_max,
      max_step_mechanical_residual_J_m2=mech_max,max_global_thermal_residual_J_m2=float(np.max(np.abs(energy_global))),
      min_retained_fraction=float(frac.min()),max_retained_fraction=float(frac.max()),
      min_temperature_K=float(T.min()),max_temperature_K=float(T.max()),
      global_thermal_relative_residual=float(np.max(np.abs(energy_global))/max(float(np.max(np.abs(Q_total))),1.))))
