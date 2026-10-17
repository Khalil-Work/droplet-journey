# Model and method

## Problem

A single spherical droplet evaporating in a constant, quiescent gas
environment — 0-D, lumped, no relative convection. State: diameter D [m]
and droplet temperature T_d [K].

## Equations

**Mass — the d² law with Spalding driving force.** Conduction-limited
evaporation with the **unit-Lewis-number closure** (mass diffusivity taken
as ρ_g D_v = k_g / c_pg):

    dD/dt = −K / (2 D),   K = 8 k_g ln(1 + B_M) / (ρ_l c_pg),

with gas conductivity k_g [W/mK], liquid density ρ_l [kg/m³], gas specific
heat c_pg [J/kgK], and the Spalding mass-transfer number

    B_M = (Y_s − Y_g) / (1 − Y_s),

where Y_s, Y_g are vapor mass fractions at the surface and in the gas.
This is not an unrestricted classical mass-transfer formula: the Le = 1
simplification is inherited from the conduction coefficient, and a separate
diffusivity would replace k_g/c_pg. With constant K the law integrates
exactly to D(t) = √(D0² − Kt) with analytical model extinction at
t_e = D0²/K.

**Surface state.** Y_s(T_d) from the integrated Clapeyron relation with a
reference anchor (p_ref, T_ref) and constant latent heat h_vap:

    ln(p_sat/p_ref) = −(h_vap/R_vap)(1/T_d − 1/T_ref),  Y_s from p_sat.

When p_sat reaches the ambient pressure the model RAISES — the boiling
regime is outside its validity and is rejected, never clamped.

**Thermal — lumped capacitance with evaporation cooling.**

    m_l c_pl dT_d/dt = h A_s (T_g − T_d) − Q_evap,
    A_s = π D² (sphere surface),  m_l = ρ_l π D³ / 6,

with the convective coefficient h [W/m²K] supplied as a callback h(D); the
demo case uses the Nu = 2 conduction-limit closure h(D) = 2 k_g / D.
Q_evap = ṁ·h_fg is the latent sink with ṁ = ρ_l π D K / 4 consistent with
the d² law. The sensible enthalpy of the departing vapor is neglected —
part of this closure's stated assumptions.

## Numerical method

The coupled system is integrated in the state (q, T_d) with q = D²:

    dq/dt = −K(T_d)  — computed directly, with no 1/D arithmetic,

using DOP853 (rtol 1e-11) and a dense solution for uniform post-sampling.

**Documented positive stopping threshold.** The lumped thermal time
constant scales with D² and the Clapeyron surface state steepens as mass
vanishes — a genuinely singular terminal stage of this closure. The
coupled run therefore TERMINATES at q_stop = ε·q0 (ε = 1e-3 by default) and
reports the stopping time, remaining mass fraction (ε^1.5), and reason.
This is a threshold stop, never complete evaporation: the analytical
constant-K helper `extinction_time` carries the model-extinction quantity
separately. Solver failures raise `RuntimeError`; they never appear as
survival.

## Event semantics

Threshold stopping times and the analytical model extinction time are
**model predictions under the stated constant conditions** — well-defined
ODE quantities — and are never presented as validated physical evaporation
lifetimes.
