"""Single-droplet evaporation and heating model (0-D, documented stopping).

Governing equations (constant gas environment, no relative convection):

Mass (d^2 law):      q = D^2,   dq/dt = -K(T_d),
                     K = 8 k_g ln(1 + B_M) / (rho_l c_pg)
Thermal (lumped):    m_l c_pl dT_d/dt = h A_s (T_g - T_d) - Q_evap,
                     A_s = pi D^2 (sphere surface)

The K formula assumes unit Lewis number (rho_g D_v = k_g / c_pg): the
mass-transfer closure inherits the conduction coefficient. This
simplification is stated here and in the docs; a separate diffusivity would
replace k_g/c_pg.

With constant K the law integrates exactly to D(t) = sqrt(D0^2 - K t) with
analytical model extinction at t_e = D0^2/K. The COUPLED system is stopped
at a documented positive threshold q_stop = epsilon * q0: near zero mass the
lumped thermal time constant collapses (it scales with D^2) and the
Clapeyron surface state steepens — a genuinely singular terminal stage of
this closure. The stop is reported with its time, remaining mass fraction,
and reason, and is NEVER presented as complete evaporation or a validated
lifetime. Solver failures raise; they never masquerade as survival.

All symbols carry SI units; the model is 0-D (no droplet internals, no
relative velocity, no variable gas state).
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

__all__ = [
    "saturation_pressure_clapeyron",
    "saturation_vapor_mass_fraction",
    "spalding_mass_number",
    "d2_law_constant",
    "d2_law_diameter",
    "extinction_time",
    "lumped_heating_rhs",
    "droplet_rhs",
    "integrate_droplet",
]


def saturation_pressure_clapeyron(T_d, p_ref, T_ref, h_vap, R_vap):
    """Saturation pressure [Pa] from the integrated Clapeyron relation.

    ln(p_sat/p_ref) = -(h_vap/R_vap)(1/T_d - 1/T_ref), constant h_vap.
    Validity is the caller's responsibility; this function does NOT clamp.
    """
    if not np.isfinite(T_d) or T_d <= 0:
        raise ValueError(f"temperature must be finite positive, got {T_d}")
    return p_ref * np.exp(-h_vap / R_vap * (1.0 / T_d - 1.0 / T_ref))


def saturation_vapor_mass_fraction(T_d, p_ref, T_ref, h_vap, R_vap,
                                   p_ambient, mw_ratio):
    """Surface vapor mass fraction from the Clapeyron fit.

    Raises ValueError when p_sat reaches p_ambient — the boiling regime is
    outside this model's validity and is rejected, never clamped.
    """
    p_sat = saturation_pressure_clapeyron(T_d, p_ref, T_ref, h_vap, R_vap)
    if p_sat >= p_ambient:
        raise ValueError(
            f"saturation pressure {p_sat:.0f} Pa reached ambient {p_ambient:.0f} Pa "
            f"at T_d = {T_d:.1f} K: boiling regime is outside this model's validity"
        )
    chi = p_sat / p_ambient
    return chi * mw_ratio / (chi * mw_ratio + (1.0 - chi))


def spalding_mass_number(Y_surface, Y_gas):
    """Spalding mass-transfer number B_M = (Y_s - Y_g)/(1 - Y_s)."""
    if not 0.0 <= Y_surface < 1.0:
        raise ValueError(f"surface vapor mass fraction must be in [0,1), got {Y_surface}")
    if not 0.0 <= Y_gas < 1.0:
        raise ValueError(f"gas vapor mass fraction must be in [0,1), got {Y_gas}")
    return (Y_surface - Y_gas) / (1.0 - Y_surface)


def d2_law_constant(rho_l, k_g, cp_g, B_M):
    """d^2-law evaporation constant K [m^2/s] (unit-Lewis-number closure).

    K = 8 k_g ln(1 + B_M) / (rho_l c_pg): the mass diffusivity is taken as
    rho_g D_v = k_g / c_pg (Le = 1), inheriting the conduction coefficient.
    """
    if not (np.isfinite(rho_l) and np.isfinite(k_g) and np.isfinite(cp_g)):
        raise ValueError("rho_l, k_g, cp_g must be finite")
    if rho_l <= 0 or k_g <= 0 or cp_g <= 0:
        raise ValueError("rho_l, k_g, cp_g must be positive")
    if not np.isfinite(B_M) or B_M <= 0:
        raise ValueError("B_M must be finite positive for evaporation (Y_s > Y_g)")
    return 8.0 * k_g * np.log1p(B_M) / (rho_l * cp_g)


def d2_law_diameter(t, D0, K):
    """Exact constant-K d^2-law diameter D(t) = sqrt(D0^2 - K t) [m]."""
    if not np.isfinite(D0) or D0 <= 0 or not np.isfinite(K) or K <= 0:
        raise ValueError("D0 and K must be finite positive")
    return np.sqrt(np.maximum(D0 * D0 - K * t, 0.0))


def extinction_time(D0, K):
    """Analytical model extinction time t_e = D0^2/K [s] (constant K)."""
    if not np.isfinite(D0) or D0 <= 0 or not np.isfinite(K) or K <= 0:
        raise ValueError("D0 and K must be finite positive")
    return D0 * D0 / K


def lumped_heating_rhs(T_d, D, T_g, h, m_l, cp_l, Q_evap):
    """Lumped droplet temperature rate dT_d/dt [K/s].

    m_l c_pl dT_d/dt = h * A_s * (T_g - T_d) - Q_evap with the SPHERE
    surface area A_s = pi D^2 (h in W/(m^2 K)). Q_evap [W] is the latent
    sink; the sensible enthalpy of the departing vapor is neglected
    (documented assumption of this closure).
    """
    if not np.isfinite(m_l) or m_l <= 0:
        raise ValueError("droplet mass must be finite positive")
    if not np.isfinite(cp_l) or cp_l <= 0:
        raise ValueError("cp_l must be finite positive")
    return (h * np.pi * D * D * (T_g - T_d) - Q_evap) / (m_l * cp_l)


def droplet_rhs(D, T_d, rho_l, k_g, cp_g, cp_l, h_fn, T_g, Y_surface_fn,
                Y_gas, Q_evap_fn):
    """Right-hand side d(D, T_d)/dt of the coupled 0-D system.

    ``h_fn(D)`` returns the convective coefficient [W/(m^2 K)] — pass a
    constant callable or the Nu = 2 conduction closure h(D) = 2 k_g / D.
    Returns (0, 0) at D <= 0.
    """
    if D <= 0.0:
        return 0.0, 0.0
    m_l = rho_l * np.pi * D**3 / 6.0
    B_M = spalding_mass_number(Y_surface_fn(T_d), Y_gas)
    K = d2_law_constant(rho_l, k_g, cp_g, B_M)
    dD = -K / (2.0 * D)
    dT = lumped_heating_rhs(T_d, D, T_g, h_fn(D), m_l, cp_l, Q_evap_fn(T_d, D))
    if not (np.isfinite(dD) and np.isfinite(dT)):
        raise FloatingPointError(
            f"non-finite rate at D={D:.3e} m, T_d={T_d:.1f} K; callback outputs "
            "left the supported regime"
        )
    return dD, dT


def _validate_inputs(D0, T_d0, t_end, rho_l, k_g, cp_g, cp_l, T_g, Y_gas):
    if not np.isfinite(D0) or D0 <= 0:
        raise ValueError(f"D0 must be finite positive, got {D0}")
    if not np.isfinite(T_d0) or T_d0 <= 0:
        raise ValueError(f"T_d0 must be finite positive, got {T_d0}")
    if not np.isfinite(t_end) or t_end <= 0:
        raise ValueError(f"t_end must be finite positive, got {t_end}")
    for name, value in (("rho_l", rho_l), ("k_g", k_g), ("cp_g", cp_g), ("cp_l", cp_l)):
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite positive, got {value}")
    if not np.isfinite(T_g) or T_g <= 0:
        raise ValueError(f"T_g must be finite positive, got {T_g}")
    if not 0.0 <= Y_gas < 1.0:
        raise ValueError(f"Y_gas must be in [0,1), got {Y_gas}")


def integrate_droplet(D0, T_d0, t_end, rho_l, k_g, cp_g, cp_l, h_fn, T_g,
                      Y_surface_fn, Y_gas, Q_evap_fn, epsilon=1e-3,
                      max_step=None, dense=True):
    """Integrate the coupled (q, T_d) system to a documented positive stop.

    State: q = D^2 with dq/dt = -K(T_d) computed directly (no 1/D
    arithmetic). The run TERMINATES at the documented threshold
    q_stop = epsilon * q0 — a threshold-stopping event near the singular
    endpoint of this closure — and reports the stopping time, remaining
    mass fraction, and reason. It is not complete evaporation; the
    analytical constant-K helper ``extinction_time`` carries that model
    quantity separately.

    Returns a dict with ``t``, ``D``, ``T_d``, ``dense`` (a callable
    sol(t) -> (q, T_d) when dense=True), and the stopping record:
    ``stop_time`` (None if t_end was reached), ``remaining_mass_fraction``,
    ``stop_reason`` ("threshold" or "t_end"), ``solver_success``.
    Raises RuntimeError when the integrator fails before either stop.
    """
    _validate_inputs(D0, T_d0, t_end, rho_l, k_g, cp_g, cp_l, T_g, Y_gas)
    if not 0.0 < epsilon < 1.0:
        raise ValueError(f"epsilon must be in (0,1), got {epsilon}")
    q0 = D0 * D0
    q_stop = epsilon * q0

    def rhs(t, y):
        q, T_d = y
        D = np.sqrt(max(q, 0.0))
        _, dT = droplet_rhs(D, T_d, rho_l, k_g, cp_g, cp_l, h_fn, T_g,
                            Y_surface_fn, Y_gas, Q_evap_fn)
        B_M = spalding_mass_number(Y_surface_fn(T_d), Y_gas)
        dq = -d2_law_constant(rho_l, k_g, cp_g, B_M)
        return [dq, dT]

    def threshold(t, y):
        return y[0] - q_stop

    threshold.terminal = True
    threshold.direction = -1

    sol = solve_ivp(
        rhs,
        (0.0, t_end),
        [q0, T_d0],
        method="DOP853",
        rtol=1e-11,
        atol=1e-16,
        events=threshold,
        dense_output=dense,
        **({"max_step": max_step} if max_step is not None else {}),
    )
    if not sol.success and not len(sol.t_events[0]):
        raise RuntimeError(
            f"droplet integration failed at t={sol.t[-1]:.6e} s: {sol.message}"
        )
    hit = len(sol.t_events[0]) > 0
    stop_time = float(sol.t_events[0][0]) if hit else None
    q_final = float(sol.y[0, -1])
    return {
        "t": sol.t,
        "D": np.sqrt(np.maximum(sol.y[0], 0.0)),
        "T_d": sol.y[1],
        "dense": sol.sol if dense else None,
        "stop_time": stop_time,
        "remaining_mass_fraction": float((q_final / q0) ** 1.5),
        "stop_reason": "threshold" if hit else "t_end",
        "solver_success": bool(sol.success or hit),
    }
