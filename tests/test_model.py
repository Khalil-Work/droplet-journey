"""Verification tests for droplet_journey.

Each test pins a distinct property:
- Spalding number limits and rejection of invalid fractions;
- the d^2-law constant against a hand-computed value (unit-Lewis closure);
- the exact constant-K solution, its linearity in d^2, and extinction time;
- Clapeyron fit reproduces its anchor; boiling regime is REJECTED, not clamped;
- SPHERE surface area scaling: at fixed h the initial heating rate halves
  when D doubles (mass ~ D^3, area ~ D^2);
- the coupled integrator reproduces the exact constant-K d^2 law when the
  surface state and heating are frozen;
- threshold stopping: documented positive stop with remaining mass and
  reason — never a silent zero-diameter completion;
- independent tight-tolerance reference (Radau, rtol 1e-13) agrees with the
  production integrator on the coupled case over the valid interval;
- input validation and solver-failure propagation.
"""

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from droplet_journey import (
    d2_law_constant,
    d2_law_diameter,
    droplet_rhs,
    extinction_time,
    integrate_droplet,
    lumped_heating_rhs,
    saturation_pressure_clapeyron,
    saturation_vapor_mass_fraction,
    spalding_mass_number,
)


def frozen_run(D0, t_end, epsilon=1e-6):
    """Frozen-surface run: constant Y_s, no thermal coupling."""
    return integrate_droplet(
        D0=D0, T_d0=300.0, t_end=t_end, rho_l=1000.0, k_g=0.026, cp_g=1010.0,
        cp_l=4180.0, h_fn=lambda D: 0.0, T_g=300.0,
        Y_surface_fn=lambda T: 0.2, Y_gas=0.0, Q_evap_fn=lambda T, D: 0.0,
        epsilon=epsilon,
    )


def test_spalding_number_values_and_limits():
    assert spalding_mass_number(0.3, 0.0) == pytest.approx(0.3 / 0.7)
    assert spalding_mass_number(0.2, 0.2) == 0.0
    with pytest.raises(ValueError, match="surface"):
        spalding_mass_number(1.0, 0.0)
    with pytest.raises(ValueError, match="gas"):
        spalding_mass_number(0.2, -0.1)


def test_d2_law_constant_matches_hand_computation():
    # rho_l=1000, k_g=0.026, cp_g=1010, B_M=e-1 => ln(1+B)=1 exactly.
    K = d2_law_constant(1000.0, 0.026, 1010.0, np.e - 1.0)
    assert K == pytest.approx(8.0 * 0.026 * 1.0 / (1000.0 * 1010.0), rel=1e-14)
    with pytest.raises(ValueError, match="B_M"):
        d2_law_constant(1000.0, 0.026, 1010.0, 0.0)
    with pytest.raises(ValueError, match="finite"):
        d2_law_constant(1000.0, 0.026, 1010.0, float("nan"))


def test_exact_constant_k_solution_and_extinction():
    D0, K = 100e-6, 1.0e-8
    assert extinction_time(D0, K) == pytest.approx(D0 * D0 / K)
    ts = np.linspace(0.0, extinction_time(D0, K), 101)
    D = d2_law_diameter(ts, D0, K)
    np.testing.assert_allclose(D * D, D0 * D0 - K * ts, rtol=1e-12, atol=1e-18)
    with pytest.raises(ValueError, match="positive"):
        extinction_time(-1.0, K)


def test_clapeyron_anchor_and_boiling_rejection():
    p_ref, T_ref, h_vap, R_vap = 3533.0, 300.0, 2.4e6, 461.5
    assert saturation_pressure_clapeyron(T_ref, p_ref, T_ref, h_vap, R_vap) == pytest.approx(p_ref)
    # Boiling regime is rejected, never clamped.
    with pytest.raises(ValueError, match="boiling"):
        saturation_vapor_mass_fraction(500.0, p_ref, T_ref, h_vap, R_vap,
                                       p_ambient=101325.0, mw_ratio=0.622)
    with pytest.raises(ValueError, match="positive"):
        saturation_pressure_clapeyron(0.0, p_ref, T_ref, h_vap, R_vap)


def test_sphere_area_scaling_of_initial_heating_rate():
    # At fixed h the initial dT/dt must halve when D doubles
    # (area ~ D^2, mass ~ D^3 => time constant ~ D^2).
    h, T_g, cp_l = 100.0, 320.0, 4180.0
    rho_l = 1000.0
    rates = []
    for D in (50e-6, 100e-6):
        m = rho_l * np.pi * D**3 / 6.0
        rates.append(lumped_heating_rhs(300.0, D, T_g, h, m, cp_l, 0.0))
    assert rates[0] == pytest.approx(2.0 * rates[1], rel=1e-12)
    with pytest.raises(ValueError, match="mass"):
        lumped_heating_rhs(300.0, 50e-6, 320.0, 100.0, 0.0, cp_l, 0.0)


def test_coupled_integrator_reproduces_frozen_d2_law():
    D0 = 80e-6
    K = d2_law_constant(1000.0, 0.026, 1010.0, spalding_mass_number(0.2, 0.0))
    t_e = extinction_time(D0, K)
    result = frozen_run(D0, t_e * 0.7)
    exact = d2_law_diameter(result["t"], D0, K)
    np.testing.assert_allclose(result["D"], exact, rtol=1e-8, atol=1e-15)
    np.testing.assert_allclose(result["T_d"], 300.0, atol=1e-9)
    assert result["stop_reason"] == "t_end"
    assert result["solver_success"]


def test_threshold_stopping_reports_mass_and_reason():
    D0 = 60e-6
    K = d2_law_constant(1000.0, 0.026, 1010.0, spalding_mass_number(0.2, 0.0))
    t_e = extinction_time(D0, K)
    result = frozen_run(D0, t_e * 10.0, epsilon=1e-3)
    assert result["stop_reason"] == "threshold"
    assert result["stop_time"] is not None
    assert result["stop_time"] == pytest.approx(t_e * (1 - 1e-3), rel=1e-6)
    # Remaining mass fraction is epsilon^1.5 = (1e-3)^1.5 ~ 3.16e-5.
    assert result["remaining_mass_fraction"] == pytest.approx(1e-3**1.5, rel=1e-3)
    assert result["D"][-1] > 0.0  # documented positive stop, never zero


def coupled_case(t_end=0.05):
    """Shared coupled case: Clapeyron surface, Nu=2 closure h(D)=2 k_g / D."""
    rho_l, k_g, cp_g, cp_l = 997.0, 0.026, 1010.0, 4180.0
    p_ref, T_ref, h_vap, R_vap, mw_ratio = 3533.0, 300.0, 2.26e6, 461.5, 0.622

    def Y_s(T):
        return saturation_vapor_mass_fraction(T, p_ref, T_ref, h_vap, R_vap,
                                               p_ambient=101325.0, mw_ratio=mw_ratio)

    def Q_evap(T, D):
        B_M = spalding_mass_number(Y_s(T), 0.0)
        K = d2_law_constant(rho_l, k_g, cp_g, B_M)
        return rho_l * np.pi * D * K / 4.0 * h_vap

    kwargs = dict(rho_l=rho_l, k_g=k_g, cp_g=cp_g, cp_l=cp_l,
                  h_fn=lambda D: 2.0 * k_g / D, T_g=330.0,
                  Y_surface_fn=Y_s, Y_gas=0.0, Q_evap_fn=Q_evap)
    return integrate_droplet(D0=80e-6, T_d0=293.0, t_end=t_end, **kwargs), kwargs


def test_coupled_case_matches_independent_tight_reference():
    result, kwargs = coupled_case(t_end=0.05)

    # Independent reference: Radau at rtol 1e-13 on the same (q, T) system.
    from droplet_journey.model import d2_law_constant as K_of

    def rhs(t, y):
        q, T = y
        D = np.sqrt(max(q, 0.0))
        B_M = spalding_mass_number(kwargs["Y_surface_fn"](T), 0.0)
        dq = -K_of(kwargs["rho_l"], kwargs["k_g"], kwargs["cp_g"], B_M)
        m = kwargs["rho_l"] * np.pi * D**3 / 6.0
        dT = lumped_heating_rhs(T, D, kwargs["T_g"], kwargs["h_fn"](D), m,
                                kwargs["cp_l"], kwargs["Q_evap_fn"](T, D))
        return [dq, dT]

    ref = solve_ivp(rhs, (0.0, 0.05), [result["D"][0] ** 2, result["T_d"][0]],
                    method="Radau", rtol=1e-13, atol=1e-18, dense_output=True)
    assert ref.success
    t_common = result["t"][result["t"] <= ref.t[-1]]
    q_ref, T_ref = ref.sol(t_common)
    np.testing.assert_allclose(result["D"][result["t"] <= ref.t[-1]] ** 2, q_ref,
                               rtol=1e-7, atol=1e-18)
    np.testing.assert_allclose(result["T_d"][result["t"] <= ref.t[-1]], T_ref,
                               rtol=1e-7, atol=1e-9)


def test_invalid_inputs_raise_meaningful_errors():
    with pytest.raises(ValueError, match="D0"):
        integrate_droplet(D0=-1e-4, T_d0=300.0, t_end=1.0, rho_l=997.0,
                          k_g=0.026, cp_g=1010.0, cp_l=4180.0,
                          h_fn=lambda D: 0.0, T_g=330.0,
                          Y_surface_fn=lambda T: 0.1, Y_gas=0.0,
                          Q_evap_fn=lambda T, D: 0.0)
    with pytest.raises(ValueError, match="t_end"):
        integrate_droplet(D0=1e-4, T_d0=300.0, t_end=float("nan"), rho_l=997.0,
                          k_g=0.026, cp_g=1010.0, cp_l=4180.0,
                          h_fn=lambda D: 0.0, T_g=330.0,
                          Y_surface_fn=lambda T: 0.1, Y_gas=0.0,
                          Q_evap_fn=lambda T, D: 0.0)
    with pytest.raises(ValueError, match="epsilon"):
        frozen_run(80e-6, 0.1, epsilon=2.0)


def test_droplet_rhs_zero_and_nonfinite_guard():
    dD, dT = droplet_rhs(0.0, 300.0, 997.0, 0.026, 1010.0, 4180.0,
                         lambda D: 0.0, 330.0, lambda T: 0.1, 0.0,
                         lambda T, D: 0.0)
    assert dD == 0.0 and dT == 0.0
    with pytest.raises(FloatingPointError, match="non-finite"):
        droplet_rhs(50e-6, 300.0, 997.0, 0.026, 1010.0, 4180.0,
                    lambda D: 0.0, 330.0, lambda T: 0.1, 0.0,
                    lambda T, D: float("nan"))
