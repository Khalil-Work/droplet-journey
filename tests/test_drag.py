"""Tests for droplet_journey.drag (Schiller-Naumann drag, terminal velocity)."""

import pytest

from droplet_journey.drag import (
    drag_coefficient,
    reynolds,
    terminal_velocity,
)

RHO_AIR, MU_AIR, RHO_WATER = 1.2, 1.8e-5, 998.0


def test_drag_coefficient_limits():
    # Stokes limit: C_d ~ 24/Re as Re -> 0.
    assert drag_coefficient(1e-8) == pytest.approx(24.0 / 1e-8, rel=1e-6)
    # Published Schiller-Naumann value at Re = 100: 24/100*(1+0.15*100^0.687).
    assert drag_coefficient(100.0) == pytest.approx(
        0.24 * (1.0 + 0.15 * 100.0 ** 0.687), rel=1e-12)
    with pytest.raises(ValueError):
        drag_coefficient(0.0)
    with pytest.raises(ValueError):
        drag_coefficient(-10.0)


def test_terminal_velocity_satisfies_force_balance():
    v = terminal_velocity(1.0e-3, RHO_WATER, RHO_AIR, MU_AIR)
    # Hand check: 1 mm water droplet falls at ~4 m/s (textbook ballpark).
    assert 3.0 < v < 5.0
    # Verify the balance: weight-buoyancy == drag to roundoff.
    d = 1.0e-3
    re = reynolds(v, d, RHO_AIR, MU_AIR)
    drag = 0.5 * RHO_AIR * drag_coefficient(re) * (3.141592653589793 / 4 * d * d) * v * v
    weight = (RHO_WATER - RHO_AIR) * 9.81 * (3.141592653589793 / 6 * d ** 3)
    assert drag == pytest.approx(weight, rel=1e-9)


def test_terminal_velocity_scales_as_diameter():
    # Stokes regime (tiny droplets): v ~ D^2. Compare 2 and 4 micrometre.
    v2 = terminal_velocity(2.0e-6, RHO_WATER, RHO_AIR, MU_AIR)
    v4 = terminal_velocity(4.0e-6, RHO_WATER, RHO_AIR, MU_AIR)
    assert v4 / v2 == pytest.approx(4.0, rel=0.02)
    assert reynolds(v2, 2.0e-6, RHO_AIR, MU_AIR) < 1.0  # genuinely Stokes


def test_terminal_velocity_grows_monotonically_with_diameter():
    sizes = (1e-4, 3e-4, 1e-3, 3e-3)
    vs = [terminal_velocity(d, RHO_WATER, RHO_AIR, MU_AIR) for d in sizes]
    assert all(b > a for a, b in zip(vs, vs[1:]))


def test_drag_input_validation():
    with pytest.raises(ValueError):
        terminal_velocity(0.0, RHO_WATER, RHO_AIR, MU_AIR)
    with pytest.raises(ValueError):
        terminal_velocity(1e-3, 0.0, RHO_AIR, MU_AIR)
    with pytest.raises(ValueError):
        terminal_velocity(1e-3, RHO_WATER, RHO_AIR, -1e-5)
    with pytest.raises(ValueError):
        reynolds(float("nan"), 1e-3, RHO_AIR, MU_AIR)
