"""Gas-phase drag and free-fall terminal velocity for a spherical droplet.

Drag coefficient: Schiller-Naumann (1935), valid for a smooth sphere at
Re_p <~ 1000 (droplets in air without deformation/oscillation):

    C_d = 24/Re_p (1 + 0.15 Re_p^0.687),   Re_p > 0
    C_d -> 24/Re_p  (Stokes limit, Re_p -> 0)

Free-fall force balance (weight - buoyancy = drag) gives the terminal
velocity implicitly; we solve it by direct fixed-point iteration on v in
the drag law, which converges monotonically here because the RHS
v -> sqrt(2 m g_eff / (rho_g C_d(Re(v)) A)) is contracting in v.

Symbols: SI. Re_p = rho_g v D / mu_g; A = pi D^2/4 projected area.
"""

from __future__ import annotations

import math
from numbers import Real


def _positive(name, value):
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite positive scalar")
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite positive, got {value}")
    return value


def drag_coefficient(re_p):
    """Schiller-Naumann sphere drag coefficient at particle Reynolds Re_p."""
    re = _positive("Re_p", re_p)
    return 24.0 / re * (1.0 + 0.15 * re ** 0.687)


def terminal_velocity(diameter, rho_d, rho_g, mu_g, g=9.81, tol=1e-10, max_iter=200):
    """Terminal free-fall velocity [m/s] of a spherical droplet.

    Solved by fixed-point iteration on the force balance
        (rho_d - rho_g) V g = 0.5 rho_g C_d(Re) A v^2
    starting from the Stokes estimate. Raises FloatingPointError if the
    iteration fails to converge within max_iter (should not happen for
    physical inputs: the map is contracting).
    """
    d = _positive("diameter", diameter)
    rho_l = _positive("rho_d", rho_d)
    rho_fl = _positive("rho_g", rho_g)
    mu = _positive("mu_g", mu_g)
    gg = _positive("g", g)
    volume = math.pi / 6.0 * d ** 3
    area = math.pi / 4.0 * d ** 2
    # Stokes start: v = g (rho_d - rho_g) D^2 / (18 mu_g).
    v = gg * (rho_l - rho_fl) * d * d / (18.0 * mu)
    for _ in range(max_iter):
        re = rho_fl * v * d / mu
        cd = drag_coefficient(max(re, 1e-12))
        v_new = math.sqrt(
            2.0 * gg * (rho_l - rho_fl) * volume / (rho_fl * cd * area)
        )
        if abs(v_new - v) <= tol * max(1.0, v_new):
            return v_new
        v = v_new
    raise FloatingPointError("terminal-velocity iteration did not converge")


def reynolds(v, diameter, rho_g, mu_g):
    """Particle Reynolds number rho_g v D / mu_g [-]."""
    vv = _positive("v", v)
    d = _positive("diameter", diameter)
    rho_fl = _positive("rho_g", rho_g)
    mu = _positive("mu_g", mu_g)
    return rho_fl * vv * d / mu
