"""Weighted, non-interacting droplet populations under a constant d-squared law.

The analytic model is D_i(t)^2 = max(D_i(0)^2 - K*t, 0).  Multiplicities
are number weights, not mass fractions. No thermal, gas-feedback or breakup
model is implied. SI units are used throughout.
"""
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class PopulationHistory:
    """Histories on the supplied time grid; diameter has shape (time, group).

    D32 is the ratio of weighted third/second diameter moments. It is NaN
    after every group has disappeared, when no liquid size is defined.
    """
    time_s: np.ndarray
    diameter_m: np.ndarray
    liquid_mass_kg: np.ndarray
    vapor_mass_kg: np.ndarray
    surviving_number: np.ndarray
    d32_m: np.ndarray
    initial_mass_kg: float

    @property
    def mass_balance_relative_error(self):
        total = self.liquid_mass_kg + self.vapor_mass_kg
        return float(np.max(np.abs(total - self.initial_mass_kg)) / self.initial_mass_kg)


def _vector(values, name, *, positive=False, nonnegative=False):
    array = np.asarray(values, dtype=float)
    if array.ndim != 1 or array.size == 0 or not np.isfinite(array).all():
        raise ValueError(f"{name} must be a nonempty finite 1-D array")
    if positive and np.any(array <= 0):
        raise ValueError(f"{name} must be positive")
    if nonnegative and np.any(array < 0):
        raise ValueError(f"{name} must be nonnegative")
    return array


def evaporate_population(diameter0_m, number_weights, time_s, *,
                         k_m2_s, density_kg_m3=997.0):
    """Return an exact constant-K population history without time stepping.

    Parameters
    ----------
    diameter0_m : 1-D array
        Strictly positive initial diameters [m].
    number_weights : 1-D array, same length
        Nonnegative number multiplicities (fractional quadrature weights
        allowed), with at least one positive entry.
    time_s : 1-D array
        Nonnegative, strictly increasing sample times [s]; need not start0.
    k_m2_s : float
        Finite nonnegative evaporation constant [m2/s]; K=0 is stationary.
    density_kg_m3 : float
        Finite positive liquid density [kg/m3].

    A group is extinct at t >= D0^2/K, including the exact boundary. Vapor
    is cumulative evaporated mass from the INITIAL population, even when
    the first requested sample is later than0. Input arrays are not mutated.
    """
    d0 = _vector(diameter0_m, "diameter0_m", positive=True)
    weights = _vector(number_weights, "number_weights", nonnegative=True)
    times = _vector(time_s, "time_s", nonnegative=True)
    if weights.shape != d0.shape or not np.any(weights > 0):
        raise ValueError("number_weights must match diameters and have positive total")
    if times.size > 1 and np.any(np.diff(times) <= 0):
        raise ValueError("time_s must be strictly increasing")
    if not np.isscalar(k_m2_s) or not np.isfinite(k_m2_s) or k_m2_s < 0:
        raise ValueError("k_m2_s must be a finite nonnegative scalar")
    if not np.isscalar(density_kg_m3) or not np.isfinite(density_kg_m3) or density_kg_m3 <= 0:
        raise ValueError("density_kg_m3 must be a finite positive scalar")
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        q0 = d0 ** 2
        initial_mass = float(density_kg_m3 * np.pi / 6 * np.sum(weights * d0 ** 3))
    if not np.isfinite(q0).all() or np.any(q0 <= 0) or not np.isfinite(initial_mass) or initial_mass <= 0:
        raise ValueError("diameter/weight scale is outside finite floating-point mass range")
    if k_m2_s == 0:
        diameter = np.broadcast_to(d0, (times.size, d0.size)).copy()
    else:
        with np.errstate(over="ignore"):
            lifetime = q0 / k_m2_s
            q = np.maximum(q0[None, :] - k_m2_s * times[:, None], 0.0)
        # Explicit lifetime comparison makes equality at extinction unambiguous.
        diameter = np.sqrt(np.where(times[:, None] < lifetime[None, :], q, 0.0))
    second = np.sum(weights * diameter ** 2, axis=1)
    third = np.sum(weights * diameter ** 3, axis=1)
    liquid = density_kg_m3 * np.pi / 6 * third
    vapor = np.maximum(initial_mass - liquid, 0.0)
    d32 = np.divide(third, second, out=np.full(times.size, np.nan), where=second > 0)
    number = np.sum(weights * (diameter > 0), axis=1)
    return PopulationHistory(times.copy(), diameter, liquid, vapor, number, d32, initial_mass)
