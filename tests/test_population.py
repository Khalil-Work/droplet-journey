"""Independent hand values, scaling and domain checks for the population model."""
import numpy as np
import pytest

from droplet_journey.population import evaporate_population


def test_single_group_hand_solution_and_exact_extinction():
    d0, k = 100e-6, 2.5e-8
    lifetime = d0*d0/k
    result = evaporate_population([d0], [8], [0, .75*lifetime, lifetime], k_m2_s=k)
    assert result.diameter_m[1, 0] == pytest.approx(d0/2)
    assert result.liquid_mass_kg[1] == pytest.approx(result.initial_mass_kg/8)
    assert result.surviving_number.tolist() == [8, 8, 0]
    assert result.liquid_mass_kg[-1] == 0
    assert np.isnan(result.d32_m[-1])


def test_weighted_moments_mass_ledger_and_monotonicity():
    result = evaporate_population([20e-6, 40e-6], [2, 1], np.linspace(0, .08, 51), k_m2_s=2.5e-8)
    assert result.d32_m[0] == pytest.approx((2*20**3+40**3)/(2*20**2+40**2)*1e-6)
    assert np.all(np.diff(result.liquid_mass_kg) <= 0)
    assert np.all(np.diff(result.vapor_mass_kg) >= 0)
    assert np.all(np.diff(result.surviving_number) <= 0)
    assert result.mass_balance_relative_error < 3e-16


def test_diameter_time_and_mass_scaling():
    d0 = np.array([30e-6, 70e-6])
    time = np.linspace(0, .15, 21)
    original = evaporate_population(d0, [3, 2], time, k_m2_s=2e-8)
    scaled = evaporate_population(2*d0, [3, 2], 4*time, k_m2_s=2e-8)
    np.testing.assert_allclose(scaled.diameter_m, 2*original.diameter_m)
    np.testing.assert_allclose(scaled.liquid_mass_kg, 8*original.liquid_mass_kg)
    np.testing.assert_array_equal(scaled.surviving_number, original.surviving_number)


def test_stationary_case_and_late_first_sample():
    stationary = evaporate_population([30e-6], [.5], [2, 3], k_m2_s=0)
    np.testing.assert_array_equal(stationary.vapor_mass_kg, 0)
    late = evaporate_population([100e-6], [1], [.3], k_m2_s=2.5e-8)
    assert late.vapor_mass_kg[0]/late.initial_mass_kg == pytest.approx(7/8)


@pytest.mark.parametrize("change", [
    dict(diameter0_m=[0]), dict(diameter0_m=[np.nan]),
    dict(number_weights=[-1]), dict(number_weights=[0]),
    dict(number_weights=[1, 2]), dict(time_s=[1, 0]),
    dict(time_s=[0, 0]), dict(time_s=[np.inf]),
    dict(k_m2_s=-1), dict(k_m2_s=np.nan), dict(density_kg_m3=0),
])
def test_invalid_population_inputs(change):
    arguments = dict(diameter0_m=[30e-6], number_weights=[1], time_s=[0, .1], k_m2_s=2e-8)
    arguments.update(change)
    with pytest.raises(ValueError):
        evaporate_population(**arguments)
