# Verification evidence

Reproduced by `python -m pytest tests -q` (10 tests) and
`python examples/run_verification.py` (figures + GIF + CSV/JSON under
`docs/results/`). Verified on Windows 11, Python 3.13.5, NumPy 2.1.3,
SciPy 1.15.3.

## What each test establishes

| Test | Property established |
| --- | --- |
| `test_spalding_number_values_and_limits` | B_M values and rejection of invalid mass fractions. |
| `test_d2_law_constant_matches_hand_computation` | K against an independent hand computation (ln(1+B)=1 exactly); input rejection including NaN. |
| `test_exact_constant_k_solution_and_extinction` | The exact constant-K law is linear in D² and terminates at t_e = D0²/K. |
| `test_clapeyron_anchor_and_boiling_rejection` | The fit reproduces its anchor; the boiling regime RAISES rather than being clamped. |
| `test_sphere_area_scaling_of_initial_heating_rate` | Physical area/mass scaling: at fixed h, doubling D halves the initial heating rate (area ~ D², mass ~ D³) — this is the test that catches a wrong heat-transfer area. |
| `test_coupled_integrator_reproduces_frozen_d2_law` | With frozen surface state and no heating, the coupled integrator matches the exact d² law to ~1e-16 relative. |
| `test_threshold_stopping_reports_mass_and_reason` | The coupled stop is at the documented positive threshold with remaining mass ε^1.5 and reason recorded; D never silently reaches zero. |
| `test_coupled_case_matches_independent_tight_reference` | The production DOP853 solution agrees with an independent Radau (rtol 1e-13) integration of the same coupled system on the shared valid interval (1e-7 relative). |
| `test_invalid_inputs_raise_meaningful_errors` | Finite-positive validation of D0, t_end, epsilon. |
| `test_droplet_rhs_zero_and_nonfinite_guard` | Zero state semantics and rejection of non-finite callback outputs. |

## Recorded run (`docs/results/`)

- `validation_d2.png/.csv`: frozen-surface validation; max relative
  diameter error 1.7e-16 (roundoff); analytical model extinction time
  0.1389 s for the 80 µm case.
- `coupled_case.png`, `coupled_history.csv`: 80 µm water droplet in dry
  330 K air with the Nu = 2 closure h(D) = 2k_g/D and Clapeyron surface.
  Threshold stop (ε = 1e-3) at t = 1.937 s with 3.16e-5 of the mass
  remaining and T_d = 294.2 K — a mild wet-bulb-like equilibrium between
  convective heating and evaporation cooling. This is a model result under
  the stated constant conditions, not a validated lifetime.
- `droplet_journey.gif`: 41 frames sampled uniformly from the DENSE solver
  solution (not the adaptive nodes).
- `run_summary.json`: the metrics above plus runtime versions and the
  stop record.

## Limits of this evidence

- 0-D, single droplet, quiescent constant gas, no relative convection, no
  radiation, no multicomponent liquid, no droplet-internal gradients
  (lumped capacitance).
- The mass closure assumes Le = 1; the thermal closure neglects the
  departing vapor's sensible enthalpy; properties are constant; the
  Clapeyron fit uses one anchor and constant h_vap.
- Boiling states are rejected, not modeled; high-pressure or
  superheated-regime behavior is outside scope.
- The threshold stop is a property of the closure's singular endgame, not
  of the physics; the recorded stopping time must not be read as a
  complete-evaporation lifetime.
- No comparison with experimental data is claimed.
