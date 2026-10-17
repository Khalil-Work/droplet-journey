# Add reproducible benchmark suite, animation, and technical guides

## What this adds
- Recorded benchmark suite (examples/run_verification.py):
  - frozen-surface validation against the exact constant-K d² law (max relative diameter error 1.7e-16);
  - coupled heating case: 80 µm water droplet in dry 330 K air, Nu = 2 closure h(D) = 2k_g/D, Clapeyron surface with boiling rejection, threshold stop ε = 1e-3 at t = 1.937 s with 3.16e-5 mass remaining and T_d = 294.2 K;
  - 41-frame animation sampled uniformly from the DENSE solver solution (not adaptive nodes).
- Raw data under docs/results/: validation_d2.csv, coupled_history.csv, run_summary.json; white/no-grid figures with rasterized SVG exports.
- Technical guides: docs/model.md (equations incl. unit-Lewis-number closure statement and sphere-area πD² heating), docs/api.md, docs/verification.md (per-test evidence + limits).
- README with executed example, metrics table, animation, and honest stopping-event semantics.

## Physical corrections made during review (before any results were published)
- Heat-transfer area fixed to the sphere surface πD² (the draft used πD — a 1/D error now caught by a dedicated area-scaling test: at fixed h, doubling D halves the initial heating rate).
- The thermal endgame is handled by a documented POSITIVE threshold stop (q = ε q0) reporting stop time, remaining mass ε^1.5, and reason — replacing the earlier temperature-freeze near extinction, which is no longer needed with the corrected area and q-formulation dq/dt = −K computed directly.
- Solver success is checked; integration failure raises RuntimeError and can never appear as droplet survival.
- Boiling regime (p_sat ≥ p_ambient) is rejected with a clear error, never clamped.
- Unit-Lewis-number assumption of the K closure stated explicitly in code and docs.

## Verification
- 10/10 tests pass locally (Python 3.13.5 / NumPy 2.1.3 / SciPy 1.15.3), including an independent Radau (rtol 1e-13) reference agreeing with the production DOP853 solution to < 1e-7 relative on the shared valid interval.
- README example executed.
- CI runs the suite on this PR.
