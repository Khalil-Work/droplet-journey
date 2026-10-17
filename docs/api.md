# API reference

Module `droplet_journey.model` (re-exported from `droplet_journey`). SI
units throughout.

- `saturation_pressure_clapeyron(T_d, p_ref, T_ref, h_vap, R_vap) -> p_sat [Pa]`
- `spalding_mass_number(Y_surface, Y_gas) -> B_M` — raises `ValueError`
  outside [0,1).
- `d2_law_constant(rho_l, k_g, cp_g, B_M) -> K [m²/s]` — raises for
  non-positive inputs or B_M ≤ 0.
- `d2_law_diameter(t, D0, K) -> D(t) [m]` — exact constant-K solution,
  clipped at 0 for t ≥ t_e.
- `extinction_time(D0, K) -> t_e [s]`.
- `lumped_heating_rhs(T_d, D, T_g, h, m_l, cp_l, Q_evap) -> dT_d/dt [K/s]`.
- `droplet_rhs(D, T_d, rho_l, k_g, cp_g, cp_l, h, T_g, Y_surface_fn, Y_gas, Q_evap_fn) -> (dD/dt, dT_d/dt)`
  — returns (0, 0) at D ≤ 0.
- `integrate_droplet(D0, T_d0, t_end, rho_l, k_g, cp_g, cp_l, h, T_g, Y_surface_fn, Y_gas, Q_evap_fn, max_step=None) -> dict`
  with keys `t`, `D`, `T_d`, `t_extinction` (or None), `status`
  (scipy: 0 = reached t_end, 1 = terminal event, −1 = step failure).
  Uses the (D², T_d) state internally and the documented isothermal
  endgame below D0/100.

Callables: `Y_surface_fn(T_d) -> Y_s`, `Q_evap_fn(T_d, D) -> [W]`.

## Example

```python
from droplet_journey import d2_law_constant, d2_law_diameter, spalding_mass_number

B_M = spalding_mass_number(0.2, 0.0)
K = d2_law_constant(rho_l=997.0, k_g=0.026, cp_g=1010.0, B_M=B_M)
D_of_t = d2_law_diameter(t=0.05, D0=80e-6, K=K)
```
