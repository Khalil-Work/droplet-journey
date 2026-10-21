# 💠 Khalil Droplet Population Atlas

A bespoke addition to Droplet Journey for KHALIL: an analytic,
weighted population experiment, with a liquid/vapor ledger and an animated
view of each size group's disappearance. This addition was written for this
repository; it contains no imported third-party solver implementation.

▶ **[Open the animated atlas](https://github.com/khalil-work/droplet-journey/raw/refs/heads/main/docs/results/population_atlas/khalil_population_atlas.gif)**

![Computed weighted droplet population and mass transfer](results/population_atlas/khalil_population_atlas.gif)

## 🧮 Equations and units

For a non-interacting group i, its constant-condition model is

```text
D_i(t) = sqrt(max(D_i0² - K t, 0))       D [m], K [m²/s], t [s]
t_ext,i = D_i0² / K                     K > 0
m_liq(t) = (rho*pi/6) sum_i n_i D_i(t)³  rho [kg/m³], m [kg]
m_vap(t) = m_liq(0) - m_liq(t)
N_surviving(t) = sum_i n_i 1[D_i(t) > 0]
D32(t) = sum_i n_i D_i(t)³ / sum_i n_i D_i(t)²
```

The weights n_i are number multiplicities, not mass fractions; fractional
quadrature weights are allowed. The model has exact, finite extinction at
the constant-K lifetime. This differs from the **positive stopping threshold
of the thermally coupled single-droplet model** elsewhere in the project.
K=0 is a stationary population. When every droplet has disappeared, D32 is
undefined: the API returns NaN, the CSV leaves that field blank, and the
plot ends the size curve rather than drawing a fictitious zero size.

The ledger verifies arithmetic bookkeeping; it does not independently
validate the evaporation physics. Individual diameters use an analytic
solution, so this experiment has no time-stepping error or convergence claim.

## 🔬 The comparison

The documented population has 16 logarithmically spaced initial diameters
from 20 to 140 µm and deterministic lognormal-shaped number weights. K is
2.5e-8 m²/s and rho is 997 kg/m³. These are illustrative constant inputs,
not experimentally fitted operating conditions. Exact values are saved in
[`case.json`](results/population_atlas/case.json).

The dashed comparison has one initial diameter equal to the population's
initial D32. Its number weight is adjusted to match the same **initial total
liquid mass**. It consequently also matches initial total surface area.
It does not match every size moment or the original number count. The
different histories demonstrate what that single-size approximation loses.

Circle diameters use one consistent scale. Their positions are schematic;
the picture does not simulate trajectories. Each circle represents a size
group, not one particle or an equal mass. Dotted outlines are initial-size
references. The right panels show the actual computed mass and size histories.

## ⚡ Reproduce and use

```bash
python -m pip install -r requirements-dev.txt
python -m pytest tests -q
python examples/khalil_population_atlas.py
```

```python
from droplet_journey.population import evaporate_population

history = evaporate_population(
    [20e-6, 60e-6, 100e-6], [200, 100, 20], [0, .1, .2, .4],
    k_m2_s=2.5e-8, density_kg_m3=997.0,
)
print(history.liquid_mass_kg)
print(history.mass_balance_relative_error)
```

The returned `PopulationHistory` contains sampled diameters `(time, group)`,
liquid/vapor masses, surviving number weight, D32 and initial mass. Inputs
must be finite 1-D arrays of matching sizes; diameters and density must be
positive, weights and times nonnegative, time strictly increasing, K>=0.
At least one number weight must be positive. Times need not start at zero;
vapor mass still refers to the original initial population.

The generator writes a 61-frame GIF, a comparison PNG, 61 sampled mass-history
rows, 976 sampled group-history rows (plus CSV headers), case inputs and
runtime/metric JSON. All files live
under [`docs/results/population_atlas/`](results/population_atlas/).

## ✅ Verification and limits

The added tests cover a hand-computed half-diameter/eighth-mass state,
exact extinction, weighted D32, monotonic histories and ledger closure,
diameter/time/mass scaling, K=0, late first sampling, and invalid inputs.
Existing thermal-model tests remain part of the suite.

This is a constant-K population benchmark. It includes no gas feedback,
droplet interaction, heating, breakup, coalescence, particle motion,
spatial spray field or experimental validation. The source is provided
under this repository's MIT license; standard mathematical relations are
not a claim of a newly discovered physical law.
