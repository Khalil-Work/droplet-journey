"""Reproduce the bespoke population atlas, raw histories and animated ledger.

Run from any directory: python examples/khalil_population_atlas.py
All liquid sizes are analytic constant-K predictions; glyph locations are
schematic, not simulated particle positions.
"""
import csv
import json
from pathlib import Path
import platform
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import PillowWriter
from matplotlib.patches import Circle
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from droplet_journey.population import evaporate_population

BLUE, ORANGE, INK = "#087fbc", "#d75d23", "#172e46"


def main():
    out = ROOT / "docs/results/population_atlas"
    out.mkdir(parents=True, exist_ok=True)
    diameters = np.geomspace(20, 140, 16) * 1e-6
    weights = np.maximum(1, np.rint(1000*np.exp(-.5*(np.log(diameters/70e-6)/.36)**2)))
    k, density = 2.5e-8, 997.0
    times = np.linspace(0, diameters[-1]**2/k, 61)
    history = evaporate_population(diameters, weights, times, k_m2_s=k, density_kg_m3=density)
    d32_initial = history.d32_m[0]
    mono_weight = history.initial_mass_kg / (density*np.pi/6*d32_initial**3)
    mono = evaporate_population([d32_initial], [mono_weight], times, k_m2_s=k, density_kg_m3=density)
    with (out/"mass_history.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(["time_s", "liquid_kg", "vapor_kg", "surviving_number_weight", "d32_um", "equal_mass_mono_liquid_kg"])
        for i, t in enumerate(times):
            writer.writerow([t, history.liquid_mass_kg[i], history.vapor_mass_kg[i], history.surviving_number[i],
                             history.d32_m[i]*1e6 if np.isfinite(history.d32_m[i]) else "", mono.liquid_mass_kg[i]])
    with (out/"group_history.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(["time_s", "group", "initial_diameter_um", "diameter_um", "number_weight"])
        for i, t in enumerate(times):
            for j in range(diameters.size):
                writer.writerow([t, j+1, diameters[j]*1e6, history.diameter_m[i, j]*1e6, weights[j]])
    case = dict(name="Khalil droplet population atlas", k_m2_s=k, density_kg_m3=density,
                diameter0_m=diameters.tolist(), number_weights=weights.tolist(),
                comparison="Same initial liquid mass and D32; monodisperse number weight adjusted accordingly",
                scope="Independent analytic constant-K droplets; no thermal coupling, gas feedback, motion, breakup or experiments")
    (out/"case.json").write_text(json.dumps(case, indent=2), encoding="utf-8")
    plt.rcParams.update({"figure.facecolor":"white", "axes.facecolor":"white", "savefig.facecolor":"white",
                         "axes.grid":False, "font.size":10, "text.color":INK, "axes.labelcolor":INK,
                         "axes.spines.top":False, "axes.spines.right":False})
    fig = plt.figure(figsize=(11.2, 7.6), layout="constrained")
    grid = fig.add_gridspec(2, 2, width_ratios=[1.05, 1])
    glyphs = fig.add_subplot(grid[:, 0])
    mass = fig.add_subplot(grid[0, 1])
    size = fig.add_subplot(grid[1, 1])
    fig.suptitle("KHALIL  |  DROPLET POPULATION ATLAS", fontsize=18, weight="bold")
    positions = [(j%4, 3-j//4) for j in range(16)]
    glyphs.set(xlim=(-.6, 3.6), ylim=(-.75, 3.75), aspect="equal")
    glyphs.axis("off")
    disks=[]
    for j, (x, y) in enumerate(positions):
        radius = .39*diameters[j]/diameters[-1]
        glyphs.add_patch(Circle((x,y), radius, fill=False, edgecolor="#bbc6cf", linestyle=":", linewidth=1))
        disk = Circle((x,y), radius, color=BLUE, alpha=.88)
        glyphs.add_patch(disk); disks.append(disk)
        glyphs.text(x, y-.48, f"{diameters[j]*1e6:.0f} µm", ha="center", fontsize=8, color=INK)
    glyphs.text(.5, -.02, "Each circle represents one diameter group.\nCircle diameter follows D(t); positions are schematic.\nDotted outlines show initial diameters; labels give D₀.",
                transform=glyphs.transAxes, ha="center", va="top", fontsize=9, color=INK)
    glyph_title = glyphs.set_title("", fontsize=12, pad=16)
    fraction = 100*history.liquid_mass_kg/history.initial_mass_kg
    vapor = 100*history.vapor_mass_kg/history.initial_mass_kg
    mono_fraction = 100*mono.liquid_mass_kg/mono.initial_mass_kg
    mass.set(xlim=(0, times[-1]), ylim=(-3,103), xlabel="Time [s]", ylabel="Initial liquid mass [%]", title="Liquid → vapor: a closed mass ledger")
    line_liq, = mass.plot([], [], color=BLUE, lw=2.5, marker="o", markevery=[-1], ms=4, label="Population liquid")
    line_vap, = mass.plot([], [], color=ORANGE, lw=2, marker="o", markevery=[-1], ms=4, label="Cumulative vapor")
    line_mono, = mass.plot([], [], color="#778696", ls="--", label="Matched-D₃₂ mono liquid")
    mass.legend(frameon=False, fontsize=8, loc="center right")
    size.set(xlim=(0,times[-1]), ylim=(0,155), xlabel="Time [s]", ylabel="D₃₂ [µm]", title="Size of the surviving liquid")
    line_size, = size.plot([], [], color=BLUE, lw=2.5, marker="o", markevery=[-1], ms=4, label="Weighted population")
    line_size_mono, = size.plot([], [], color="#778696", ls="--", label="Monodisperse")
    size.legend(frameon=False, fontsize=8, loc="lower right")
    size.text(.02,.97,"No liquid left → D₃₂ undefined", transform=size.transAxes, va="top", fontsize=8, color=INK)
    def draw(i):
        for j, disk in enumerate(disks):
            disk.set_radius(.39*history.diameter_m[i,j]/diameters[-1])
        glyph_title.set_text(f"t = {times[i]:.3f} s  ·  {np.count_nonzero(history.diameter_m[i]):d}/16 groups remain\nLiquid mass = {fraction[i]:.1f}% of initial")
        stop=slice(0,i+1)
        line_liq.set_data(times[stop], fraction[stop]); line_vap.set_data(times[stop], vapor[stop])
        line_mono.set_data(times[stop], mono_fraction[stop])
        line_size.set_data(times[stop], 1e6*history.d32_m[stop]); line_size_mono.set_data(times[stop], 1e6*mono.d32_m[stop])
    writer = PillowWriter(fps=10)
    with writer.saving(fig, str(out/"khalil_population_atlas.gif"), dpi=100):
        for i in range(times.size):
            draw(i); writer.grab_frame(facecolor="white")
    draw(times.size-1)
    fig.savefig(out/"population_comparison.png", dpi=160, facecolor="white")
    plt.close(fig)
    summary = dict(python=platform.python_version(), numpy=np.__version__, matplotlib=matplotlib.__version__,
                   frames=int(times.size), initial_mass_kg=history.initial_mass_kg,
                   initial_d32_um=float(d32_initial*1e6), initial_number_weight=float(weights.sum()),
                   final_liquid_mass_kg=float(history.liquid_mass_kg[-1]),
                   maximum_relative_mass_ledger_error=history.mass_balance_relative_error,
                   time_range_s=[float(times[0]), float(times[-1])])
    (out/"run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
