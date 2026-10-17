"""Produce all verification results, figures, animation, and raw data.

Runs from the repository root:

    python examples/run_verification.py

Cases:
1. Frozen-surface validation: constant B_M, no heating — the coupled
   integrator must reproduce the exact constant-K d^2 law.
2. Coupled heating case: water droplet in dry air, Clapeyron surface
   (boiling regime rejected), Nu = 2 conduction closure h(D) = 2 k_g / D,
   latent sink consistent with the d^2 law. Stopped at the documented
   positive threshold q = 1e-3 q0 with remaining mass and reason recorded.

Outputs under ``docs/results/``: figures (white background, no grid lines),
a 41-frame animation sampled from the DENSE solution, CSVs, run_summary.json.
"""

from __future__ import annotations

import json
import pathlib
import platform
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.animation as mpl_animation
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from droplet_journey import (
    d2_law_constant,
    d2_law_diameter,
    extinction_time,
    integrate_droplet,
    saturation_vapor_mass_fraction,
    spalding_mass_number,
)

RESULTS = pathlib.Path(__file__).resolve().parents[1] / "docs" / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

BLUE, RED, GREEN, GRAY = "#1f4e79", "#b23a48", "#0b6e4f", "#555555"
plt.rcParams.update(
    {
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.size": 10,
    }
)

# Water-in-air reference properties (SI); constant-property assumptions.
RHO_L, CP_L, CP_G, K_G, H_FG, R_VAP = 997.0, 4180.0, 1010.0, 0.026, 2.26e6, 461.5
P_REF, T_REF, P_AMB, MW_RATIO = 3533.0, 300.0, 101325.0, 0.622
EPSILON = 1.0e-3  # documented threshold: stop at q = 1e-3 q0


def surface_mass_fraction(T_d):
    """Clapeyron surface vapor mass fraction; boiling regime REJECTED."""
    return saturation_vapor_mass_fraction(T_d, P_REF, T_REF, H_FG, R_VAP,
                                          P_AMB, MW_RATIO)


def mdot(T_d, D):
    """Evaporation rate [kg/s] consistent with the d^2 law at (T_d, D)."""
    B_M = spalding_mass_number(surface_mass_fraction(T_d), 0.0)
    K = d2_law_constant(RHO_L, K_G, CP_G, B_M)
    return RHO_L * np.pi * D * K / 4.0


def case_validation(D0=80e-6):
    B_M = spalding_mass_number(0.2, 0.0)
    K = d2_law_constant(RHO_L, K_G, CP_G, B_M)
    t_e = extinction_time(D0, K)
    result = integrate_droplet(
        D0=D0, T_d0=300.0, t_end=0.7 * t_e, rho_l=RHO_L, k_g=K_G, cp_g=CP_G,
        cp_l=CP_L, h_fn=lambda D: 0.0, T_g=300.0,
        Y_surface_fn=lambda T: 0.2, Y_gas=0.0, Q_evap_fn=lambda T, D: 0.0,
        epsilon=EPSILON,
    )
    exact = d2_law_diameter(result["t"], D0, K)
    rel_err = np.abs(result["D"] - exact) / D0
    rows = ["time_s,D_numeric_m,D_exact_m,relative_error"]
    for t, d, de, e in zip(result["t"], result["D"], exact, rel_err):
        rows.append(f"{t:.9e},{d:.12e},{de:.12e},{e:.12e}")
    (RESULTS / "validation_d2.csv").write_text("\n".join(rows), encoding="utf-8")
    return result, exact, rel_err, t_e


def case_coupled(D0=80e-6, T_d0=293.0, T_g=330.0, t_end=30.0, n_frames=41):
    result = integrate_droplet(
        D0=D0, T_d0=T_d0, t_end=t_end, rho_l=RHO_L, k_g=K_G, cp_g=CP_G,
        cp_l=CP_L, h_fn=lambda D: 2.0 * K_G / D, T_g=T_g,
        Y_surface_fn=surface_mass_fraction, Y_gas=0.0,
        Q_evap_fn=lambda T, D: mdot(T, D) * H_FG,
        epsilon=EPSILON,
    )
    # Dense uniform sampling for the animation (uses sol.sol, not the
    # adaptive nodes).
    t_grid = np.linspace(0.0, result["t"][-1], n_frames)
    q_grid, T_grid = result["dense"](t_grid)
    D_grid = np.sqrt(np.maximum(q_grid, 0.0))
    rows = ["time_s,D_m,T_d_K"]
    for t, d, td in zip(result["t"], result["D"], result["T_d"]):
        rows.append(f"{t:.9e},{d:.12e},{td:.9f}")
    (RESULTS / "coupled_history.csv").write_text("\n".join(rows), encoding="utf-8")
    return result, (t_grid, D_grid, T_grid)


def make_cover(val, exact, rel_err):
    fig, ax = plt.subplots(1, 3, figsize=(13.0, 4.0))
    ax[0].plot(val["t"] * 1e3, val["D"] * 1e6, "-", color=RED, lw=1.8, label="coupled integrator")
    ax[0].plot(val["t"] * 1e3, exact * 1e6, "--", color=GRAY, lw=1.6, label="exact d² law")
    ax[0].set_xlabel("t [ms]")
    ax[0].set_ylabel("D [µm]")
    ax[0].set_title("Diameter: frozen-surface validation")
    ax[0].legend(frameon=False, fontsize=8)

    ax[1].plot(val["t"] * 1e3, (val["D"] ** 2) * 1e12, "-", color=BLUE, lw=1.8)
    ax[1].plot(val["t"] * 1e3, (exact**2) * 1e12, "--", color=GRAY, lw=1.6)
    ax[1].set_xlabel("t [ms]")
    ax[1].set_ylabel("D² [µm²]")
    ax[1].set_title("D² is linear in time (d² law)")

    ax[2].semilogy(val["t"] * 1e3, rel_err, color=GREEN, lw=1.8)
    ax[2].set_xlabel("t [ms]")
    ax[2].set_ylabel("|D_num − D_exact| / D0 [-]")
    ax[2].set_title(f"Integration error (max {rel_err.max():.1e})")
    fig.suptitle("Droplet Journey — frozen-surface validation against the exact d² law", fontsize=13)
    fig.tight_layout()
    fig.savefig(RESULTS / "validation_d2.png", dpi=150, bbox_inches="tight")
    fig.savefig(RESULTS / "validation_d2.svg", dpi=180, bbox_inches="tight")
    plt.close(fig)


def make_coupled_figure(coupled):
    fig, ax = plt.subplots(1, 2, figsize=(10.4, 4.0))
    ax[0].plot(coupled["t"], coupled["D"] * 1e6, color=BLUE, lw=1.8)
    ax[0].set_xlabel("t [s]")
    ax[0].set_ylabel("D [µm]")
    ax[0].set_title("Coupled case: evaporating diameter\n(stopped at documented threshold)")
    ax[1].plot(coupled["t"], coupled["T_d"], color=RED, lw=1.8)
    ax[1].axhline(330.0, color=GRAY, ls="--", lw=1.0, label="gas temperature")
    ax[1].set_xlabel("t [s]")
    ax[1].set_ylabel("T_d [K]")
    ax[1].set_title("Droplet temperature (heating + evaporation cooling)")
    ax[1].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(RESULTS / "coupled_case.png", dpi=150, bbox_inches="tight")
    fig.savefig(RESULTS / "coupled_case.svg", dpi=180, bbox_inches="tight")
    plt.close(fig)


def make_gif(grids):
    t_grid, D_grid, T_grid = grids
    fig, ax = plt.subplots(2, 1, figsize=(8.0, 7.0), dpi=140)
    (line_D,) = ax[0].plot([], [], color=BLUE, lw=1.8, marker="o", markevery=[-1], ms=4)
    ax[0].set_xlim(0.0, t_grid[-1])
    ax[0].set_ylim(0.0, D_grid[0] * 1e6 * 1.05)
    ax[0].set_xlabel("t [s]")
    ax[0].set_ylabel("D [µm]")
    title_D = ax[0].set_title("Diameter at t = 0.00 s")
    (line_T,) = ax[1].plot([], [], color=RED, lw=1.8, marker="o", markevery=[-1], ms=4)
    ax[1].set_xlim(0.0, t_grid[-1])
    margin = max(0.2, 0.12 * float(np.ptp(T_grid)))
    ax[1].set_ylim(T_grid.min() - margin, T_grid.max() + margin)
    ax[1].text(0.98, 0.08, "Gas temperature: 330 K; droplet close-up",
               transform=ax[1].transAxes, ha="right", fontsize=9)
    ax[1].set_xlabel("t [s]")
    ax[1].set_ylabel("T_d [K]")
    title_T = ax[1].set_title("Droplet temperature at t = 0.00 s")
    fig.tight_layout(pad=2.0)

    writer = mpl_animation.PillowWriter(fps=8)
    with writer.saving(fig, RESULTS / "droplet_journey.gif", dpi=140):
        for i, t in enumerate(t_grid):
            line_D.set_data(t_grid[: i + 1], D_grid[: i + 1] * 1e6)
            line_T.set_data(t_grid[: i + 1], T_grid[: i + 1])
            title_D.set_text(f"Diameter at t = {t:.2f} s")
            title_T.set_text(f"Droplet temperature at t = {t:.2f} s")
            writer.grab_frame()
    plt.close(fig)
    return len(t_grid)


def main():
    val, exact, rel_err, t_e_val = case_validation()
    make_cover(val, exact, rel_err)

    coupled, grids = case_coupled()
    make_coupled_figure(coupled)
    frames = make_gif(grids)

    summary = {
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": __import__("scipy").__version__,
            "system": platform.platform(),
        },
        "validation": {
            "case": "frozen surface, constant B_M, no heating",
            "max_relative_diameter_error": float(rel_err.max()),
            "analytical_model_extinction_time_s": t_e_val,
        },
        "coupled": {
            "case": "water in dry air, Clapeyron surface (boiling rejected), "
                    "Nu=2 closure h(D)=2 k_g/D, Le=1 mass closure",
            "D0_um": 80.0,
            "stop_reason": coupled["stop_reason"],
            "threshold_stopping_time_s": coupled["stop_time"],
            "remaining_mass_fraction": coupled["remaining_mass_fraction"],
            "T_d_final_K": float(coupled["T_d"][-1]),
            "D_final_um": float(coupled["D"][-1] * 1e6),
            "solver_success": coupled["solver_success"],
            "note": "threshold stop near the singular endpoint of the 0-D "
                    "closure; a model result under stated constant "
                    "conditions, not a validated lifetime",
        },
        "animation": {"file": "droplet_journey.gif", "frames": frames},
    }
    (RESULTS / "run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
