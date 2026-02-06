#!/usr/bin/env python3
"""
Quick plots for Emanuel RCE single-column model outputs.

Usage:
  python plot_rce_outputs.py                 # reads ./output/
  python plot_rce_outputs.py --dir output1   # reads ./output1/

Requires:
  numpy, matplotlib
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import matplotlib
import numpy as np
import matplotlib.pyplot as plt


def load_txt(path: Path) -> np.ndarray:
    try:
        return np.loadtxt(path)
    except OSError as e:
        raise SystemExit(f"Missing file: {path}") from e


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="output", help="run directory containing *.out files (default: output)")
    ap.add_argument(
        "--save",
        default="",
        help="directory to save PNGs (headless-friendly). If omitted, shows interactive windows.",
    )
    args = ap.parse_args()

    # If we're on a headless machine (no DISPLAY) and user didn't ask to show,
    # default to saving PNGs into the run directory.
    headless = (os.environ.get("DISPLAY", "") == "") and (os.environ.get("WAYLAND_DISPLAY", "") == "")
    save_dir = Path(args.save) if args.save else (Path(args.dir) if headless else None)
    if save_dir is not None:
        matplotlib.use("Agg", force=True)

    outdir = Path(args.dir)
    time = load_txt(outdir / "time.out")
    prof = load_txt(outdir / "profile.out")

    # --- time.out columns (from Users_guide.pdf) ---
    # 1 time (days)
    # 2 precip (mm/day)
    # 3 evap (mm/day)
    # 4 Ta lowest level (C)
    # 5 SST (C)
    # 8 TOA SW (W/m^2)
    # 9 TOA LW (W/m^2)
    t_days = time[:, 0]
    precip = time[:, 1]
    evap = time[:, 2]
    ta0 = time[:, 3]
    sst = time[:, 4]
    toa_sw = time[:, 7]
    toa_lw = time[:, 8]

    # --- profile.out columns (see Users_guide.pdf) ---
    # We only need the first few + last few for HW plots.
    p_hpa = prof[:, 0]
    T_c = prof[:, 1]
    q_gkg = prof[:, 2]
    rh_frac = prof[:, 4]  # appears 0-1 in your file; convert to %
    cloud_frac = prof[:, -2]
    cloud_w_gkg = prof[:, -1]

    # --- Plot: precip/evap ---
    fig1 = plt.figure()
    plt.plot(t_days, precip, label="Precip (mm/day)")
    plt.plot(t_days, evap, label="Evap (mm/day)", linestyle="--")
    plt.xlabel("Time (days)")
    plt.ylabel("mm/day")
    plt.title("Precipitation and evaporation")
    plt.grid(True, alpha=0.3)
    plt.legend()

    # --- Plot: SST and surface air temperature ---
    fig2 = plt.figure()
    plt.plot(t_days, sst, label="SST (C)")
    plt.plot(t_days, ta0, label="Ta lowest level (C)", linestyle="--")
    plt.xlabel("Time (days)")
    plt.ylabel("C")
    plt.title("Surface temperatures")
    plt.grid(True, alpha=0.3)
    plt.legend()

    # --- Plot: TOA SW/LW ---
    fig3 = plt.figure()
    plt.plot(t_days, toa_sw, label="TOA SW (W/m^2)")
    plt.plot(t_days, toa_lw, label="TOA LW (W/m^2)", linestyle="--")
    plt.xlabel("Time (days)")
    plt.ylabel("W/m^2")
    plt.title("Top-of-atmosphere fluxes")
    plt.grid(True, alpha=0.3)
    plt.legend()

    # --- Plot: equilibrium profiles ---
    fig4 = plt.figure()
    plt.plot(T_c, p_hpa, label="T (C)")
    plt.gca().invert_yaxis()
    plt.xlabel("Temperature (C)")
    plt.ylabel("Pressure (hPa)")
    plt.title("Equilibrium temperature profile")
    plt.grid(True, alpha=0.3)

    fig5 = plt.figure()
    plt.plot(rh_frac * 100.0, p_hpa, label="RH (%)")
    plt.gca().invert_yaxis()
    plt.xlabel("Relative humidity (%)")
    plt.ylabel("Pressure (hPa)")
    plt.title("Equilibrium relative humidity profile")
    plt.grid(True, alpha=0.3)

    fig6 = plt.figure()
    plt.plot(cloud_w_gkg, p_hpa, label="Cloud water (g/kg)")
    plt.gca().invert_yaxis()
    plt.xlabel("Cloud condensed water (g/kg)")
    plt.ylabel("Pressure (hPa)")
    plt.title("Equilibrium cloud water profile")
    plt.grid(True, alpha=0.3)

    fig7 = plt.figure()
    plt.plot(cloud_frac, p_hpa, label="Cloud fraction")
    plt.gca().invert_yaxis()
    plt.xlabel("Cloud fraction")
    plt.ylabel("Pressure (hPa)")
    plt.title("Equilibrium cloud fraction profile")
    plt.grid(True, alpha=0.3)

    if save_dir is not None:
        save_dir.mkdir(parents=True, exist_ok=True)
        figs = {
            "time_precip_evap.png": fig1,
            "time_sst_ta.png": fig2,
            "time_toa_fluxes.png": fig3,
            "profile_temperature.png": fig4,
            "profile_rh.png": fig5,
            "profile_cloud_water.png": fig6,
            "profile_cloud_fraction.png": fig7,
        }
        for name, fig in figs.items():
            fig.savefig(save_dir / name, dpi=200, bbox_inches="tight")
        print(f"Saved {len(figs)} figures to: {save_dir.resolve()}")
    else:
        plt.show()


if __name__ == "__main__":
    main()

