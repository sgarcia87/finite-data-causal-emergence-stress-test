#!/usr/bin/env python3
"""Regenerate the triadic resolved-mode figure from committed CSV results."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "triadic" / "triadic_resolved_modes_summary.csv"
OUTPUT = ROOT / "assets" / "fig4_triadic_resolved_modes.png"


def main() -> None:
    data = pd.read_csv(RESULTS)
    figure, axis = plt.subplots(figsize=(7.2, 4.5))
    for tau, group in data.groupby("tau"):
        group = group.sort_values("samples_per_state")
        axis.plot(
            group["samples_per_state"],
            group["median_resolved"],
            marker="o",
            linewidth=2,
            label=rf"$\tau={tau:g}$",
        )
        axis.fill_between(
            group["samples_per_state"],
            group["q05_resolved"],
            group["q95_resolved"],
            alpha=0.15,
        )
    axis.axhline(8, color="#555555", linestyle="--", linewidth=1, label="exact 8-state quotient")
    axis.axhline(2, color="#999999", linestyle=":", linewidth=1.3, label="exact 2-state quotient")
    axis.set_xscale("log")
    axis.set_xlabel("transition samples per microstate")
    axis.set_ylabel("resolved spectral modes")
    axis.set_title("Resolved modes vary while the exact 48→8→2 hierarchy is fixed")
    axis.set_ylim(0, 50)
    axis.grid(alpha=0.2)
    axis.legend(frameon=False, fontsize=9)
    figure.tight_layout()
    figure.savefig(OUTPUT, dpi=180)
    plt.close(figure)


if __name__ == "__main__":
    main()
