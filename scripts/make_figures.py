#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
RESULTS=ROOT/"results"
ASSETS=ROOT/"assets"
ASSETS.mkdir(exist_ok=True)

# Figure 1
rep=pd.read_csv(RESULTS/"final_benchmark_representation_invariance.csv")
vals=[rep.base_clear_CE.mean(),rep.clone_clear_CE.mean()]
plt.figure(figsize=(6,4))
plt.bar(["Original 4-state","Equal 4x clone\n(16 observed states)"],vals)
plt.ylabel("Fraction classified as clear CE")
plt.ylim(0,1.05)
plt.title("Clear-CE status under exact behavior-preserving refinement")
plt.tight_layout()
plt.savefig(ASSETS/"fig1_representation_refinement.png",dpi=180)
plt.close()

# Figure 2
low=pd.read_csv(RESULTS/"final_benchmark_lowrank_summary.csv")
plt.figure(figsize=(7,4.5))
for fam,q in low.groupby("family"):
    rank=fam.split("_")[-1]
    plt.plot(q.samples_per_row,q.median_resolved_fraction,marker="o",label=f"true rank {rank}")
plt.xscale("log")
plt.ylim(0,1.05)
plt.xlabel("Transition samples per TPM row")
plt.ylabel("Median fraction of true modes resolved")
plt.title("Positive evidence for resolved modes grows with data")
plt.legend()
plt.tight_layout()
plt.savefig(ASSETS/"fig2_lowrank_resolved_fraction.png",dpi=180)
plt.close()

# Figure 3
full=pd.read_csv(RESULTS/"final_benchmark_fullrank_summary.csv")
plt.figure(figsize=(8,4.8))
for fam,q in full.groupby("family"):
    plt.plot(q.samples_per_row,q.median_resolved_lower_bound,marker="o",label=fam)
plt.xscale("log")
plt.xlabel("Transition samples per TPM row")
plt.ylabel("Median resolved-mode count")
plt.title("Full rank can coexist with only a few currently resolved modes")
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig(ASSETS/"fig3_fullrank_resolved_modes.png",dpi=180)
plt.close()
