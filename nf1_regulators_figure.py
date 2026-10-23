"""Two-panel figure of candidate NF1 exon 23a splicing regulators.

(A) Heatmap of log2(TPM + 1): GTEx EBV-transformed lymphocytes and whole blood
    (median TPM) beside NGN2 neurons (mean TPM across lines, per day).
(B) Time course per gene across NGN2 differentiation (mean +/- SD per line).

Inputs (written by the R scripts in this repository):
    nf1_regulators_gtex_summary.csv      from nf1_regulators_gtex_recount3.R
    nf1_regulators_neurons_summary.csv   from nf1_regulators_neurons_nfcore.R

Usage:
    pip install matplotlib numpy
    python nf1_regulators_figure.py [input_folder]

Writes nf1_regulators_figure.png (600 dpi) to the input folder.
Uses the Lato font if installed, otherwise the matplotlib default.
"""
import csv
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

GENES = ["NF1", "REST", "CELF3", "CELF4", "CELF5", "CELF6", "ELAVL2", "ELAVL3", "ELAVL4"]
TISSUES = [("Cells - EBV-transformed lymphocytes", "EBV\nlymphocytes"),
           ("Whole Blood", "Whole\nblood")]

folder = sys.argv[1] if len(sys.argv) > 1 else "."
gtex_file = os.path.join(folder, "nf1_regulators_gtex_summary.csv")
neur_file = os.path.join(folder, "nf1_regulators_neurons_summary.csv")
out_file = os.path.join(folder, "nf1_regulators_figure.png")

if any(f.name == "Lato" for f in fm.fontManager.ttflist):
    plt.rcParams["font.family"] = "Lato"

# Read summaries -------------------------------------------------------------------------
with open(gtex_file, newline="") as f:
    gtex = {(r["tissue"], r["gene"]): float(r["median_tpm"]) for r in csv.DictReader(f)}
with open(neur_file, newline="") as f:
    neur = {(r["gene"], r["line"], int(float(r["day"]))):
            (float(r["mean_tpm"]), float(r["sd_tpm"]) if r["sd_tpm"] not in ("", "NA") else 0.0)
            for r in csv.DictReader(f)}

genes = [g for g in GENES if (TISSUES[0][0], g) in gtex and any(k[0] == g for k in neur)]
lines = sorted({k[1] for k in neur})
days = sorted({k[2] for k in neur})

def neuron_mean(g, d):
    v = [neur[(g, l, d)][0] for l in lines if (g, l, d) in neur]
    return np.mean(v) if v else np.nan

M = np.array([[gtex[(t, g)] for t, _ in TISSUES] + [neuron_mean(g, d) for d in days]
              for g in genes])
cols = [lab for _, lab in TISSUES] + [f"Day {d}" for d in days]

fig = plt.figure(figsize=(13, 6.2))

# Panel A: heatmap -----------------------------------------------------------------------
ax = fig.add_axes([0.07, 0.13, 0.4, 0.74])
Z = np.log2(M + 1)
im = ax.imshow(Z, cmap="Blues", aspect="auto", vmin=0, vmax=max(9.5, np.nanmax(Z)))
ax.set_xticks(range(len(cols)))
ax.set_xticklabels(cols, fontsize=10)
ax.set_yticks(range(len(genes)))
ax.set_yticklabels(genes, style="italic", fontsize=11)
for i in range(len(genes)):
    for j in range(len(cols)):
        v = M[i, j]
        s = f"{v:.0f}" if v >= 10 else (f"{v:.1f}" if v >= 0.1 else "<0.1")
        ax.text(j, i, s, ha="center", va="center", fontsize=9,
                color="white" if Z[i, j] > 5.5 else "black")
nt = len(TISSUES)
ax.axvline(nt - 0.5, color="black", lw=1.5)
ax.text((nt - 1) / 2, -0.9, "GTEx", ha="center", fontsize=11, fontweight="bold")
ax.text(nt + (len(days) - 1) / 2, -0.9, "NGN2 neurons", ha="center", fontsize=11, fontweight="bold")
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
cb.set_label("log$_2$(TPM + 1)", fontsize=10)
cb.outline.set_visible(False)
fig.text(0.01, 0.95, "A", fontsize=16, fontweight="bold")

# Panel B: time course -------------------------------------------------------------------
fig.text(0.53, 0.95, "B", fontsize=16, fontweight="bold")
ncol = 3
nrow = int(np.ceil(len(genes) / ncol))
gs = fig.add_gridspec(nrow, ncol, left=0.56, right=0.99, bottom=0.1, top=0.92,
                      hspace=0.55, wspace=0.35)
styles = [dict(color="black", marker="o", mfc="black"),
          dict(color="grey", marker="s", mfc="white")]
for k, g in enumerate(genes):
    a = fig.add_subplot(gs[k // ncol, k % ncol])
    for li, l in enumerate(lines):
        dd = [d for d in days if (g, l, d) in neur]
        a.errorbar(dd, [neur[(g, l, d)][0] for d in dd], yerr=[neur[(g, l, d)][1] for d in dd],
                   capsize=2, lw=1.2, ms=4, label=l.replace("_NGN2", ""),
                   **styles[li % len(styles)])
    a.set_title(g, style="italic", fontsize=11)
    a.set_ylim(bottom=0)
    a.set_xticks(days)
    a.tick_params(labelsize=8)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
    if k % ncol == 0:
        a.set_ylabel("TPM", fontsize=9)
    if k // ncol == nrow - 1:
        a.set_xlabel("Day", fontsize=9)
    if k == ncol - 1:
        a.legend(frameon=False, fontsize=8, loc="lower right")

plt.savefig(out_file, dpi=600)
print("Wrote", out_file)
