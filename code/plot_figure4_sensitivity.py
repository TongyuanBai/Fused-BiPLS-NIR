# -*- coding: utf-8 -*-
"""
Figure 4: sensitivity of Fused-BiPLS to the fusion parameter lambda.

Reads ./results/lambda_rmsecv_intervals_<dataset>.csv written by
code/sensitivity_analysis.R and plots, against log(lambda), the number of fused
intervals (left axis, dashed) and the RMSECV after BiPLS (right axis, solid).
The star marks the lambda selected by Mallows' Cp.

Usage (from the repository root):
    python code/plot_figure4_sensitivity.py

Output: ./results/lambda_sensitivity_<dataset>_panel.pdf/.png (single panels)
        ./results/Figure_4.pdf and ./results/Figure_4.png      (3 x 2 grid)
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

# ---------- Style ----------
FONT_FAMILY = "DejaVu Sans"
plt.rcParams.update({
    "font.family":      FONT_FAMILY,
    "font.size":        11,
    "axes.linewidth":   1.0,
    "axes.edgecolor":   "#333333",
    "xtick.color":      "#333333",
    "ytick.color":      "#333333",
    "xtick.major.width": 0.9,
    "ytick.major.width": 0.9,
    "axes.labelcolor":  "#1A1A1A",
    "savefig.dpi":      2000,
})

# Colours
COLOR_INTV   = "#27517F"
COLOR_RMSECV = "#BD4A30"
COLOR_OPT    = "#1A1A1A"   # near-black

DATASETS = ["corn", "diesel", "meat", "milk", "soil", "tablet"]
REPO_ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = REPO_ROOT / "results"
OUTPUT_DIR = REPO_ROOT / "results"
os.makedirs(OUTPUT_DIR, exist_ok=True)


def load(ds):
    df = pd.read_csv(os.path.join(INPUT_DIR, f"lambda_rmsecv_intervals_{ds}.csv"))
    df = df.sort_values("Log_Lambda").reset_index(drop=True)
    log_lam = df["Log_Lambda"].values
    rmsecv  = df["RMSECV"].values
    intvs   = df["Segments_Fused"].values
    is_opt  = df["Is_Optimal_Cp"].values == "YES"
    opt_idx = int(np.where(is_opt)[0][0])
    return log_lam, rmsecv, intvs, opt_idx


def draw_panel(ax_left, log_lam, rmsecv, intvs, opt_idx,
               show_legend=True, label_fontsize=12, tick_fontsize=10.5):
    log_opt    = log_lam[opt_idx]
    rmsecv_opt = rmsecv[opt_idx]

    # Number of intervals (left)
    ax_left.plot(log_lam, intvs,
                 color=COLOR_INTV, linestyle="--", linewidth=1.5,
                 marker="^", markersize=4.0, markeredgecolor=COLOR_INTV,
                 markerfacecolor="white", markeredgewidth=1.0,
                 label="Number of intervals", zorder=2)
    ax_left.set_xlabel(r"$\log(\lambda)$", fontsize=label_fontsize)
    ax_left.set_ylabel("Number of intervals", color=COLOR_INTV, fontsize=label_fontsize)
    ax_left.tick_params(axis="y", colors=COLOR_INTV, labelsize=tick_fontsize)
    ax_left.tick_params(axis="x", labelsize=tick_fontsize)
    ax_left.spines["left"].set_color(COLOR_INTV)
    ax_left.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=6))

    # RMSECV (right)
    ax_right = ax_left.twinx()
    ax_right.plot(log_lam, rmsecv,
                  color=COLOR_RMSECV, linestyle="-", linewidth=1.6,
                  marker="o", markersize=4.0, markeredgecolor=COLOR_RMSECV,
                  markerfacecolor="white", markeredgewidth=1.0,
                  label="RMSECV", zorder=3)
    ax_right.set_ylabel("RMSECV", color=COLOR_RMSECV, fontsize=label_fontsize)
    ax_right.tick_params(axis="y", colors=COLOR_RMSECV, labelsize=tick_fontsize)
    ax_right.spines["right"].set_color(COLOR_RMSECV)

    # Optimal lambda: vertical line + star marker
    ax_left.axvline(log_opt, color=COLOR_OPT, linestyle=":", linewidth=0.9, zorder=1)
    ax_right.scatter([log_opt], [rmsecv_opt],
                     marker="*", s=160, color=COLOR_OPT,
                     edgecolor="white", linewidth=0.7, zorder=4,
                     label=r"Optimal $\lambda^\star$")

    ax_left.grid(axis="y", color="#EAEAEA", linewidth=0.6, zorder=0)
    ax_left.set_axisbelow(True)
    ax_left.set_xlim(log_lam.min() - 0.5, log_lam.max() + 0.5)

    if show_legend:
        h1, l1 = ax_left.get_legend_handles_labels()
        h2, l2 = ax_right.get_legend_handles_labels()
        legend = ax_left.legend(h1 + h2, l1 + l2,
                                loc="upper right", frameon=True,
                                framealpha=0.92, edgecolor="#999999",
                                fontsize=9.5, handlelength=2.0,
                                borderpad=0.4, labelspacing=0.4)
        legend.get_frame().set_linewidth(0.5)
    return ax_right


# =====================================================
# 1) Six individual figures
# =====================================================
for ds in DATASETS:
    log_lam, rmsecv, intvs, opt_idx = load(ds)
    fig, ax_left = plt.subplots(figsize=(6.4, 4.4))
    draw_panel(ax_left, log_lam, rmsecv, intvs, opt_idx,
               show_legend=True, label_fontsize=12.5, tick_fontsize=10.5)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f"lambda_sensitivity_{ds}_panel.pdf"),
                dpi=600, bbox_inches="tight", facecolor="white")
    plt.savefig(os.path.join(OUTPUT_DIR, f"lambda_sensitivity_{ds}_panel.png"),
                dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  individual panel: {ds}")


# =====================================================
# 2) Combined 3 x 2 grid (all six datasets, one figure)
# =====================================================
PANEL_LABELS = ["(a) Corn", "(b) Diesel", "(c) Meat",
                "(d) Milk", "(e) Soil",   "(f) Tablet"]

fig, axes = plt.subplots(3, 2, figsize=(13.5, 13.0))
axes = axes.flatten()

for i, ds in enumerate(DATASETS):
    log_lam, rmsecv, intvs, opt_idx = load(ds)
    draw_panel(axes[i], log_lam, rmsecv, intvs, opt_idx,
               show_legend=False, label_fontsize=12.5, tick_fontsize=11)
    axes[i].text(0.015, 1.02, PANEL_LABELS[i], transform=axes[i].transAxes,
                 fontsize=13, weight="bold", color="#1A1A1A",
                 ha="left", va="bottom")

# Single shared legend at the top of the whole figure
legend_handles = [
    plt.Line2D([0], [0], color=COLOR_INTV, linestyle="--", linewidth=1.5,
               marker="^", markersize=6, markerfacecolor="white",
               markeredgecolor=COLOR_INTV, markeredgewidth=1.0,
               label="Number of intervals"),
    plt.Line2D([0], [0], color=COLOR_RMSECV, linestyle="-", linewidth=1.6,
               marker="o", markersize=6, markerfacecolor="white",
               markeredgecolor=COLOR_RMSECV, markeredgewidth=1.0,
               label="RMSECV"),
    plt.Line2D([0], [0], marker="*", color="none", markerfacecolor=COLOR_OPT,
               markeredgecolor="white", markersize=14, markeredgewidth=0.7,
               label=r"Optimal $\lambda^\star$"),
]
fig.legend(handles=legend_handles, loc="upper center",
           bbox_to_anchor=(0.5, 0.995), ncol=3, frameon=False,
           fontsize=12.5, handlelength=2.5, columnspacing=2.5)

plt.subplots_adjust(left=0.06, right=0.94, top=0.95, bottom=0.05,
                    hspace=0.35, wspace=0.30)

plt.savefig(os.path.join(OUTPUT_DIR, "Figure_4.pdf"),
            dpi=600, bbox_inches="tight", facecolor="white")
plt.savefig(os.path.join(OUTPUT_DIR, "Figure_4.png"),
            dpi=300, bbox_inches="tight", facecolor="white")
plt.close(fig)
print("  combined 3 x 2 grid saved.")

print("\nAll figures saved to", OUTPUT_DIR)