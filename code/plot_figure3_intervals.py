# -*- coding: utf-8 -*-
"""
Figure 3: spectra and selected variables of five variable-selection methods
(LASSO, MWPLS, MC-UVE, SCARS and Fused-BiPLS) on the six datasets at k = 0.7.

For each dataset and method, the repetition with the lowest test-set RMSEP at
k = 0.7 is taken from ./results/<dataset>_All_Methods_Detailed.csv, the file
written by code/benchmark_models.R. The selected runs are also saved to
./results/best_performance_k0.7_selected_methods.csv.

Usage (from the repository root):
    python code/plot_figure3_intervals.py

Output: ./results/Figure_3.pdf and ./results/Figure_3.png
"""

import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io as sio
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.lines import Line2D

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = Path(os.environ.get("FUSED_BIPLS_DATA_DIR", REPO_ROOT / "data"))
RESULTS_DIR = REPO_ROOT / "results"
K_PROP = 0.7

FUSED_METHOD = "Ours (Fused-biPLS)"          # method label used in benchmark_models.R
POINT_METHODS = ["LASSO", "MWPLS", "MC-UVE", "SCARS"]


# ==========================================
# 1. Data loader (same files and columns as the R scripts)
# ==========================================
def data_file(name):
    """Return the path of a data file; the match is case-insensitive."""
    for f in DATA_DIR.iterdir():
        if f.name.lower() == name.lower():
            return f
    raise FileNotFoundError(f"'{name}' not found in {DATA_DIR}. See data/README.md.")


def load_dataset(ds_name):
    if ds_name == "corn":
        mat = sio.loadmat(data_file("corn.mat"))
        y = mat["propvals"]["data"][0][0][:, 0]
        X = mat["m5spec"]["data"][0][0]
    elif ds_name == "diesel":
        df = pd.read_csv(data_file("Filtered_Diesel_Spec_Data.csv"))
        y = df.iloc[:, 0].values
        X = df.iloc[:, 1:].values
    elif ds_name == "meat":
        df = pd.read_excel(data_file("Tecator.xlsx"))
        y = df.iloc[:, 0].values
        X = df.iloc[:, 1:].values
    elif ds_name == "milk":
        mat = sio.loadmat(data_file("milk.mat"))
        y = mat["protein"].flatten()
        X = mat["X"]
    elif ds_name == "soil":
        mat = sio.loadmat(data_file("NIRsoil.mat"))
        y = mat["soilref"]["data"][0][0][:, 0]
        X = mat["soil"]["data"][0][0]
    elif ds_name == "tablet":
        mat = sio.loadmat(data_file("NIRdata_tablets.mat"))
        y = mat["Matrix"][:, 0]
        X = mat["Matrix"][:, 3:407]          # columns 4:407 in R (404 channels)
    else:
        raise ValueError(ds_name)
    return X, y


# ==========================================
# 2. Best run per method at k = 0.7
# ==========================================
def best_runs(ds_name):
    detailed = pd.read_csv(RESULTS_DIR / f"{ds_name}_All_Methods_Detailed.csv")
    sub = detailed[np.isclose(detailed["k_prop"], K_PROP)]
    rows = []
    for method in POINT_METHODS + [FUSED_METHOD]:
        m = sub[sub["Method"] == method]
        if m.empty:
            continue
        rows.append(m.loc[m["RMSE_Test"].idxmin()])
    out = pd.DataFrame(rows)
    out.insert(0, "Dataset", ds_name)
    return out


def parse_selection(best):
    """1-based indices in the CSV are converted to 0-based indices."""
    selected_wavelengths, fused_intervals = {}, []
    for _, row in best.iterrows():
        method = row["Method"]
        if method == FUSED_METHOD:
            intervals = str(row.get("Selected_Intervals", ""))
            for interval in intervals.split(";"):
                if "-" in interval:
                    s, e = interval.split("-")
                    fused_intervals.append((int(s.strip()) - 1, int(e.strip())))
        else:
            vars_str = str(row["Selected_Vars"])
            selected_wavelengths[method] = [int(v) - 1 for v in vars_str.split(";") if v.strip().isdigit()]
    return selected_wavelengths, fused_intervals


# ==========================================
# 3. Global style
# ==========================================
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "axes.linewidth": 1.1,
    "axes.labelsize": 11,
    "xtick.labelsize": 9.5,
    "ytick.labelsize": 9.5,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.size": 4,
    "ytick.major.size": 4,
    "legend.fontsize": 11,
})

dataset_names = ["corn", "diesel", "meat", "milk", "soil", "tablet"]
dataset_display = {"corn": "Corn", "diesel": "Diesel", "meat": "Meat",
                   "milk": "Milk", "soil": "Soil", "tablet": "Tablet"}
dataset_info = {
    "corn":   {"start": 1100, "end": 2498,  "unit": "Wavelength (nm)"},
    "diesel": {"start": 750,  "end": 1550,  "unit": "Wavelength (nm)"},
    "meat":   {"start": 850,  "end": 1050,  "unit": "Wavelength (nm)"},
    "milk":   {"start": 1000, "end": 7224,  "unit": r"Wavenumber (cm$^{-1}$)"},
    "soil":   {"start": 400,  "end": 2500,  "unit": "Wavelength (nm)"},
    "tablet": {"start": 7400, "end": 10500, "unit": r"Wavenumber (cm$^{-1}$)"},
}
panel_labels = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]
markers = ["o", "s", "^", "D", "v", "p", "X", "*"]
method_color_cycle = plt.cm.tab10.colors
FUSED_COLOR = "#D62728"


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # ==========================================
    # 4. 3 x 2 grid: spectra (top) and selected variables (bottom) per dataset
    # ==========================================
    fig = plt.figure(figsize=(15, 16))
    outer = GridSpec(3, 2, figure=fig, left=0.06, right=0.985, top=0.945, bottom=0.045,
                     hspace=0.38, wspace=0.20)
    global_method_color = {}
    all_best = []

    for idx, ds in enumerate(dataset_names):
        r, c = divmod(idx, 2)
        inner = GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[r, c],
                                        height_ratios=[2.5, 1.2], hspace=0.08)
        ax_spec = fig.add_subplot(inner[0])
        ax_sel = fig.add_subplot(inner[1], sharex=ax_spec)

        X, _ = load_dataset(ds)
        info = dataset_info[ds]
        wavelengths = np.linspace(info["start"], info["end"], X.shape[1])

        best = best_runs(ds)
        all_best.append(best)
        sel_idx, fused_intervals = parse_selection(best)
        sel_dict = {m: wavelengths[i] for m, i in sel_idx.items()}

        # Top panel: spectra
        cmap = LinearSegmentedColormap.from_list("blue_orange", ["#4C72B0", "#DD8452"])
        colors = cmap(np.linspace(0, 1, X.shape[0]))
        for i in range(X.shape[0]):
            ax_spec.plot(wavelengths, X[i, :], color=colors[i], linewidth=0.6, alpha=0.75)
        ax_spec.set_ylabel("Absorbance", fontweight="bold")
        ax_spec.set_xlim(wavelengths[0], wavelengths[-1])
        ax_spec.tick_params(labelbottom=False)
        ax_spec.text(0.0, 1.02, f"{panel_labels[idx]} {dataset_display[ds]}",
                     transform=ax_spec.transAxes, fontsize=13, fontweight="bold",
                     ha="left", va="bottom")

        # Bottom panel: one lane per method; Fused-BiPLS on the lowest lane
        method_names = list(sel_dict.keys())
        n_methods = len(method_names)
        for row_idx, method in enumerate(method_names):
            y_lane = n_methods - row_idx + 1
            wl = sel_dict[method]
            if method not in global_method_color:
                global_method_color[method] = method_color_cycle[len(global_method_color) % 10]
            ax_sel.scatter(wl, np.full_like(wl, y_lane, dtype=float),
                           marker=markers[row_idx % len(markers)], s=22,
                           color=global_method_color[method], edgecolor="black",
                           linewidth=0.3, zorder=3, alpha=0.9)

        for (s, e) in fused_intervals:
            s_safe = min(s, len(wavelengths) - 1)
            e_safe = min(e, len(wavelengths) - 1)
            ax_sel.hlines(y=1, xmin=wavelengths[s_safe], xmax=wavelengths[e_safe],
                          colors=FUSED_COLOR, linewidth=5.0, zorder=4, capstyle="butt")

        ax_sel.set_yticks(list(range(1, n_methods + 2)))
        ax_sel.set_yticklabels(["Fused-BiPLS"] + list(reversed(method_names)),
                               fontsize=10.5, fontweight="bold")
        ax_sel.tick_params(axis="y", pad=3)
        ax_sel.set_ylim(0.4, n_methods + 1.6)
        ax_sel.set_xlabel(info["unit"], fontweight="bold")
        ax_sel.set_ylabel("Method", fontweight="bold")
        ax_sel.grid(axis="x", color="#EEEEEE", linewidth=0.6, zorder=0)
        ax_sel.set_axisbelow(True)

    # ==========================================
    # 5. Shared legend
    # ==========================================
    handles = [
        Line2D([0], [0], marker=markers[i % len(markers)], linestyle="None", markersize=8,
               markerfacecolor=col, markeredgecolor="black", markeredgewidth=0.3, label=method)
        for i, (method, col) in enumerate(global_method_color.items())
    ]
    handles.append(Line2D([0], [0], color=FUSED_COLOR, linewidth=5, label="Fused-BiPLS"))
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.995),
               ncol=min(len(handles), 7), frameon=False, fontsize=11.5,
               handletextpad=0.5, columnspacing=1.6)

    # ==========================================
    # 6. Save
    # ==========================================
    pd.concat(all_best, ignore_index=True).to_csv(
        RESULTS_DIR / "best_performance_k0.7_selected_methods.csv", index=False)
    fig.savefig(RESULTS_DIR / "Figure_3.pdf", dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(RESULTS_DIR / "Figure_3.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved {RESULTS_DIR / 'Figure_3.pdf'}")


if __name__ == "__main__":
    main()
