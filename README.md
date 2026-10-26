# Fused-BiPLS: Adaptive Interval Selection for Near-Infrared Spectra

This repository contains the code for the following article:

> T. Bai, X. Peng, B. Xu, and P. He, "Adaptive Interval Selection for Near-Infrared Spectra via Fused Lasso and Backward Interval PLS," *Journal of Chemometrics* 40, no. 7 (2026): e70157. https://doi.org/10.1002/cem.70157

Fused-BiPLS is a two-stage wavelength-interval selection method for multivariate calibration of near-infrared (NIR) spectra. The first stage determines contiguous spectral segments; the second stage selects a subset of these segments for a partial least squares (PLS) model.

## Method

1. **Segmentation by the fused lasso.** A fused lasso penalty on the first-order differences of adjacent regression coefficients is solved along its exact solution path (`genlasso`). The regularization parameter λ is the value on the path that minimizes Mallows' Cp, with the noise variance estimated from the cross-validated error of a full-spectrum PLS model. The resulting piecewise-constant coefficient profile defines contiguous segments; segments narrower than three channels are merged with a neighbouring segment.
2. **Backward interval PLS (BiPLS).** Starting from all segments, the segment whose removal gives the largest decrease in 10-fold RMSECV is removed at each step. The procedure stops when no removal reduces RMSECV. The number of latent variables is chosen by the one-sigma rule, with the minimum-RMSECV value used as a fallback.

Spectra are corrected by multiplicative scatter correction (MSC) and mean-centred. The MSC reference spectrum and the centring vectors are estimated on the training set only and then applied to the test set.

## Repository structure

```text
├── code/
│   ├── benchmark_models.R          # Table 3: Fused-BiPLS and six baselines on six datasets
│   ├── effect_size.R               # Table 4: Cohen's d of Fused-BiPLS against each baseline
│   ├── sensitivity_analysis.R      # Data for Figure 4: sensitivity to the fusion parameter λ
│   ├── simulation_figure2.py       # Figure 2: coefficient estimates on simulated collinear spectra
│   ├── plot_figure3_intervals.py   # Figure 3: variables selected by each method at k = 0.7
│   └── plot_figure4_sensitivity.py # Figure 4: number of intervals and RMSECV against log(λ)
├── data/
│   └── README.md                   # Sources and expected file names of the datasets
├── requirements.txt                # Python packages for the figure scripts
└── LICENSE
```

## Requirements

- R (version 4.x) with the packages `genlasso`, `prospectr`, `pls`, `glmnet`, `R.matlab`, `readxl`, `ggplot2`, `scales` and `pbapply`:

  ```r
  install.packages(c("genlasso", "prospectr", "pls", "glmnet", "R.matlab",
                     "readxl", "ggplot2", "scales", "pbapply"))
  ```

  On Windows, `extrafont` is used, if installed, to register the Times New Roman font for the plots of `sensitivity_analysis.R`.

- Python 3 with the packages in `requirements.txt` (only for the figure scripts):

  ```bash
  pip install -r requirements.txt
  ```

## Data

The six datasets (Corn, Diesel, Meat, Milk, Soil and Tablet) are publicly available and are not included in this repository. `data/README.md` lists their sources, the file names expected by the scripts and the variables that are read. Place the files in `data/`.

## Usage

Run all scripts from the repository root. Outputs are written to `results/`.

```bash
# Table 3: 30 random 80/20 partitions x k in {0.5, 0.7, 1.0} for each dataset
Rscript code/benchmark_models.R

# Table 4: Cohen's d effect sizes
Rscript code/effect_size.R

# Figure 4: sensitivity analysis (k = 0.7, seed = 42), then the combined figure
Rscript code/sensitivity_analysis.R
python code/plot_figure4_sensitivity.py

# Figure 3: requires the output of benchmark_models.R
python code/plot_figure3_intervals.py

# Figure 2: simulated data, no input files required
python code/simulation_figure2.py
```

Notes:

- `benchmark_models.R` evaluates Fused-BiPLS, PLS, PCR, LASSO, MWPLS, MC-UVE and SCARS on identical training and test sets. Within each repetition, 80% of the samples form the calibration pool and 20% form the test set; the training subset, of size ⌊k · n<sub>pool</sub>⌋, is drawn from the pool with the Kennard–Stone algorithm. The hyperparameter ranges are given in Table 2 of the article. The runs are distributed over the available CPU cores. The exact fused lasso path is the most time-consuming step; in the article, the calibration time ranged from a few seconds for the smallest spectral matrix to about ninety minutes for the dataset with the largest number of channels.
- `effect_size.R` uses, by default, the RMSEP means and standard deviations reported in Table 3, which reproduces Table 4. Set `USE_BENCHMARK_OUTPUT <- TRUE` to compute the effect sizes from the summary files written by `benchmark_models.R`.
- `plot_figure3_intervals.py` plots, for each method, the repetition with the lowest test-set RMSEP at k = 0.7.
- `sensitivity_analysis.R` also writes a single-panel PDF for each dataset; `plot_figure4_sensitivity.py` produces the six-panel figure.

## Citation

```bibtex
@article{Bai2026FusedBiPLS,
  author  = {Bai, Tongyuan and Peng, Xiaoling and Xu, Baojun and He, Ping},
  title   = {Adaptive Interval Selection for Near-Infrared Spectra via Fused Lasso and Backward Interval {PLS}},
  journal = {Journal of Chemometrics},
  year    = {2026},
  volume  = {40},
  number  = {7},
  pages   = {e70157},
  doi     = {10.1002/cem.70157}
}
```

## License

This project is released under the MIT License (see `LICENSE`).
