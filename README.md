# Adaptive Interval Selection for Near-Infrared Spectra via Fused Lasso and Backward Interval PLS

This repository contains the R implementation of the Fused-BiPLS framework, a two-stage variable selection method designed for Near-Infrared (NIR) spectral analysis.

## Description
The Fused-BiPLS framework addresses the contiguous band structure and high collinearity inherent in NIR spectra. The method partitions the spectrum into intervals and selects predictive features through the following process:
1.  **Stage 1: Structural Segmentation**: Fused Lasso partitions the full spectrum into contiguous segments. The regularization parameter $\lambda$ is determined by Mallows' $C_p$ criterion.
2.  **Stage 2: Backward Elimination**: A backward elimination strategy using Partial Least Squares (BiPLS) removes intervals that do not reduce the cross-validation error.

## Methodology and Implementation
The implementation includes data preprocessing, model training, and comparative evaluation against established chemometric baselines.

### Preprocessing
Data undergoes Multiplicative Scatter Correction (MSC) and mean centering. Preprocessing parameters are derived exclusively from the calibration set to prevent data leakage.

### Comparative Baselines
The framework is evaluated against six models:
* Partial Least Squares (PLS)
* Principal Component Regression (PCR)
* Least Absolute Shrinkage and Selection Operator (LASSO)
* Moving-Window PLS (MWPLS)
* Monte-Carlo Uninformative Variable Elimination (MC-UVE)
* Stability Competitive Adaptive Reweighted Sampling (SCARS)

### Performance Metrics
Model performance is quantified using:
* Root Mean Square Error of Prediction (RMSEP)
* Coefficient of Determination ($R^2$)
* Cohen’s $d$ effect size for statistical comparison across 30 independent repetitions

## Repository Structure
* `benchmark_models.R`: Executes 30 independent random splits (80% calibration, 20% test) for all 7 methods across 6 benchmark datasets.
* `sensitivity_analysis.R`: Conducts a parallelized grid search on $\lambda$ to evaluate the trade-off between interval count and RMSECV.
* `/data`: Contains spectral datasets in `.mat`, `.csv`, and `.xlsx` formats.
* `/results`: Stores detailed CSV logs and visualization plots.

## Usage
1.  **Environment**: Requires R 4.x with packages `genlasso`, `prospectr`, `pls`, `glmnet`, `R.matlab`, and `ggplot2`.
2.  **Execution**: 
    * Run `benchmark_models.R` to reproduce comparative results for Corn, Diesel, Meat, Milk, Soil, and Tablet datasets.
    * Run `sensitivity_analysis.R` to generate $\lambda$ sensitivity curves for specific datasets.

## Affiliation
* **Author**: Tongyuan Bai
* **Institution**: Department of Statistics and Data Science, Beijing Normal University-Hong Kong Baptist University United International College / Hong Kong Baptist University
* **Supervisor**: Prof. Ping He
