###############################################################################
# Fused-BiPLS: Cohen's d effect sizes against each baseline (Table 4)
#
#   d = (mean RMSEP of baseline - mean RMSEP of Fused-BiPLS) / pooled SD,
#   pooled SD = sqrt((sd_baseline^2 + sd_FusedBiPLS^2) / 2),
# computed from the mean and standard deviation of the test-set RMSEP over the
# 30 repetitions at each calibration proportion k.
#   d > 0 : Fused-BiPLS has the lower RMSEP
#   d < 0 : the baseline has the lower RMSEP
#
# Two input sources are supported:
#   USE_BENCHMARK_OUTPUT <- FALSE  (default) use the RMSEP mean and SD reported
#                                  in Table 3 of the paper, which reproduces
#                                  Table 4 exactly;
#   USE_BENCHMARK_OUTPUT <- TRUE   read ./results/<dataset>_All_Methods_Summary.csv
#                                  written by code/benchmark_models.R.
#
# Usage (from the repository root):
#   Rscript code/effect_size.R
#
# Output: ./results/effect_size_results.csv  (long format)
#         ./results/effect_size_table.tex    (LaTeX table with shading by magnitude)
#
# Shading in the LaTeX table:
#   d >= 0.8           deep blue    + two up-triangles   (large, Fused-BiPLS lower)
#   0.5 <= d < 0.8     medium blue  + one up-triangle    (medium)
#   0.2 <= d < 0.5     light blue                        (small)
#   |d| < 0.2          no shading                        (negligible)
#  -0.5 < d <= -0.2    light red                         (small, baseline lower)
#  -0.8 < d <= -0.5    medium red   + one down-triangle  (medium)
#   d <= -0.8          deep red     + two down-triangles (large)
###############################################################################

USE_BENCHMARK_OUTPUT <- FALSE
RESULTS_DIR <- "./results"   # location of the benchmark summary files
OUTPUT_DIR  <- "./results"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- RMSEP mean and SD reported in Table 3 of the paper ----------------------
# Format: dataset -> data.frame with columns Method, k, Mean, SD
manuscript_table <- list(
  Corn = data.frame(
    Method = rep(c("PLS","PCR","LASSO","MWPLS","MC-UVE","SCARS","Fused-BiPLS"), each = 3),
    k      = rep(c(0.5, 0.7, 1.0), times = 7),
    Mean   = c(0.1830, 0.1742, 0.1655,   # PLS
               0.1717, 0.1645, 0.1582,   # PCR
               0.1575, 0.1666, 0.1476,   # LASSO
               0.2584, 0.2573, 0.2428,   # MWPLS
               0.1806, 0.1648, 0.1545,   # MC-UVE
               0.1782, 0.1656, 0.1527,   # SCARS
               0.1932, 0.1574, 0.1466),  # Fused-BiPLS
    SD     = c(0.0338, 0.0276, 0.0318,
               0.0337, 0.0329, 0.0356,
               0.0319, 0.0339, 0.0249,
               0.0581, 0.0556, 0.0380,
               0.0386, 0.0319, 0.0335,
               0.0367, 0.0334, 0.0308,
               0.0432, 0.0320, 0.0283),
    stringsAsFactors = FALSE
  ),
  Diesel = data.frame(
    Method = rep(c("PLS","PCR","LASSO","MWPLS","MC-UVE","SCARS","Fused-BiPLS"), each = 3),
    k      = rep(c(0.5, 0.7, 1.0), times = 7),
    Mean   = c(0.0018, 0.0017, 0.0016,
               0.0030, 0.0029, 0.0028,
               0.0032, 0.0032, 0.0032,
               0.0033, 0.0032, 0.0032,
               0.0026, 0.0025, 0.0026,
               0.0026, 0.0025, 0.0026,
               0.0015, 0.0014, 0.0013),
    SD     = c(0.0003, 0.0002, 0.0002,
               0.0003, 0.0003, 0.0003,
               0.0003, 0.0003, 0.0003,
               0.0002, 0.0002, 0.0003,
               0.0005, 0.0004, 0.0005,
               0.0004, 0.0004, 0.0005,
               0.0001, 0.0001, 0.0001),
    stringsAsFactors = FALSE
  ),
  Meat = data.frame(
    Method = rep(c("PLS","PCR","LASSO","MWPLS","MC-UVE","SCARS","Fused-BiPLS"), each = 3),
    k      = rep(c(0.5, 0.7, 1.0), times = 7),
    Mean   = c(2.3644, 2.3308, 2.3194,
               2.6302, 2.4717, 2.4477,
               2.4153, 2.3998, 2.4346,
               2.5879, 2.3496, 2.2834,
               2.3271, 2.2516, 2.3091,
               2.3577, 2.2457, 2.2555,
               2.3045, 2.2021, 2.1946),
    SD     = c(0.2636, 0.3019, 0.3394,
               0.4169, 0.3686, 0.3685,
               0.3411, 0.3955, 0.3101,
               0.3039, 0.2739, 0.2819,
               0.2919, 0.2533, 0.2594,
               0.2334, 0.2679, 0.2586,
               0.2709, 0.2727, 0.2689),
    stringsAsFactors = FALSE
  ),
  Milk = data.frame(
    Method = rep(c("PLS","PCR","LASSO","MWPLS","MC-UVE","SCARS","Fused-BiPLS"), each = 3),
    k      = rep(c(0.5, 0.7, 1.0), times = 7),
    Mean   = c(0.3052, 0.2897, 0.2650,
               0.2626, 0.2363, 0.2043,
               0.3026, 0.2986, 0.2803,
               0.3316, 0.3030, 0.3067,
               0.2668, 0.2343, 0.1997,
               0.2654, 0.2367, 0.1986,
               0.2552, 0.2435, 0.2435),
    SD     = c(0.2061, 0.1951, 0.1912,
               0.1657, 0.1136, 0.0995,
               0.1512, 0.1562, 0.1666,
               0.1416, 0.1410, 0.1454,
               0.1171, 0.1178, 0.1053,
               0.1209, 0.1176, 0.1059,
               0.1301, 0.1748, 0.1388),
    stringsAsFactors = FALSE
  ),
  Soil = data.frame(
    Method = rep(c("PLS","PCR","LASSO","MWPLS","MC-UVE","SCARS","Fused-BiPLS"), each = 3),
    k      = rep(c(0.5, 0.7, 1.0), times = 7),
    Mean   = c(3.6692, 3.1961, 2.8409,
               3.2094, 3.0413, 2.8243,
               3.5799, 3.4823, 3.4327,
               3.0545, 2.8453, 2.2716,
               3.5962, 3.2661, 3.3900,
               3.5456, 3.2826, 3.4016,
               2.8913, 2.7463, 2.4609),
    SD     = c(0.7586, 0.7047, 0.6210,
               0.6324, 0.5984, 0.5979,
               0.5051, 0.5504, 0.5629,
               0.5991, 0.6044, 0.4309,
               0.9198, 0.9648, 0.8514,
               1.0411, 0.8831, 0.9019,
               0.8297, 0.9723, 0.6095),
    stringsAsFactors = FALSE
  ),
  Tablet = data.frame(
    Method = rep(c("PLS","PCR","LASSO","MWPLS","MC-UVE","SCARS","Fused-BiPLS"), each = 3),
    k      = rep(c(0.5, 0.7, 1.0), times = 7),
    Mean   = c(0.3686, 0.3568, 0.3580,
               0.3660, 0.3597, 0.3550,
               0.3485, 0.3446, 0.3396,
               0.3535, 0.3432, 0.3383,
               0.3575, 0.3431, 0.3340,
               0.3567, 0.3435, 0.3338,
               0.3480, 0.3427, 0.3404),
    SD     = c(0.0308, 0.0328, 0.0335,
               0.0290, 0.0322, 0.0276,
               0.0338, 0.0294, 0.0290,
               0.0310, 0.0301, 0.0343,
               0.0272, 0.0240, 0.0256,
               0.0284, 0.0235, 0.0244,
               0.0284, 0.0284, 0.0305),
    stringsAsFactors = FALSE
  )
)

# --- Optional: build the same table from the benchmark output ---------------
read_benchmark_summary <- function(ds_label) {
  f <- file.path(RESULTS_DIR, paste0(tolower(ds_label), "_All_Methods_Summary.csv"))
  if (!file.exists(f)) stop(sprintf("Missing %s; run code/benchmark_models.R first.", f))
  sm <- read.csv(f, stringsAsFactors = FALSE, check.names = TRUE)
  sm$Method[sm$Method == "Ours (Fused-biPLS)"] <- "Fused-BiPLS"
  data.frame(Method = sm$Method, k = sm$k_prop,
             Mean = sm$RMSE_Test.mean, SD = sm$RMSE_Test.sd,
             stringsAsFactors = FALSE)
}

if (USE_BENCHMARK_OUTPUT) {
  manuscript_table <- setNames(lapply(names(manuscript_table), read_benchmark_summary),
                               names(manuscript_table))
}

BASELINES <- c("PLS", "PCR", "LASSO", "MWPLS", "MC-UVE", "SCARS")
REFERENCE <- "Fused-BiPLS"
KS        <- c(0.5, 0.7, 1.0)

# --- Helpers ----------------------------------------------------------------
cohens_d <- function(mu_b, sd_b, mu_r, sd_r) {
  pooled_sd <- sqrt((sd_b^2 + sd_r^2) / 2)
  if (pooled_sd == 0) return(NA_real_)
  (mu_b - mu_r) / pooled_sd
}

# Map d to: (color command, symbol). Uses xcolor + cellcolor.
# We use 5 saturation levels per direction.
d_to_format <- function(d) {
  if (is.na(d)) return(list(cell = "white", text = "--", symbol = ""))
  ad <- abs(d)
  if (d >= 0.8) {
    list(cell = "FBdeep",   text = sprintf("\\textbf{%.2f}", d), symbol = "$\\blacktriangle\\blacktriangle$")
  } else if (d >= 0.5) {
    list(cell = "FBmedium", text = sprintf("\\textbf{%.2f}", d), symbol = "$\\blacktriangle$")
  } else if (d >= 0.2) {
    list(cell = "FBlight",  text = sprintf("%.2f", d), symbol = "")
  } else if (d >  -0.2) {
    list(cell = "white",    text = sprintf("%.2f", d), symbol = "")
  } else if (d >  -0.5) {
    list(cell = "BLlight",  text = sprintf("%.2f", d), symbol = "")
  } else if (d >  -0.8) {
    list(cell = "BLmedium", text = sprintf("\\textbf{%.2f}", d), symbol = "$\\blacktriangledown$")
  } else {
    list(cell = "BLdeep",   text = sprintf("\\textbf{%.2f}", d), symbol = "$\\blacktriangledown\\blacktriangledown$")
  }
}

# --- Compute all d values ----------------------------------------------------
results <- list()
for (ds in names(manuscript_table)) {
  df <- manuscript_table[[ds]]
  for (kv in KS) {
    ref <- df[df$Method == REFERENCE & df$k == kv, ]
    mu_r <- ref$Mean
    sd_r <- ref$SD
    for (m in BASELINES) {
      base <- df[df$Method == m & df$k == kv, ]
      mu_b <- base$Mean
      sd_b <- base$SD
      d <- cohens_d(mu_b, sd_b, mu_r, sd_r)
      results[[length(results) + 1]] <- data.frame(
        Dataset = ds, k = kv, Baseline = m,
        Mean_Base = mu_b, SD_Base = sd_b,
        Mean_FBPLS = mu_r, SD_FBPLS = sd_r,
        Cohens_d = round(d, 3),
        stringsAsFactors = FALSE
      )
    }
  }
}
result_df <- do.call(rbind, results)
result_df <- result_df[order(factor(result_df$Dataset,
                                    levels = c("Corn","Diesel","Meat","Milk","Soil","Tablet")),
                             result_df$k), ]

# Save raw CSV
write.csv(result_df, file.path(OUTPUT_DIR, "effect_size_results.csv"), row.names = FALSE)

# --- Build LaTeX table -------------------------------------------------------
# We organize the table as: rows = (Dataset, k); cols = each baseline.
# Each cell shows: shaded background + d value + optional symbol (above the number).

# Required LaTeX preamble snippet (we'll print it as a comment at the top of the .tex)
preamble_note <- c(
  "% ---------------------------------------------------------------------------",
  "% LaTeX preamble requirements (add to your main .tex if not already present):",
  "%   \\usepackage[table]{xcolor}",
  "%   \\usepackage{colortbl}",
  "%   \\usepackage{booktabs}",
  "%   \\usepackage{multirow}",
  "%   \\usepackage{makecell}     % for line breaks inside cells",
  "%",
  "% Color palette (also defined inside the table for self-containment):",
  "%   FBdeep / FBmedium / FBlight  : blue scale -- Fused-BiPLS advantage",
  "%   BLdeep / BLmedium / BLlight  : red scale  -- Baseline advantage",
  "% ---------------------------------------------------------------------------",
  ""
)

# Color definitions block (place once, can be moved to preamble of main doc)
color_defs <- c(
  "% --- Color palette (you may move these to the document preamble) ---",
  "\\definecolor{FBdeep}  {RGB}{47, 84,150}    % Fused-BiPLS large advantage  (deep blue)",
  "\\definecolor{FBmedium}{RGB}{157,195,230}   % medium advantage             (medium blue)",
  "\\definecolor{FBlight} {RGB}{222,235,247}   % small advantage              (pale blue)",
  "\\definecolor{BLlight} {RGB}{252,228,214}   % small baseline advantage     (pale red)",
  "\\definecolor{BLmedium}{RGB}{248,203,173}   % medium baseline advantage    (medium red)",
  "\\definecolor{BLdeep}  {RGB}{197, 90, 17}   % large baseline advantage     (deep red)",
  ""
)

# Body
tex_lines <- c(
  preamble_note,
  color_defs,
  "\\begin{table*}[htbp]",
  "\\centering",
  "\\caption{Cohen's $d$ effect sizes of Fused-BiPLS against each baseline on the held-out test set RMSEP. Positive $d$ values indicate that Fused-BiPLS achieves lower RMSEP. Background colour intensity scales with effect magnitude; symbols flag medium and large effects.}",
  "\\label{tab:effect_size}",
  "\\footnotesize",
  "\\setlength{\\tabcolsep}{4pt}",
  "\\renewcommand{\\arraystretch}{1.35}",
  "\\begin{tabular}{l c c c c c c c}",
  "\\toprule",
  "Dataset & $k$ & vs PLS & vs PCR & vs LASSO & vs MWPLS & vs MC-UVE & vs SCARS \\\\",
  "\\midrule"
)

# Build rows
for (ds in c("Corn","Diesel","Meat","Milk","Soil","Tablet")) {
  for (i in seq_along(KS)) {
    kv <- KS[i]
    row_d <- result_df[result_df$Dataset == ds & result_df$k == kv, ]
    # First column: dataset name only on middle row of each block
    ds_label <- if (i == 2) sprintf("\\textbf{%s}", ds) else ""
    cells <- c()
    for (m in BASELINES) {
      d_val <- row_d$Cohens_d[row_d$Baseline == m]
      f <- d_to_format(d_val)
      # cell content: small symbol on top, d value below (if symbol exists)
      if (nchar(f$symbol) > 0) {
        content <- sprintf("\\makecell{\\scriptsize %s\\\\[-1pt] %s}", f$symbol, f$text)
      } else {
        content <- f$text
      }
      cells <- c(cells, sprintf("\\cellcolor{%s}%s", f$cell, content))
    }
    tex_lines <- c(tex_lines,
                   sprintf("%s & %.1f & %s \\\\",
                           ds_label, kv, paste(cells, collapse = " & ")))
  }
  if (ds != "Tablet") {
    tex_lines <- c(tex_lines, "\\addlinespace[2pt]")
  }
}

tex_lines <- c(
  tex_lines,
  "\\bottomrule",
  "\\end{tabular}",
  "",
  "\\vspace{1ex}",
  "{\\raggedright \\footnotesize \\textbf{Note:} Effect sizes follow the conventional interpretation of Cohen \\cite{cohen1988statistical}. Cell colour intensity (blue $\\to$ Fused-BiPLS advantage; red $\\to$ baseline advantage) and symbols flag effect magnitude:",
  "$\\blacktriangle\\blacktriangle$ large advantage of Fused-BiPLS ($d \\geq 0.8$);",
  "$\\blacktriangle$ medium advantage ($0.5 \\leq d < 0.8$);",
  "no symbol with light-blue cell: small advantage ($0.2 \\leq d < 0.5$);",
  "white cell: negligible difference ($|d| < 0.2$);",
  "light-red cell: small baseline advantage ($-0.5 < d \\leq -0.2$);",
  "$\\blacktriangledown$ medium baseline advantage ($-0.8 < d \\leq -0.5$);",
  "$\\blacktriangledown\\blacktriangledown$ large baseline advantage ($d \\leq -0.8$).\\par}",
  "\\end{table*}"
)

writeLines(tex_lines, file.path(OUTPUT_DIR, "effect_size_table.tex"))

# --- Console summary ---------------------------------------------------------
cat("\n=== Cohen's d summary ===\n\n")
n_total      <- nrow(result_df)
n_large_pos  <- sum(result_df$Cohens_d >=  0.8, na.rm = TRUE)
n_med_pos    <- sum(result_df$Cohens_d >=  0.5 & result_df$Cohens_d < 0.8, na.rm = TRUE)
n_small_pos  <- sum(result_df$Cohens_d >=  0.2 & result_df$Cohens_d < 0.5, na.rm = TRUE)
n_neg_neg    <- sum(abs(result_df$Cohens_d) < 0.2, na.rm = TRUE)
n_small_neg  <- sum(result_df$Cohens_d <= -0.2 & result_df$Cohens_d > -0.5, na.rm = TRUE)
n_med_neg    <- sum(result_df$Cohens_d <= -0.5 & result_df$Cohens_d > -0.8, na.rm = TRUE)
n_large_neg  <- sum(result_df$Cohens_d <= -0.8, na.rm = TRUE)

cat(sprintf("Total cells:                              %d\n", n_total))
cat(sprintf("Large F-BiPLS advantage  (d >= 0.8):     %d\n", n_large_pos))
cat(sprintf("Medium F-BiPLS advantage (0.5 to 0.8):   %d\n", n_med_pos))
cat(sprintf("Small F-BiPLS advantage  (0.2 to 0.5):   %d\n", n_small_pos))
cat(sprintf("Negligible difference    (|d| < 0.2):    %d\n", n_neg_neg))
cat(sprintf("Small baseline advantage (-0.5 to -0.2): %d\n", n_small_neg))
cat(sprintf("Medium baseline advantage(-0.8 to -0.5): %d\n", n_med_neg))
cat(sprintf("Large baseline advantage (d <= -0.8):    %d\n", n_large_neg))

cat(sprintf("\nFiles saved:\n  %s\n  %s\n",
            file.path(OUTPUT_DIR, "effect_size_results.csv"),
            file.path(OUTPUT_DIR, "effect_size_table.tex")))