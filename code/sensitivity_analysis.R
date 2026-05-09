# ==============================================================================
# Fused-BiPLS Sensitivity Analysis - Parallelized Edition
# ==============================================================================

library(parallel)
library(genlasso)
library(prospectr)
library(pls)
library(R.matlab)
library(ggplot2)
library(scales)
library(readxl)
library(pbapply)   # For parallel progress bar
library(extrafont) # For publication-quality plot fonts

# ================= 1. Global Setup & Fonts =================
# Note: Run font_import(prompt = FALSE) once if 'Times New Roman' is not recognized.
loadfonts(device = "win", quiet = TRUE) 

MAX_LVS <- 20    # Maximum number of latent variables for BiPLS
out_dir <- "./results"
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

dataset_names <- c("corn", "diesel", "meat", "milk", "soil", "tablet")

# ================= 2. Data Loading Function =================
load_dataset <- function(ds_name) {
  data_dir <- "./data"
  
  if (ds_name == "corn") {
    mat_data <- readMat(file.path(data_dir, "corn.mat"))
    y <- as.numeric(mat_data$propvals$data[, 1])
    X <- as.matrix(mat_data$m5spec$data)
  } else if (ds_name == "diesel") {
    df <- read.csv(file.path(data_dir, "Filtered_Diesel_Spec_Data.csv"), header = TRUE)
    y  <- as.numeric(df[,1])
    X  <- as.matrix(df[,-1])
  } else if (ds_name == "meat") {
    df <- read_excel(file.path(data_dir, "Tecator.xlsx"))
    y  <- as.numeric(df[[1]])
    X  <- as.matrix(df[, -1])
  } else if (ds_name == "milk") {
    mat_data <- readMat(file.path(data_dir, "milk.mat"))
    y  <- as.numeric(mat_data$protein)
    X  <- as.matrix(mat_data$X)
  } else if (ds_name == "soil") {
    df <- readMat(file.path(data_dir, "NIRsoil.mat"))
    y  <- as.numeric(df$soilref$data[,1])
    X  <- as.matrix(df$soil$data)
  } else if (ds_name == "tablet") {
    mat_data <- readMat(file.path(data_dir, "NIRdata_tablets.mat"))
    y  <- as.numeric(mat_data$Matrix[, 1])
    X  <- as.matrix(mat_data$Matrix[, 4:407])
  }
  return(list(X = X, y = y))
}

# ================= 3. Core Sensitivity Analysis Function =================
run_sensitivity_analysis <- function(ds_name, k = 0.7, seed = 42) {
  cat(sprintf("\n=======================================================\n"))
  cat(sprintf(">>> Generating sensitivity analysis for [%s] (k=%.1f, seed=%d)...\n", toupper(ds_name), k, seed))
  cat(sprintf("=======================================================\n"))
  
  # ---------- 3.1 Load & Partition Data ----------
  data <- load_dataset(ds_name)
  X <- data$X
  y <- data$y
  p <- ncol(X)
  
  set.seed(seed)
  n_cal_pool <- floor(0.80 * nrow(X))
  cal_pool_ix <- sample(seq_len(nrow(X)), size = n_cal_pool)
  global_test_ix  <- setdiff(seq_len(nrow(X)), cal_pool_ix)
  
  X_cal_pool <- X[cal_pool_ix, , drop = FALSE]
  n_train <- floor(k * nrow(X_cal_pool))
  ks_sub <- kenStone(X_cal_pool, k = n_train, metric = "euclid")
  train_ix <- cal_pool_ix[ks_sub$model]
  test_ix <- global_test_ix
  
  # ---------- 3.2 Preprocessing (MSC + Mean Centering) ----------
  ref      <- colMeans(X[train_ix, , drop = FALSE])
  Xtr_msc  <- msc(X[train_ix, , drop = FALSE], ref = ref)
  Xte_msc  <- msc(X[test_ix,  , drop = FALSE], ref = ref)
  
  mu_x   <- colMeans(Xtr_msc)
  mu_y   <- mean(y[train_ix])
  Xtr <- scale(Xtr_msc, center = mu_x, scale = FALSE)
  Xte <- scale(Xte_msc, center = mu_x, scale = FALSE)
  ytr <- y[train_ix] - mu_y
  yte <- y[test_ix]  - mu_y
  n_tr <- length(ytr)
  
  # ---------- 3.3 Fused Lasso Path & Mallows' Cp ----------
  D <- getD1dSparse(p)
  fuse_res <- fusedlasso(ytr, Xtr, D, gamma = 0) 
  
  # Estimate sigma^2 using full-spectrum PLS with maximum allowed LVs
  pls_ref <- plsr(ytr ~ Xtr, ncomp = min(MAX_LVS, n_tr - 1), validation = "CV")
  sigma2 <- (min(RMSEP(pls_ref, estimate = "CV")$val[1, 1, -1]))^2
  
  cp_vals <- sapply(seq_along(fuse_res$lambda), function(i) {
    beta <- coef(fuse_res, lambda = fuse_res$lambda[i])$beta
    rss  <- sum((ytr - Xtr %*% beta)^2)
    rss - n_tr * sigma2 + 2 * sigma2 * fuse_res$df[i]
  })
  
  opt_lambda_cp <- fuse_res$lambda[which.min(cp_vals)]
  lambda_vals <- fuse_res$lambda
  df_vals <- fuse_res$df
  
  # ---------- 3.4 Lambda Grid & Parallel Environment Setup ----------
  n_grid <- min(50, length(lambda_vals))
  lambda_idx_seq <- unique(pmax(1, as.integer(round(seq(1, length(lambda_vals), length.out = n_grid)))))
  inner_segments <- cvsegments(n_tr, k = 10, type = "random")
  
  cv_rmsep_once <- function(ivals) {
    if (length(ivals) == 0) return(Inf)
    cols <- sort(unlist(lapply(ivals, function(iv) seq(iv[1], iv[2]))))
    X_sub <- Xtr[, cols, drop = FALSE]
    ncomp_safe <- min(ncol(X_sub), n_tr - 1, MAX_LVS)
    if (ncomp_safe < 1) return(Inf)
    
    fit <- plsr(ytr ~ X_sub, ncomp = ncomp_safe, validation = "CV", segments = inner_segments)
    nc <- selectNcomp(fit, method = "onesigma", plot = FALSE)
    if (nc == 0) nc <- which.min(RMSEP(fit, estimate = "CV")$val[1, 1, -1])
    return(RMSEP(fit, estimate = "CV")$val[1, 1, nc + 1])
  }
  
  n_cores <- max(1, parallel::detectCores() - 1)
  cl <- makeCluster(n_cores)
  clusterEvalQ(cl, {
    library(pls)
    library(genlasso)
  })
  
  clusterExport(cl, varlist = c(
    "fuse_res", "lambda_vals", "lambda_idx_seq", "p", 
    "Xtr", "ytr", "inner_segments", "MAX_LVS", "n_tr", "cv_rmsep_once"
  ), envir = environment())
  
  cat(sprintf("Evaluating %d Lambda nodes using %d parallel cores...\n", length(lambda_idx_seq), n_cores))
  
  # ---------- 3.5 Parallel BiPLS with Progress Bar ----------
  res_list <- pblapply(seq_along(lambda_idx_seq), function(kidx) {
    idx      <- lambda_idx_seq[kidx]
    lambda_k <- lambda_vals[idx]
    
    beta_k <- coef(fuse_res, lambda = lambda_k)$beta
    cuts_k <- which(abs(diff(beta_k)) > 1e-6)
    starts_k <- c(1, cuts_k + 1)
    ends_k   <- c(cuts_k, p)
    intervals_k <- Map(function(a,b) c(a,b), starts_k, ends_k)
    
    # Merge extremely narrow intervals (width < 3)
    min_width <- 3
    if (length(intervals_k) > 1) {
      i <- 1
      while (i < length(intervals_k)) {
        width <- intervals_k[[i]][2] - intervals_k[[i]][1] + 1
        if (width < min_width) {
          intervals_k[[i+1]][1] <- intervals_k[[i]][1]
          intervals_k[[i]] <- NULL
        } else {
          i <- i + 1
        }
      }
      last_idx <- length(intervals_k)
      if (last_idx > 1 && (intervals_k[[last_idx]][2] - intervals_k[[last_idx]][1] + 1) < min_width) {
        intervals_k[[last_idx-1]][2] <- intervals_k[[last_idx]][2]
        intervals_k[[last_idx]] <- NULL
      }
    }
    
    segment_count <- length(intervals_k)
    selected_k    <- intervals_k
    best_rmsep_k  <- cv_rmsep_once(selected_k)
    
    # BiPLS Backward Elimination
    while (length(selected_k) > 1) {
      cand_rmse_k <- sapply(seq_along(selected_k), function(i) cv_rmsep_once(selected_k[-i]))
      idx_min_k <- which.min(cand_rmse_k)
      if (cand_rmse_k[idx_min_k] < best_rmsep_k) {
        best_rmsep_k <- cand_rmse_k[idx_min_k]
        selected_k   <- selected_k[-idx_min_k]
      } else {
        break
      }
    }
    
    cols_sel_k <- sort(unlist(lapply(selected_k, function(iv) seq(iv[1], iv[2]))))
    
    return(list(
      rmsep          = best_rmsep_k,
      interval_count = length(selected_k),
      variable_count = length(cols_sel_k),
      segment_count  = segment_count
    ))
  }, cl = cl)
  
  stopCluster(cl) 
  
  # ---------- 3.6 Extract Results ----------
  rmsep_list          <- sapply(res_list, function(x) x$rmsep)
  interval_count_list <- sapply(res_list, function(x) x$interval_count)
  variable_count_list <- sapply(res_list, function(x) x$variable_count)
  segment_count_list  <- sapply(res_list, function(x) x$segment_count)
  
  # ---------- 3.7 Prepare Plot Data ----------
  range_rmsecv   <- range(rmsep_list, na.rm = TRUE)
  range_segcount <- range(segment_count_list, na.rm = TRUE)
  RMSECV_scaled  <- rescale(rmsep_list, to = range_segcount) 
  
  df_plot <- data.frame(
    log_lambda    = log(lambda_vals[lambda_idx_seq]),
    RMSECV        = rmsep_list,           
    RMSECV_scaled = RMSECV_scaled,        
    Segment_Count = segment_count_list
  )
  
  log_opt_lambda_cp <- log(opt_lambda_cp)
  opt_k_idx <- which.min(abs(lambda_vals[lambda_idx_seq] - opt_lambda_cp))
  
  # ---------- 3.8 Generate Publication-Ready Plot ----------
  p_plot <- ggplot(df_plot, aes(x = log_lambda)) +
    geom_line(aes(y = Segment_Count, color = "Number of Intervals"), linetype = "dashed", size = 0.8) +
    geom_point(aes(y = Segment_Count, color = "Number of Intervals", shape = "Number of Intervals"), size = 2) +
    
    geom_line(aes(y = RMSECV_scaled, color = "RMSECV"), size = 0.8) +
    geom_point(aes(y = RMSECV_scaled, color = "RMSECV", shape = "RMSECV"), size = 2) +
    
    annotate("point", x = log_opt_lambda_cp, y = RMSECV_scaled[opt_k_idx], color = "red", shape = 8, size = 3) +
    annotate("text", x = log_opt_lambda_cp, y = RMSECV_scaled[opt_k_idx], 
             label = "optimal value", vjust = -1, hjust = 0.5, size = 3.5, color = "red", family = "Times New Roman") +
    geom_vline(xintercept = log_opt_lambda_cp, linetype = "dotted", color = "gray50", size = 0.8) +
    
    scale_color_manual(name = NULL, values = c("Number of Intervals" = "#4C72B0", "RMSECV" = "#DD8452")) +
    scale_shape_manual(name = NULL, values = c("Number of Intervals" = 17, "RMSECV" = 16)) +
    labs(x = expression(log(lambda)), y = "Number of Intervals") +
    scale_y_continuous(
      sec.axis = sec_axis(~ rescale(., from = range_segcount, to = range_rmsecv), name = "RMSECV (after BiPLS)") 
    ) +
    theme_classic(base_size = 12, base_family = "Times New Roman") +
    theme(
      axis.title.x = element_text(size = 13, color = "black", margin = margin(t = 5)),
      axis.title.y = element_text(size = 13, color = "black", margin = margin(r = 5)),
      axis.title.y.right = element_text(size = 13, color = "black", margin = margin(l = 5)),
      axis.text    = element_text(size = 11, color = "black"),
      legend.position = c(0.75, 0.85),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.text = element_text(size = 11, color = "black", family = "Times New Roman"),
      panel.border = element_rect(color = "black", fill = NA, size = 0.8),
      axis.line = element_blank() 
    )
  
  # ---------- 3.9 Save Output ----------
  pdf_path <- file.path(out_dir, paste0("lambda_sensitivity_", ds_name, ".pdf"))
  ggsave(filename = pdf_path, plot = p_plot, device = cairo_pdf, width = 7, height = 5, dpi = 600)
  
  df_csv <- data.frame(
    Dataset             = ds_name,
    Lambda              = signif(lambda_vals[lambda_idx_seq], 6),
    Log_Lambda          = round(log(lambda_vals[lambda_idx_seq]), 6),
    RMSECV              = round(rmsep_list, 6),
    Intervals_Selected  = interval_count_list,
    Segments_Fused      = segment_count_list,
    Variables_Selected  = variable_count_list,
    Cp                  = round(cp_vals[lambda_idx_seq], 6),
    Is_Optimal_Cp       = ifelse(seq_len(n_grid) == opt_k_idx, "YES", "")
  )
  csv_path <- file.path(out_dir, paste0("lambda_rmsecv_intervals_", ds_name, ".csv"))
  write.csv(df_csv, csv_path, row.names = FALSE)
  
  cat(sprintf("Output saved to: %s\n", out_dir))
}

# ================= 4. Batch Execution =================
for (ds in dataset_names) {
  run_sensitivity_analysis(ds_name = ds, k = 0.7, seed = 42)
}

cat("\nProcess completed successfully.\n")