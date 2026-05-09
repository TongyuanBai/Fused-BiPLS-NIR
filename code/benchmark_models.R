library(parallel)
library(genlasso)
library(prospectr)
library(pls)
library(R.matlab)
library(ggplot2)
library(scales)
library(readxl)
library(glmnet)

# ======= 1. Global Parameters and Paths =======
seeds <- 1:30    # 30 independent random splits for robust evaluation
k_list <- c(0.5, 0.7, 1.0)
MAX_LVS <- 20    # Maximum number of latent variables for full-spectrum PLS/PCR/BiPLS

# Relative path for outputs (ensures reproducibility across different machines)
out_dir <- "./results"
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

dataset_names <- c("corn", "diesel", "meat", "milk", "soil", "tablet")

# ======= 2. Data Loading Function =======
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

# ======= 3. Single Execution Function (Evaluates all methods) =======
run_once <- function(seed, k, X, y) {
  set.seed(seed)
  
  # ---------- 3.1 Random Global Test Set Partitioning (80/20) ----------
  n_cal_pool <- floor(0.80 * nrow(X))
  cal_pool_ix <- sample(seq_len(nrow(X)), size = n_cal_pool)
  global_test_ix  <- setdiff(seq_len(nrow(X)), cal_pool_ix)
  
  # ---------- 3.2 Calibration Sampling via Kennard-Stone ----------
  X_cal_pool <- X[cal_pool_ix, , drop = FALSE]
  n_train <- floor(k * nrow(X_cal_pool))
  
  if (n_train >= nrow(X_cal_pool)) {
    train_ix <- cal_pool_ix
  } else {
    ks_sub <- kenStone(X_cal_pool, k = n_train, metric = "euclid")
    train_ix <- cal_pool_ix[ks_sub$model]
  }
  test_ix <- global_test_ix
  
  # ---------- 3.3 Preprocessing: MSC + Mean Centering ----------
  # Parameters derived exclusively from the training set
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
  p <- ncol(Xtr)
  
  calc_metrics <- function(y_true, y_pred) {
    rmse <- sqrt(mean((y_true - y_pred)^2))
    r2 <- 1 - sum((y_true - y_pred)^2) / sum((y_true - mean(y_true))^2)
    return(c(RMSE = rmse, R2 = r2))
  }
  
  results <- list() 
  
  # ========================================================
  # METHOD 1: Fused Lasso + BiPLS (Proposed)
  # ========================================================
  t0 <- Sys.time()
  
  # Stage 1: Structural Spectral Segmentation via Fused Lasso
  D <- getD1dSparse(p)
  fuse_res <- fusedlasso(ytr, Xtr, D, gamma = 0) 
  
  # Estimate noise variance (sigma^2) using cross-validated PLS
  pls_ref <- plsr(ytr ~ Xtr, ncomp = min(MAX_LVS, n_tr - 1), validation = "CV")
  sigma2 <- (min(RMSEP(pls_ref, estimate = "CV")$val[1, 1, -1]))^2
  
  # Optimal lambda selection via Mallows' Cp
  cp_vals <- sapply(seq_along(fuse_res$lambda), function(i) {
    beta <- coef(fuse_res, lambda = fuse_res$lambda[i])$beta
    rss  <- sum((ytr - Xtr %*% beta)^2)
    rss - n_tr * sigma2 + 2 * sigma2 * fuse_res$df[i]
  })
  opt_lambda <- fuse_res$lambda[which.min(cp_vals)]
  beta_opt <- coef(fuse_res, lambda = opt_lambda)$beta
  
  # Identify piecewise segments
  cuts   <- which(abs(diff(beta_opt)) > 1e-6)
  starts <- c(1, cuts + 1)
  ends   <- c(cuts, p)
  intervals <- Map(function(a, b) c(a, b), starts, ends)
  
  # Merge extremely narrow intervals (width < 3)
  min_width <- 3
  if (length(intervals) > 1) {
    i <- 1
    while (i < length(intervals)) {
      width <- intervals[[i]][2] - intervals[[i]][1] + 1
      if (width < min_width) {
        intervals[[i+1]][1] <- intervals[[i]][1]
        intervals[[i]] <- NULL
      } else {
        i <- i + 1
      }
    }
    last_idx <- length(intervals)
    if (last_idx > 1) {
      if ((intervals[[last_idx]][2] - intervals[[last_idx]][1] + 1) < min_width) {
        intervals[[last_idx-1]][2] <- intervals[[last_idx]][2]
        intervals[[last_idx]] <- NULL
      }
    }
  }
  
  # Stage 2: Backward Interval Elimination via PLS
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
  
  curr_ivs <- intervals
  best_ivs <- curr_ivs
  best_cv_rmse <- cv_rmsep_once(curr_ivs)
  
  while (length(curr_ivs) > 1) {
    rmses <- sapply(seq_along(curr_ivs), function(i) cv_rmsep_once(curr_ivs[-i]))
    idx_drop <- which.min(rmses)
    
    if (rmses[idx_drop] < best_cv_rmse) {
      best_cv_rmse <- rmses[idx_drop]
      curr_ivs <- curr_ivs[-idx_drop]
      best_ivs <- curr_ivs
    } else { 
      break 
    }
  }
  
  # Final model calibration and prediction
  cols_final <- sort(unlist(lapply(best_ivs, function(iv) seq(iv[1], iv[2]))))
  final_fit  <- plsr(ytr ~ Xtr[, cols_final, drop = FALSE], ncomp = min(length(cols_final), n_tr - 1, MAX_LVS), validation = "CV")
  
  final_nc <- selectNcomp(final_fit, method = "onesigma", plot = FALSE)
  if (final_nc == 0) final_nc <- which.min(RMSEP(final_fit, estimate = "CV")$val[1, 1, -1])
  
  pred_te_ours <- predict(final_fit, newdata = Xte[, cols_final, drop = FALSE], ncomp = final_nc)[, , 1]
  mets_ours <- calc_metrics(yte, pred_te_ours)
  time_ours <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  
  # Format intervals for logging
  str_init_intervals <- paste(sapply(intervals, function(iv) paste0(iv[1], "-", iv[2])), collapse = "; ")
  str_best_intervals <- paste(sapply(best_ivs, function(iv) paste0(iv[1], "-", iv[2])), collapse = "; ")
  
  results[[1]] <- data.frame(
    Method = "Ours (Fused-biPLS)", 
    seed = seed, 
    k_prop = k, 
    RMSE_Test = mets_ours["RMSE"], 
    R2_Test = mets_ours["R2"], 
    Model_Info = paste0("Vars:", length(cols_final), "; Intvs:", length(best_ivs)), 
    Time = time_ours,
    Opt_Lambda = opt_lambda,                      
    Initial_Intervals = str_init_intervals,       
    Selected_Intervals = str_best_intervals,      
    Selected_Vars = paste(cols_final, collapse = ";")
  )
  
  # ========================================================
  # METHOD 2: PLS
  # ========================================================
  t0 <- Sys.time()
  fit_pls <- plsr(ytr ~ Xtr, ncomp = min(MAX_LVS, n_tr-1), validation = "CV", segments=10)
  nc_pls <- selectNcomp(fit_pls, method = "onesigma", plot = FALSE)
  if (nc_pls == 0) nc_pls <- which.min(RMSEP(fit_pls, estimate = "CV")$val[1,1,-1])
  pred_pls <- predict(fit_pls, newdata = Xte, ncomp = nc_pls)[,,1]
  mets_pls <- calc_metrics(yte, pred_pls)
  
  results[[2]] <- data.frame(Method = "PLS", seed = seed, k_prop = k, RMSE_Test = mets_pls["RMSE"], 
                             R2_Test = mets_pls["R2"], Model_Info = paste0(nc_pls, "LV"), 
                             Time = as.numeric(difftime(Sys.time(), t0, units="secs")),
                             Opt_Lambda = NA, Initial_Intervals = NA, Selected_Intervals = NA, 
                             Selected_Vars = "All")
  
  # ========================================================
  # METHOD 3: PCR
  # ========================================================
  t0 <- Sys.time()
  pca_res <- prcomp(Xtr, center = FALSE, scale. = FALSE)
  max_pc <- min(MAX_LVS, n_tr-1)
  best_pc <- 1; min_rmse_pc <- Inf
  folds <- sample(rep(1:10, length.out = n_tr))
  
  for (pc in 1:max_pc) {
    preds_cv <- numeric(n_tr)
    for (f in 1:10) {
      idx_val <- which(folds == f); idx_tr <- which(folds != f)
      X_c <- pca_res$x[idx_tr, 1:pc, drop = FALSE]
      fit_lm <- lm(ytr[idx_tr] ~ X_c - 1)
      preds_cv[idx_val] <- pca_res$x[idx_val, 1:pc, drop = FALSE] %*% coef(fit_lm)
    }
    rmse_cv <- sqrt(mean((ytr - preds_cv)^2))
    if (rmse_cv < min_rmse_pc) { min_rmse_pc <- rmse_cv; best_pc <- pc }
  }
  
  fit_lm_final <- lm(ytr ~ pca_res$x[, 1:best_pc, drop=FALSE] - 1)
  Xte_pca <- predict(pca_res, newdata = Xte)
  pred_pca <- Xte_pca[, 1:best_pc, drop=FALSE] %*% coef(fit_lm_final)
  mets_pca <- calc_metrics(yte, pred_pca)
  
  results[[3]] <- data.frame(Method = "PCR", seed = seed, k_prop = k, RMSE_Test = mets_pca["RMSE"], 
                             R2_Test = mets_pca["R2"], Model_Info = paste0(best_pc, "PC"), 
                             Time = as.numeric(difftime(Sys.time(), t0, units="secs")),
                             Opt_Lambda = NA, Initial_Intervals = NA, Selected_Intervals = NA, 
                             Selected_Vars = "All")
  
  # ========================================================
  # METHOD 4: LASSO
  # ========================================================
  t0 <- Sys.time()
  lasso_fit <- cv.glmnet(Xtr, ytr, alpha = 1, nfolds = 10, intercept = FALSE)
  pred_lasso <- as.numeric(predict(lasso_fit, newx = Xte, s = "lambda.min"))
  
  lasso_coef <- as.numeric(coef(lasso_fit, s = "lambda.min"))[-1]
  lasso_idx <- which(lasso_coef != 0)
  lasso_vars_str <- if(length(lasso_idx) > 0) paste(lasso_idx, collapse = ";") else "None"
  
  mets_lasso <- calc_metrics(yte, pred_lasso)
  results[[4]] <- data.frame(Method = "LASSO", seed = seed, k_prop = k, RMSE_Test = mets_lasso["RMSE"], 
                             R2_Test = mets_lasso["R2"], Model_Info = paste0(length(lasso_idx), "vars"), 
                             Time = as.numeric(difftime(Sys.time(), t0, units="secs")),
                             Opt_Lambda = NA, Initial_Intervals = NA, Selected_Intervals = NA, 
                             Selected_Vars = lasso_vars_str)
  
  # ========================================================
  # METHOD 5: MWPLS
  # ========================================================
  t0 <- Sys.time()
  win_size <- 15
  best_s <- 1; min_rmse_mw <- Inf; best_nc_mw <- 1
  
  for(s in 1:(p - win_size + 1)) {
    e <- s + win_size - 1
    Xw <- Xtr[, s:e, drop = FALSE]
    max_c <- min(ncol(Xw), n_tr-1, 10) # Max 10 LVs for window
    fit_w <- plsr(ytr ~ Xw, ncomp = max_c, validation = "CV", segments = 10)
    rmses <- RMSEP(fit_w, estimate = "CV")$val[1,1,-1]
    nc_w <- which.min(rmses)
    if(rmses[nc_w] < min_rmse_mw) {
      min_rmse_mw <- rmses[nc_w]; best_s <- s; best_nc_mw <- nc_w
    }
  }
  
  best_e <- best_s + win_size - 1
  fit_mw_final <- plsr(ytr ~ Xtr[, best_s:best_e, drop = FALSE], ncomp = best_nc_mw)
  pred_mw <- predict(fit_mw_final, ncomp = best_nc_mw, newdata = Xte[, best_s:best_e, drop = FALSE])[,,1]
  mets_mw <- calc_metrics(yte, pred_mw)
  
  results[[5]] <- data.frame(Method = "MWPLS", seed = seed, k_prop = k, RMSE_Test = mets_mw["RMSE"], 
                             R2_Test = mets_mw["R2"], Model_Info = paste0("Win", best_s, "-", best_e), 
                             Time = as.numeric(difftime(Sys.time(), t0, units="secs")),
                             Opt_Lambda = NA, Initial_Intervals = NA, Selected_Intervals = NA, 
                             Selected_Vars = paste(best_s:best_e, collapse = ";"))
  
  # ========================================================
  # METHOD 6 & 7: MC-UVE & SCARS
  # ========================================================
  mc_nrep <- 50
  B_mat <- matrix(0, nrow = mc_nrep, ncol = p)
  for(i in 1:mc_nrep) {
    idx <- sample(1:n_tr, size = floor(0.8 * n_tr))
    fit_b <- plsr(ytr[idx] ~ Xtr[idx, , drop=FALSE], ncomp = min(10, length(idx)-1))
    B_mat[i,] <- as.vector(coef(fit_b, ncomp = fit_b$ncomp, intercept = FALSE))
  }
  
  vi_uve <- abs(colMeans(B_mat)) / (apply(B_mat, 2, sd) + 1e-8)
  vi_scars <- - (apply(B_mat, 2, sd) / (colMeans(abs(B_mat)) + 1e-8))
  
  eval_var_select <- function(vi_scores, method_name) {
    t_start <- Sys.time()
    sorted_idx <- order(vi_scores, decreasing = TRUE)
    best_k <- 5; min_rmse_sel <- Inf
    
    for(k_vars in seq(5, min(length(sorted_idx), 100), by = 5)) {
      sel <- sorted_idx[1:k_vars]
      max_lv_sel <- min(10, k_vars, n_tr-1) # Max 10 LVs for subsets
      fit_cv <- plsr(ytr ~ Xtr[, sel, drop=FALSE], ncomp = max_lv_sel, validation = "CV", segments = 10)
      rmse_cv <- min(RMSEP(fit_cv, estimate = "CV")$val[1,1,-1])
      if(rmse_cv < min_rmse_sel) { min_rmse_sel <- rmse_cv; best_k <- k_vars }
    }
    
    sel_final <- sorted_idx[1:best_k]
    fit_final <- plsr(ytr ~ Xtr[, sel_final, drop=FALSE], ncomp = min(10, best_k, n_tr-1))
    pred_sel <- predict(fit_final, ncomp = fit_final$ncomp, newdata = Xte[, sel_final, drop=FALSE])[,,1]
    mets_sel <- calc_metrics(yte, pred_sel)
    
    return(data.frame(Method = method_name, seed = seed, k_prop = k, 
                      RMSE_Test = mets_sel["RMSE"], R2_Test = mets_sel["R2"], 
                      Model_Info = paste0(best_k, "vars"), 
                      Time = as.numeric(difftime(Sys.time(), t_start, units="secs")),
                      Opt_Lambda = NA, Initial_Intervals = NA, Selected_Intervals = NA, 
                      Selected_Vars = paste(sel_final, collapse = ";"))) 
  }
  
  results[[6]] <- eval_var_select(vi_uve, "MC-UVE")
  results[[7]] <- eval_var_select(vi_scars, "SCARS")
  
  # Status Output
  cat(sprintf("Seed:%d | k:%.1f | %s (Ours) RMSE:%.4f | PLS RMSE:%.4f\n", 
              seed, k, "Fused-BiPLS", mets_ours["RMSE"], mets_pls["RMSE"]))
  
  return(do.call(rbind, results))
}

# ======= 4. Main Execution Loop =======

for (ds in dataset_names) {
  cat("\n=======================================================\n")
  cat(sprintf("Processing Dataset: %s", toupper(ds)), "\n")
  cat("=======================================================\n")
  
  data <- load_dataset(ds)
  X_current <- data$X
  y_current <- data$y
  
  grid <- expand.grid(seed = seeds, k = k_list)
  cl <- makeCluster(max(1, parallel::detectCores() - 1), outfile = "")
  
  clusterEvalQ(cl, { 
    library(genlasso); library(prospectr); library(pls); library(glmnet) 
  })
  clusterExport(cl, varlist = c("grid", "run_once", "MAX_LVS", "X_current", "y_current"), envir = environment())
  
  res_list <- parLapply(cl, 1:nrow(grid), function(i) {
    run_once(grid$seed[i], grid$k[i], X_current, y_current)
  })
  stopCluster(cl)
  
  results_all <- do.call(rbind, res_list)
  
  agg_results <- do.call(data.frame, aggregate(
    cbind(RMSE_Test, R2_Test) ~ k_prop + Method, 
    data = results_all, 
    FUN = function(x) c(mean = mean(x), sd = sd(x))
  ))
  
  print(paste("Aggregated Results for", ds, ":"))
  print(head(agg_results, 10))
  
  csv_detail_path <- file.path(out_dir, paste0(ds, "_All_Methods_Detailed.csv"))
  csv_agg_path    <- file.path(out_dir, paste0(ds, "_All_Methods_Summary.csv"))
  plot_path       <- file.path(out_dir, paste0(ds, "_Learning_Curve.png"))
  
  write.csv(results_all, csv_detail_path, row.names = FALSE)
  write.csv(agg_results, csv_agg_path, row.names = FALSE)
  
  p <- ggplot(agg_results, aes(x = k_prop, y = RMSE_Test.mean, color = Method, group = Method)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    scale_x_continuous(breaks = k_list) +
    labs(title = paste("RMSEP Learning Curve -", toupper(ds)),
         subtitle = "Comparison of Baseline Methods vs. Ours",
         x = "Calibration Sampling Ratio (k)", 
         y = "RMSEP on Random Global Test Set") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave(plot_path, plot = p, width = 8, height = 6, dpi = 300)
  
  cat(sprintf("Dataset %s completed. Results saved to %s\n", toupper(ds), out_dir))
}

cat("\nAll 6 datasets processed successfully.\n")