# ==============================================================================
# GreedyComBat Covariate Selection Functions
# ==============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

#' Calculate Cross-Validated Fraction of Variance Explained (CV-FVE)
#'
#' Evaluates the amount of variance in a set of features collectively explained by a set of covariates using cross-validation. 
#' Gracefully handles rank-deficiency and perfectly collinear covariates using QR decomposition.
#'
#' @param data A data frame containing features and covariates. All categorical 
#' variables MUST be pre-converted to factors. Optionally includes a `subid` column 
#' for subject-level fold blocking (repeated measures).
#' @param feat_cols A character vector specifying the names of the continuous feature columns (e.g., brain ROIs).
#' @param covariates A character vector specifying the names of the covariate columns to include in the model.
#' @param seed An integer used to set the random seed for reproducible cross-validation fold generation. Default is 42.
#' @param cv_folds An integer specifying the number of cross-validation folds. Default is 5.
#' @param min_n An integer specifying the minimum number of rows required in `data`. Default is 20.
#'
#' @return A list containing two numeric values: \code{R2} (the training fraction of variance explained, in percent) 
#' and \code{R2_cv} (the cross-validated fraction of variance explained, in percent).
calc_cv_fve <- function(data, feat_cols, covariates, seed = 42, cv_folds = 5, min_n = 20) {
  
  # Determine the columns that MUST exist i.e. the covariates and features
  all_cols <- c(feat_cols, covariates)
  
  # Immediately error if any required column does not exist
  if (!all(all_cols %in% names(data))) {
    stop(paste("ERROR: Columns missing from data:", paste(setdiff(all_cols, names(data)), collapse=", ")))
  }
  
  # Immediately error if there are ANY NAs in either the covariate or feature columns
  if (any(sapply(data[, all_cols, drop=FALSE], function(x) any(!is.finite(x))))) {
    stop("ERROR: Missing (NA) or non-finite (Inf/NaN) data detected.")
  }
  
  # Check sample size against min
  if (nrow(data) < min_n) {
    stop(sprintf("ERROR: Dataset has too few rows (%d). Minimum required: %d.", nrow(data), min_n))
  }
  
  # Check sample size against number of folds
  if (nrow(data) < cv_folds) {
    stop(sprintf("ERROR: Cannot perform %d-fold CV on only %d rows.", cv_folds, nrow(data)))
  }

  # Default return value if no covariates remain
  if (length(covariates) == 0) return(list(R2 = 0, R2_cv = 0))
  
  # Define the model formula and numeric matrix
  fmla <- as.formula(paste0("~ 1 + ", paste(covariates, collapse = " + ")))
  mm <- model.matrix(fmla, data = data)
  
  # Extract the subid if it exists, so we can subject block for CV folds (if necessary)
  # This is necessary for selection for repeated measures harmonisation e.g. LongComBat
  grp <- if ("subid" %in% names(data)) as.character(data$subid) else NULL
  
  set.seed(seed + 1001) 
  
  # Define the folds
  if (!is.null(grp)) {
    
    # Subject IDs present
    u <- unique(grp)
    
    # Define the lookup table for subid to fold
    fold_map <- setNames(sample(rep(seq_len(cv_folds), length.out = length(u))), u)
    
    # Apply the lookup table to determine the folds for each scan
    fid <- as.integer(fold_map[grp])
  } else {
    
    # Subject IDs not present, randomly assign folds to each scan
    fid <- sample(rep(seq_len(cv_folds), length.out = nrow(data)))
  }
  
  # Running total for SSE and SST
  global_sse <- 0; global_sst <- 0
  
  # Training
  sst_tr <- 0; ssr_tr <- 0
  
  # Iterate through features
  for (f in feat_cols) {
    
    # Extract data corresponding to each feature
    y <- data[[f]]
    
    # Predicted values
    yhat <- rep(NA, length(y))
    
    # Iterate through folds
    for (fold in sort(unique(fid))) {
      
      # Identify the test and train folds
      te <- which(fid == fold)
      tr <- which(fid != fold)
      
      # QR decomposition to calculate coefficients
      # Also deals with collinearity (column dropping) and 0 variance columns
      b <- qr.coef(qr(mm[tr, , drop=FALSE]), y[tr]); b[!is.finite(b)] <- 0
      
      # Apply those coefficients to predict thicknesses for the test data
      yhat[te] <- drop(mm[te, , drop=FALSE] %*% b)
    }
    
    # Update the SSE and SST totals
    global_sse <- global_sse + sum((y - yhat)^2)
    global_sst <- global_sst + sum((y - mean(y))^2)
    
    # Training values
    b_tr <- qr.coef(qr(mm), y); b_tr[!is.finite(b_tr)] <- 0
    sst_tr <- sst_tr + sum((y - mean(y))^2)
    ssr_tr <- ssr_tr + sum((drop(mm %*% b_tr) - mean(y))^2)
  }
  
  # Calculate the CV FVE values
  r2_cv <- if (global_sst > 0) 100 * (1 - global_sse / global_sst) else 0
  r2_tr <- if(sst_tr > 0) 100 * ssr_tr / sst_tr else 0
  
  # Return the training and CV FVE values
  return(list(R2 = r2_tr, R2_cv = r2_cv))
}

#' Select Covariates using a Greedy CV-FVE Algorithm
#'
#' Performs variable selection by first applying a marginal CV-FVE screen, followed by a greedy forward 
#' selection algorithm to iteratively build a covariate set that maximizes out-of-sample explained variance.
#'
#' @param df_for_cv A data frame containing the features and candidate covariates.
#' @param feat_cols A character vector specifying the names of the feature columns.
#' @param candidates A character vector specifying the names of the candidate covariate columns to evaluate.
#' @param seed_base An integer used as the base random seed for reproducibility across iterations. Default is 42.
#' @param keep_threshold A numeric value specifying the minimum CV-FVE (in percent) required for a covariate 
#' to pass the initial marginal screening. Default is 0.5.
#' @param min_gain A numeric value specifying the minimum increase in CV-FVE (in percent) required to add 
#' a covariate during the greedy forward selection phase. Default is 0.5.
#' @param max_k An integer specifying the maximum number of covariates the algorithm is permitted to select. Default is 30.
#'
#' @return A character vector of the selected covariate names.
select_greedy_vars <- function(df_for_cv, feat_cols, candidates,
                               seed_base = 42, keep_threshold = 0.5,
                               min_gain = 0.5, max_k = 30) {
  
  # If no candidates provided, stop
  if (!length(candidates)) return(character(0))
  
  # Run the CV-FVE algorithm for each candidate covariate individually
  screen_df <- bind_rows(lapply(candidates, function(v) {
    res <- calc_cv_fve(df_for_cv, feat_cols, covariates = c(v), seed = seed_base)
    tibble(variable = v, cv = res$R2_cv)
  })) %>% arrange(desc(cv))
  
  # Keep the variables with marginal CV-FVE above the min threshold
  keep_vars <- screen_df %>% filter(is.finite(cv), cv >= keep_threshold) %>% pull(variable)
  if (!length(keep_vars)) return(character(0))
  
  selected <- character(0)
  current_cv <- 0
  
  # Loop through each variable addition consideration
  for (k in seq_len(max_k)) {
    if (!length(keep_vars)) break
    iter_seed <- seed_base + (1000 * k)
    
    # Check the CV-FVE scores after addition of each allowable covariate
    scores <- bind_rows(lapply(keep_vars, function(v) {
      res <- calc_cv_fve(df_for_cv, feat_cols, covariates = c(selected, v), seed = iter_seed)
      tibble(var = v, cv = res$R2_cv)
    })) %>% arrange(desc(cv))
    
    # Select the best covariate to add
    best <- scores[1, ]
    
    # If gain less than the threshold, stop
    if ((best$cv - current_cv) < min_gain) break
    
    # Otherwise include the best covariate and continue
    selected <- c(selected, best$var)
    current_cv <- best$cv
    keep_vars <- setdiff(keep_vars, best$var)
  }
  unique(selected)
}

