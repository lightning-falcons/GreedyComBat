suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(tidyr); library(tibble); library(ComBatFamily); library(patchwork)
})

# Source the GreedyComBat Functions
source("GreedyComBat.R")

# ==============================================================================
# ADNI Analysis Specific Globals
# ==============================================================================

# Global core covariates (reuse everywhere)
CORE_VARS <- c("baselineAge", "SEX", "DIAGNOSIS")

# Slope helper
compute_slope <- function(x, y) {
  cov(x, y) / var(y)
}

compute_cohen_d <- function(x, y, level0 = NULL, level1 = NULL) {
  
  # Null check for levels
  if (is.null(level0) || is.null(level1)) {
    stop("Must provide level0 and level1")
  }
  
  # Determine the groups
  g0 <- x[y == level0]
  g1 <- x[y == level1]
  
  # Calculate pooled standard deviation
  n0 <- length(g0); n1 <- length(g1)
  sd_p <- sqrt(((n0 - 1) * var(g0) + (n1 - 1) * var(g1)) / (n0 + n1 - 2))
  
  # Return Cohen's d
  (mean(g1) - mean(g0)) / sd_p
}

# ==============================================================================
# Data Preparation
# ==============================================================================
ants_fp   <- "ADNI_ANTsSST_protocol_2019_03_14.csv"
covar_fp  <- "combat_outputs/combat_covariates_v7_FINAL_IMPUTED.csv"
schema_fp <- "adni_covariates_curated_schema_v5.csv"

# Load data, covariates and schema
# Load the data and covariates exactly as is (characters)
raw_harm <- read_csv(ants_fp, col_types = cols(.default = col_character()), show_col_types = FALSE)
raw_cov  <- read_csv(covar_fp, col_types = cols(.default = col_character()), show_col_types = FALSE)
raw_sch  <- read_csv(schema_fp, show_col_types = FALSE)

# Create the schema lookup table mapping the variables to their types
raw_sch <- raw_sch %>% 
  transmute(variable = trimws(variable), vic_type = trimws(vic_type)) %>% 
  filter(nzchar(variable))
vmap <- as.list(setNames(raw_sch$vic_type, raw_sch$variable))

# Sets default data types for baselineAge, SEX, DIAGNOSIS and Siemens
defaults <- list(baselineAge="numeric", SEX="factor", DIAGNOSIS="factor", Siemens="factor")
for(n in names(defaults)) if(is.null(vmap[[n]])) vmap[[n]] <- defaults[[n]]

# Take only the first scan for each subid
df_harm_base <- raw_harm %>% 
  mutate(EXAM_DATE = as.Date(EXAM_DATE)) %>% 
  filter(!is.na(subid), !is.na(EXAM_DATE)) %>% 
  arrange(subid, EXAM_DATE, IMAGE_ID) %>% 
  group_by(subid) %>% 
  slice(1) %>% 
  ungroup()

# Filter out duplicate scans for covariate information
df_cov_unique <- raw_cov %>% 
  distinct(IMAGE_ID, .keep_all = TRUE)

# Keep only the covariate rows corresponding to scans of interest
df_cov_valid <- df_cov_unique %>% 
  semi_join(df_harm_base, by = "IMAGE_ID")

# Get the list of valid scanners, which are those with at least 2 scans
valid_scanners <- df_cov_valid %>% 
  filter(nzchar(manufac.model.coil.strength.site)) %>% 
  count(manufac.model.coil.strength.site) %>% 
  filter(n >= 2) %>% pull(manufac.model.coil.strength.site)

# Keep only those scans from valid scanners and the corresponding covariate information
df_cov_final  <- df_cov_valid %>% filter(manufac.model.coil.strength.site %in% valid_scanners)
df_harm_final <- df_harm_base %>% semi_join(df_cov_final, by = "IMAGE_ID")

# Look for repeated column names in the covariates and data and drop those from df_harm_final
overlap_cols <- setdiff(intersect(names(df_harm_final), names(df_cov_final)), "IMAGE_ID")
if(length(overlap_cols) > 0) {
  df_harm_final <- df_harm_final %>% select(-all_of(overlap_cols))
}

# Join the data and covariates info to create the master table
df_model <- inner_join(df_harm_final, df_cov_final, by = "IMAGE_ID")

# Identify feature columns
feat_cols <- grep("^thickness", names(df_model), value = TRUE)

# Store every feature as numeric 
for(f in feat_cols) df_model[[f]] <- as.numeric(df_model[[f]])

# Drop any non-schema covariates
df_model <- df_model %>% 
  select(any_of(c("IMAGE_ID", "subid", "ID", "manufac.model.coil.strength.site", feat_cols, names(vmap))))

# Identify every column that isn't an ID or a brain feature
vars_to_type <- setdiff(names(df_model), c("IMAGE_ID", "subid", "ID", feat_cols))

# Type convert df_model columns according to the schema
for (v in vars_to_type) {
  # Warning if the type can't be found in the schema
  # Search for the appropriate type based on the earlier lookup table
  if (is.null(vmap[[v]])) {
    warning(paste("Type unknown for column:", v, "- defaulting to 'factor'"))
    type <- "factor"
  } else {
    type <- tolower(vmap[[v]])
  }
  
  # Convert the types
  if (type %in% c("numeric", "integer", "double")) {
    df_model[[v]] <- as.numeric(df_model[[v]]) 
  } else {
    df_model[[v]] <- as.factor(df_model[[v]])
  }
}

# Define the batch column name
batch_col <- "manufac.model.coil.strength.site"

df_harm   <- df_model
X_raw     <- as.matrix(df_harm[, feat_cols, drop = FALSE])
batch_vec <- df_harm[[batch_col]]

message("Global Data Ready: df_model (", nrow(df_model), " rows).")

# Look at only the columns in the schema
candidates <- intersect(names(df_model), names(vmap))

# Remove all ID related columns may have been in the schema
blacklist <- c("subid", "IMAGE_ID", "ID", batch_col, feat_cols)

# Remove Siemens
siemens_cols <- "Siemens"
blacklist <- unique(c(blacklist, siemens_cols))

# Apply the blacklist to the candidates list
candidates <- setdiff(candidates, blacklist)

# Keep only covariates which have at least 2 unique values
candidates <- candidates[sapply(df_model[candidates], function(x) length(unique(x)) > 1)]

message("Screening ", length(candidates), " candidates...")

# ==============================================================================
# FIGURE 1
# ==============================================================================

# Screen each covariate by itself and then order by CV-FVE
# For the publication table
screen_df <- bind_rows(lapply(candidates, function(v) {
  res <- calc_cv_fve(df_model, feat_cols, covariates = c(v), seed = 42)
  tibble(variable = v, cv = res$R2_cv, train = res$R2)
})) %>% arrange(desc(cv))

write_csv(screen_df, "marginal_screen_CV_FVE.csv")

# Run GreedyComBat variable selection
selected <- select_greedy_vars(df_model, feat_cols, candidates, seed_base = 42,
                               keep_threshold = 0.5, min_gain = 0.5, max_k = 30)

# Create the tibble that stores the sequence of variables added
trace <- tibble(k = 0, added = "<none>", cv_fve = 0, train_fve = 0)

# Reply the CV-FVE and training values in sequence by re-calling calc_cv_fve
if (length(selected)) {
  cur <- character(0); current_cv <- 0
  for (k in seq_along(selected)) {
    v <- selected[k]
    res <- calc_cv_fve(df_model, feat_cols, covariates = c(cur, v), seed = 42 + 1000 * k)
    cur <- c(cur, v)
    trace <- bind_rows(trace, tibble(k = k, added = v, cv_fve = res$R2_cv, train_fve = res$R2))
    current_cv <- res$R2_cv
    message(sprintf("k=%d added %s (CV=%.2f%%)", k, v, res$R2_cv))
  }
}

write_csv(trace, "greedy_CV_FVE_final.csv")

# Plot FIGURE 1
# Neater label names
df_plot <- trace %>%
  mutate(
    label = gsub("_", " ", gsub("^.*__", "", ifelse(added=="<none>", "", added))),
    label = gsub("weight kg", "Weight (kg)", label, ignore.case = TRUE),
    label = gsub("entry age", "Entry Age", label, ignore.case = TRUE)
  ) %>%
  arrange(k)

# Label positioning
lab_df <- df_plot %>% filter(k > 0)
kmax <- max(lab_df$k)
lab_df <- lab_df %>% mutate(nx = ifelse(k==kmax, -0.22, -0.15), ny = ifelse(k==kmax, 0.55, -0.45))

p1 <- ggplot(df_plot, aes(x = k, y = cv_fve)) +
  geom_line(linewidth = 1.1, colour = "#0072B2") +
  geom_point(shape = 21, size = 3.2, stroke = 1.1, colour = "#0072B2", fill = "white") +
  scale_x_continuous(breaks = seq(0, max(df_plot$k), 1)) +
  geom_text(data = lab_df, aes(label = label), nudge_x = lab_df$nx, nudge_y = lab_df$ny, size = 5) +
  labs(x = "Variables added (k)", y = "CV-FVE (%)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.major = element_line(colour = "grey90", linewidth = 0.45),
        panel.grid.minor = element_blank(), axis.title = element_text(size = 18),
        axis.text = element_text(size = 12), plot.margin = margin(8, 18, 8, 8)) +
  coord_cartesian(clip = "off")

ggsave("Figure_1.pdf", p1, width = 7.2, height = 5, device = cairo_pdf)
ggsave("Figure_1.png", p1, width = 7.2, height = 5, dpi = 600)
message("Figure 1 Generated.")

# ==============================================================================
# FIGURE 2
# ==============================================================================
# Get the greedy select variables
greedy_vars <- selected

# The covariates to use for harmonisation
cov_sets <- list(
  RAW         = NULL,
  CORE        = CORE_VARS,
  GREEDY      = greedy_vars,
  CORE_GREEDY = unique(c(CORE_VARS, greedy_vars)),
  ALL         = unique(c(CORE_VARS, candidates)) 
)

# Ensure core variables exist
if(!all(CORE_VARS %in% names(df_model))) stop("Core variables missing from df_model!")

message("Covariate Sets Defined.")

# ComBat wrapper
run_combat <- function(X, bat, cov_df) {
  
  # Empty covariates case
  if (ncol(cov_df) == 0) {
    model_matrix <- data.frame(`(Intercept)` = 1)[rep(1, nrow(X)), , drop = FALSE]
    ff <- stats::as.formula("y ~ 1")
  } else {
    model_matrix <- cov_df
    ff <- stats::reformulate(names(cov_df), response = "y")
  }
  
  # Run ComBat
  out <- ComBatFamily::comfam(
    data    = X,
    bat     = bat, 
    covar   = model_matrix,
    model   = stats::lm, 
    formula = ff
  )
  
  out$dat.combat
}


message("Starting Harmonisation...")

harm_datasets <- list()
harm_datasets[["RAW"]] <- X_raw

# Harmonise everything that's not RAW
for (set_name in setdiff(names(cov_sets), "RAW")) {
  covs <- cov_sets[[set_name]]
  
  # Check for constant variables
  current_cov_df <- df_harm[, covs, drop = FALSE]
  
  bad_vars <- names(current_cov_df)[sapply(current_cov_df, function(x) length(unique(x)) < 2)]
  if (length(bad_vars) > 0) {
    stop(paste("ERROR: Covariate set", set_name, "contains constant variables:", paste(bad_vars, collapse=", ")))
  }
  
  # Harmonise
  harm_datasets[[set_name]] <- run_combat(X_raw, batch_vec, current_cov_df)
}

message("Harmonisation Complete.")

# Brain region of interest
roi_col <- "thickness.left.superior.frontal"

# Stop if not found
if (!(roi_col %in% names(df_harm))) {
  stop(paste("CRITICAL ERROR: Target ROI", roi_col, "not found in dataset!"))
}

message("Plotting ROI: ", roi_col)

plot_data <- list()
stat_data <- list()

# Go through each dataset
for (set_name in names(harm_datasets)) {
  
  # Extract the ROI
  y_vals <- harm_datasets[[set_name]][, roi_col]
  
  # Select the core covariates columns, the outcome and the scanner
  tmp_df <- df_harm %>%
    select(baselineAge, SEX, DIAGNOSIS) %>%
    mutate(outcome = y_vals, scanner = batch_vec)
  
  # Fit the model and get the residuals
  fit <- lm(outcome ~ baselineAge + SEX + DIAGNOSIS, data = tmp_df)
  tmp_df$resid <- residuals(fit)
  
  # Test for the mean and variance differences between scanners
  aov_p <- tryCatch({
    m <- lm(resid ~ scanner, data = tmp_df)
    anova(m)["scanner", "Pr(>F)"]
  }, error = function(e) NA)
  
  flg_p <- tryCatch({
    fligner.test(resid ~ scanner, data = tmp_df)$p.value
  }, error = function(e) NA)
  
  stat_data[[set_name]] <- tibble(method = set_name, anova_p = aov_p, fligner_p = flg_p)
  
  # Order scanners by their means
  scanner_ord <- tmp_df %>% 
    group_by(scanner) %>% 
    summarise(m = mean(resid), .groups="drop") %>% 
    arrange(m) %>% 
    pull(scanner)
  
  # Store the scanner order
  tmp_df$scanner_f <- factor(tmp_df$scanner, levels = scanner_ord)
  
  # Label the method
  tmp_df$key <- factor(paste0(set_name, "__", tmp_df$scanner_f), 
                       levels = paste0(set_name, "__", scanner_ord))
  tmp_df$method <- set_name
  
  # Select only the necessary information for plotting
  plot_data[[set_name]] <- tmp_df %>% select(method, key, resid)
}

# Combien into one df
df_plot2 <- bind_rows(plot_data)
df_stats <- bind_rows(stat_data)

# Factor Levels
method_levels <- c("RAW", "CORE", "GREEDY", "CORE_GREEDY", "ALL")
method_labels <- c("Raw", "Core", "Greedy", "Core + Greedy", "All")

df_plot2$method <- factor(df_plot2$method, levels = method_levels, labels = method_labels)
df_stats$method <- factor(df_stats$method, levels = method_levels, labels = method_labels)

# Format the p-values
df_stats <- df_stats %>%
  mutate(label = sprintf("ANOVA %s\nFK %s", 
                         ifelse(anova_p < 0.001, "p < 0.001", sprintf("p = %.3g", anova_p)),
                         ifelse(fligner_p < 0.001, "p < 0.001", sprintf("p = %.3g", fligner_p))))

ylim <- range(df_plot2$resid, na.rm=TRUE); pad <- 0.05 * diff(ylim)

p2 <- ggplot(df_plot2, aes(x = key, y = resid)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  geom_boxplot(fill = "grey90", color = "black", linewidth = 0.35, 
               outlier.color = "black", outlier.size = 0.6, outlier.shape = 16, width = 0.6) +
  facet_wrap(~ method, ncol = 1, scales = "free_x") +
  labs(x = "scanner", y = "residuals (mm)") +
  coord_cartesian(ylim = c(ylim[1] - pad, ylim[2] + pad)) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey90"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_line(color = "black", linewidth = 0.3),
    strip.text = element_text(face = "bold", hjust = 0.5),
    strip.background = element_blank(),
    axis.title = element_text(face = "bold")
  ) +
  geom_text(data = df_stats, aes(x = 1, y = Inf, label = label), 
            hjust = 0, vjust = 1.1, inherit.aes = FALSE, size = 3.5)

ggsave("Figure_2_Boxplots.pdf", p2, width = 10, height = 2 + 2 * length(method_levels))
ggsave("Figure_2_Boxplots.png", p2, width = 10, height = 2 + 2 * length(method_levels), dpi = 300)

message("Figure 2 Generated.")

# ==============================================================================
# FIGURE 3
# ==============================================================================
if (!requireNamespace("ggseg", quietly = TRUE)) stop("Please install ggseg.")
if (!requireNamespace("sf", quietly = TRUE)) stop("Please install sf.")

normalise_region <- function(x) {
  x %>%
    tolower() %>%
    gsub("^thickness\\.(left|right)\\.", "", .) %>%
    gsub("\\.", " ", .) %>%
    trimws()
}

plot_eta_map <- function(eta_named, title, out_stub) {
  df <- tibble(roi = names(eta_named), value = as.numeric(eta_named)) %>%
    mutate(
      hemi = ifelse(grepl("thickness\\.left\\.", roi), "left",
                    ifelse(grepl("thickness\\.right\\.", roi), "right", NA_character_)),
      region = normalise_region(roi)
    ) %>%
    filter(!is.na(hemi)) %>%
    mutate(hemi = tolower(hemi), region = tolower(region))
  
  atlas <- as.data.frame(ggseg::dk)
  atlas$region <- tolower(atlas$region)
  atlas$hemi <- tolower(atlas$hemi)
  atlas2 <- dplyr::left_join(atlas, df, by = c("region", "hemi"))
  atlas_sf <- sf::st_as_sf(atlas2)
  
  p <- ggplot(atlas_sf) +
    geom_sf(aes(fill = value), color = "grey30", linewidth = 0.1) +
    coord_sf() +
    scale_fill_viridis_c(option = "C", limits = c(0, global_max_eta),
                         oob = scales::squish, name = expression(eta^2)) +
    labs(title = title, fill = expression(eta^2)) +
    theme_void(base_size = 11) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  ggsave(paste0(out_stub, ".pdf"), p, width = 8.5, height = 6.5)
  ggsave(paste0(out_stub, ".png"), p, width = 8.5, height = 6.5, dpi = 300)
  p
}

compute_eta2 <- function(X_mat, df_use, batch_col_name) {
  
  # Fit the model with core covaraites and get the residuals
  fit_bio <- lm(as.matrix(X_mat) ~ baselineAge + SEX + DIAGNOSIS, data = df_use)
  resids  <- residuals(fit_bio)
  
  # Get the batches
  batch_factor <- as.factor(df_use[[batch_col_name]])
  
  eta_values <- apply(resids, 2, function(r) {
    # Fit the scanner model
    stats <- anova(lm(r ~ batch_factor))
    
    # Calculate Eta^2
    ss_batch <- stats["batch_factor", "Sum Sq"]
    ss_total <- sum(stats$`Sum Sq`)
    
    return(ss_batch / ss_total)
  })
  
  return(eta_values)
}

eta_results <- list()

# Loop through the datasets already harmonised from before
for (set_name in names(harm_datasets)) {
  
  message(sprintf("Calculating Eta^2 for: %s", set_name))
  
  # Get feature matrix
  X_use <- harm_datasets[[set_name]]
  
  # Remove the features from df_harm
  df_harm_pure <- df_harm %>% select(-any_of(feat_cols))
  
  # Get the Eta^2 values
  eta_vec <- compute_eta2(X_use, df_harm_pure, batch_col)
  
  # Store the results
  eta_results[[set_name]] <- eta_vec
  
  # Report the mean values
  message(sprintf("  -> Mean eta^2 = %.4f", mean(eta_vec)))
}

# Find the global max value to set the scale
global_max_eta <- max(unlist(eta_results), na.rm = TRUE)
message(sprintf("Global Max Eta^2: %.4f", global_max_eta))

# Plot and save individual maps + stacked figure
plot_list <- list()
for (nm in names(eta_results)) {
  title <- switch(nm,
                  RAW = "Unharmonised",
                  CORE = "ComBat (Core)",
                  GREEDY = "ComBat (Greedy)",
                  CORE_GREEDY = "ComBat (Core + Greedy)",
                  ALL = "ComBat (All)",
                  nm)
  plot_list[[nm]] <- plot_eta_map(eta_results[[nm]], title, paste0("Figure_3_eta2_", nm))
}

# Plot and save stacked figure using the cowplot method
if (requireNamespace("cowplot", quietly = TRUE)) {
  message("Stacking plots with cowplot...")
  
  # Extract the legend from the first plot (they all share the same scale)
  shared_legend <- cowplot::get_legend(
    plot_list[[1]] + ggplot2::theme(legend.position = "bottom")
  )
  
  # Strip the legends from all individual plots to prevent clutter
  plots_nolegend <- lapply(plot_list, function(p) {
    p + ggplot2::theme(legend.position = "none")
  })
  
  # Stack the brain maps vertically
  p_stack_main <- cowplot::plot_grid(
    plotlist = plots_nolegend, 
    ncol = 1, 
    align = "v"
  )
  
  # Append the shared legend to the very bottom
  p_stack <- cowplot::plot_grid(
    p_stack_main, 
    shared_legend, 
    ncol = 1, 
    rel_heights = c(1, 0.08) # Gives 8% of the height to the legend
  )
  
  # Calculate dynamic height based on number of plots
  stack_height <- max(6, 2.2 * length(plot_list))
  
  ggsave("Figure_3_eta2_stack.pdf", p_stack, width = 8.5, height = stack_height, bg = "white")
  ggsave("Figure_3_eta2_stack.png", p_stack, width = 8.5, height = stack_height, dpi = 300, bg = "white")
  
}

message("Figure 3 Generated.")

# ==============================================================================
# FIGURE 4
# ==============================================================================

# Targets
targets <- c("baselineAge", "All_Subjects_CDR_07Aug2025__CDRSB", "SEX", "DIAGNOSIS", "Siemens")

# Clean target data
get_target_data <- function(df, tg) {
  y <- df[[tg]]
  
  # Binarise diagnosis
  if (tg == "DIAGNOSIS") {
    y <- factor(ifelse(y == "AD", "AD", "nonAD"), levels = c("nonAD", "AD"))
  }
  
  if (tg == "Siemens") {
    y <- factor(y, levels = c("Other", "Siemens"))
  }
  
  return(y)
}

effect_rows <- list()

# Loop through each dataset
for (set_name in names(harm_datasets)) {
  
  # Get the brain data
  X_use <- harm_datasets[[set_name]]
  
  for (tg in targets) {
    
    # Get the variable of interest
    y <- get_target_data(df_harm, tg)
    
    # Confirm the type
    is_numeric_target <- is.numeric(y)
    
    # Calculate the effect sizes
    stats <- apply(X_use, 2, function(x) {
      if (is_numeric_target) {
        val <- compute_slope(x, as.numeric(y))
        type <- "slope"
      } else {
        # Handle specific reference levels
        if (tg == "DIAGNOSIS") {
          val <- compute_cohen_d(x, y, level0 = "nonAD", level1 = "AD")
        } else if (tg == "Siemens") {
          val <- compute_cohen_d(x, y, level0 = "Other", level1 = "Siemens")
        } else if (tg == "SEX") {
          val <- compute_cohen_d(x, y, level0 = "F", level1 = "M")
        }
        type <- "cohen_d"
      }
      return(c(val, type))
    })
    
    # Store results
    chunk <- tibble(
      dataset = set_name,
      target = tg,
      roi = colnames(stats),
      effect_size = as.numeric(stats[1, ]),
      effect_type = stats[2, ]
    )
    
    effect_rows[[paste(set_name, tg, sep="_")]] <- chunk
  }
}

# Bind rows together
es_df <- bind_rows(effect_rows) %>%
  filter(is.finite(effect_size)) %>%
  mutate(
    target = case_when(
      target == "All_Subjects_CDR_07Aug2025__CDRSB" ~ "CDR-SB",
      target == "baselineAge" ~ "Baseline Age",
      target == "DIAGNOSIS"   ~ "Diagnosis",
      target == "SEX"         ~ "Sex",
      TRUE ~ target
    ),
    dataset = factor(dataset, levels = c("RAW", "CORE", "GREEDY", "CORE_GREEDY", "ALL"),
                     labels = c("Raw", "Core", "Greedy", "Core + Greedy", "All"))
  )

level_order <- c("Raw", "Core", "Greedy", "Core + Greedy", "All")
es_df <- es_df %>% mutate(dataset = factor(dataset, levels = level_order))

p_es <- ggplot(es_df, aes(x = dataset, y = effect_size, fill = dataset)) +
  geom_boxplot(
    width = 0.6,
    outlier.size = 1,
    outlier.alpha = 0.3,
    color = "black",
    linewidth = 0.4
  ) +
  facet_wrap(~ target, scales = "free", ncol = 2) +
  labs(y = "Effect size", x = NULL) +
  scale_fill_discrete(guide = "none") +
  theme_classic(base_size = 14) +
  theme(
    panel.grid       = element_blank(),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 12),
    axis.text.x      = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y      = element_text(color = "black"),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave("Figure_4.pdf", p_es, width = 8, height = 10)
ggsave("Figure_4.png", p_es, width = 8, height = 10, dpi = 300)
message("Figure 4 Generated.")

# ==============================================================================
# FIGURE 5
# ==============================================================================

# Parallelise
parallelise_reps <- TRUE

# 30 random samplings for each p (number of covariates used)
reps_per_p <- 30

# Number of cores to use
max_cores_use <- max(1L, parallel::detectCores(logical = FALSE) - 6L)

# Candidate pool excluding CORE covariates
candidate_all <- setdiff(candidates, CORE_VARS)

# Find the min and max number of covariates includable
p_min <- length(CORE_VARS)
P_max <- length(unique(c(CORE_VARS, candidate_all)))

# Create evenly spaced sequence of 10 values
p_grid <- sort(unique(as.integer(round(seq(p_min, P_max, length.out = 10)))))

# Effect size calculated 
siemens_effect_values <- function(X_mat, y_raw) {
  
  y_site <- factor(y_raw, levels = c("Other", "Siemens"))
  
  d_vals <- apply(X_mat, 2, function(x) {
    compute_cohen_d(x, y_site, level0 = "Other", level1 = "Siemens")
  })
  
  return(d_vals)
}

# Store the raw value values
X_raw_full <- X_raw
y_siemens <- df_harm$Siemens

# Calculate Siemens effect sizes (RAW)
raw_d_vals <- siemens_effect_values(X_raw_full, y_siemens)

# Select the core covariates
core_Z <- df_harm[, CORE_VARS, drop = FALSE]
core_keep <- complete.cases(core_Z)

# Run ComBat with CORE covariates
X_core <- run_combat(X_raw_full[core_keep, , drop = FALSE], df_harm[[batch_col]][core_keep], core_Z[core_keep, , drop = FALSE])

# Calculate Siemens effect sizes (CORE)
core_d_vals <- siemens_effect_values(X_core, y_siemens[core_keep])

# Store results
baseline_df <- data.frame(
  p        = NA_integer_, 
  rep      = NA_integer_, 
  dataset  = c("RAW", "CORE"),
  median_d = c(median(raw_d_vals, na.rm = TRUE), median(core_d_vals, na.rm = TRUE)),
  mean_d   = c(mean(raw_d_vals, na.rm = TRUE), mean(core_d_vals, na.rm = TRUE)),
  d_vals   = I(list(raw_d_vals, core_d_vals))
)

run_one_rep <- function(p, rep_idx) {
  
  # Number of covariates to sample besides CORE
  n_extra <- p - length(CORE_VARS)
  
  # Sample them
  extra_vars <- if (n_extra > 0) sample(candidate_all, size = n_extra, replace = FALSE) else character(0)
  
  # Final p covariates selected
  all_vars_p <- unique(c(CORE_VARS, extra_vars))
  
  # Greedy via marginal screen + greedy selection (shared helper)
  greedy_vars <- select_greedy_vars(
    df_for_cv = df_harm,
    feat_cols = feat_cols,
    candidates = all_vars_p,
    seed_base = 42 + rep_idx,
    keep_threshold = 0.5,
    min_gain = 0.5,
    max_k = 30
  )
  
  # Harmonise ALL
  Z_all <- df_harm[, all_vars_p, drop = FALSE]
  keep_all <- complete.cases(Z_all)
  X_all <- run_combat(X_raw_full[keep_all, , drop = FALSE], 
                             df_harm[[batch_col]][keep_all], 
                             Z_all[keep_all, , drop = FALSE])
  all_d_vals <- siemens_effect_values(X_all, y_siemens[keep_all])
  
  all_row <- data.frame(
    p = p, rep = rep_idx, dataset = "ALL",
    median_d = median(all_d_vals, na.rm = TRUE),
    mean_d   = mean(all_d_vals, na.rm = TRUE),
    d_vals   = I(list(all_d_vals))
  )
  
  # Harmonise GREEDY
  if (length(greedy_vars) > 0) {
    Zg <- df_harm[, greedy_vars, drop = FALSE]
    keep_g <- complete.cases(Zg)
    Xg <- run_combat(X_raw_full[keep_g, , drop = FALSE], 
                            df_harm[[batch_col]][keep_g], 
                            Zg[keep_g, , drop = FALSE])
    g_d_vals <- siemens_effect_values(Xg, y_siemens[keep_g])
  } else {
    g_d_vals <- rep(NA_real_, length(raw_d_vals)) # Return vector of NAs to match shape
  }
  
  greedy_row <- data.frame(
    p = p, rep = rep_idx, dataset = "GREEDY",
    median_d = median(g_d_vals, na.rm = TRUE),
    mean_d   = mean(g_d_vals, na.rm = TRUE),
    d_vals   = I(list(g_d_vals))
  )
  
  # Only return the dynamic results
  rbind(all_row, greedy_row)
}


all_runs <- list()

# Run the grid
for (p in p_grid) {
  
  if (parallelise_reps) {
    # Cluster for this specific p
    n_cores <- min(reps_per_p, max_cores_use)
    cl <- parallel::makeCluster(n_cores)
    
    # Export environment and load libraries
    parallel::clusterExport(cl, varlist = ls(envir = environment()), envir = environment())
    parallel::clusterEvalQ(cl, { 
      library(stats)
      library(dplyr)
      library(tibble)
      library(ComBatFamily) 
    })
    
    # Run parallel apply
    rep_idx <- seq_len(reps_per_p)
    out_rep <- parallel::parLapply(cl, rep_idx, function(r) {
      set.seed(1000L + p * 100L + r)
      run_one_rep(p, r)
    })
    
    # Shut down cluster
    parallel::stopCluster(cl)
    
  } else {
    out_rep <- lapply(seq_len(reps_per_p), function(r) {
      set.seed(1000L + p * 100L + r)
      run_one_rep(p, r)
    })
  }
  
  # Bind the results for this p and finally print the message
  all_runs[[as.character(p)]] <- bind_rows(out_rep)
  message("Done p = ", p)
}

# Merge runs
sim_results <- bind_rows(all_runs)

# Duplicate RAW and CORE for every p
expanded_baselines <- expand.grid(p = p_grid, dataset = c("RAW", "CORE"), stringsAsFactors = FALSE) %>%
  left_join(baseline_df %>% select(-p, -rep), by = "dataset") %>%
  mutate(rep = NA_integer_)

# Combine everything into final_results
final_results <- bind_rows(expanded_baselines, sim_results) %>%
  mutate(
    dataset = factor(dataset, levels = c("RAW", "CORE", "ALL", "GREEDY")),
    p = as.numeric(p)
  )

# Summary
res_sum <- final_results %>%
  group_by(p, dataset) %>%
  summarise(
    mean_val = mean(mean_d, na.rm = TRUE),
    q25_val  = quantile(mean_d, 0.25, na.rm = TRUE),
    q75_val  = quantile(mean_d, 0.75, na.rm = TRUE),
    .groups  = "drop"
  )

p_vary <- ggplot() +
  geom_ribbon(
    data = res_sum,
    aes(x = p, ymin = q25_val, ymax = q75_val, fill = dataset),
    alpha = 0.2,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line(
    data = res_sum,
    aes(x = p, y = mean_val, color = dataset, group = dataset),
    linewidth = 0.7
  ) +
  geom_point(
    data = res_sum,
    aes(x = p, y = mean_val, color = dataset),
    size = 2
  ) +
  scale_color_discrete(name = "Method") +
  scale_fill_discrete(guide = "none") +
  labs(x = "Number of covariates (p)", y = "Mean Siemens effect size (Cohen's d)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("Figure_5.pdf", p_vary, width = 7.5, height = 4.5)
ggsave("Figure_5.png", p_vary, width = 7.5, height = 4.5, dpi = 300)
message("Figure 5 Generated.")

# ==============================================================================
# FIGURE 6
# ==============================================================================

# Simulation parameters
seed <- 42
delta_thickness <- 0.2
prop_rois <- 0.5
sd_eps_W <- 0.4
Delta_scanner <- 1.0

# Get unique subjects and scanners
sub_primary <- df_harm %>%
  select(subid, batch = .data[[batch_col]]) %>%
  distinct(subid, .keep_all = TRUE)

# Assign scanners to -1 and +1 groups
scanners <- sort(unique(sub_primary$batch))
n_scan <- length(scanners)
Gk <- setNames(rep(c(-1, 1), length.out = n_scan), sample(scanners))

# Generate the simulated signal
set.seed(seed)
sub_w <- sub_primary %>%
  mutate(
    # Look up the group for the scanner
    scanner_sign = Gk[batch],
    # Calculate the scanner shift
    scanner_shift = Delta_scanner * scanner_sign,
    # Add noise
    W_SIM = scanner_shift + rnorm(n(), mean = 0, sd = sd_eps_W)
  ) %>%
  select(subid, W_SIM, scanner_group = scanner_sign)

# Drop ROIs, keep covariates and join it with W_SIM
cov_sim <- df_harm %>%
  select(-all_of(feat_cols)) %>%
  left_join(sub_w, by = "subid")

# Determine randomly the affected ROIs
set.seed(seed + 1L)
m_aff <- as.integer(round(prop_rois * length(feat_cols)))
affected_rois <- if (m_aff == 0L) character(0) else sample(feat_cols, size = m_aff, replace = FALSE)

# Inject the signal
X_sim <- X_raw
if (length(affected_rois)) {
  U <- cov_sim$W_SIM
  for (r in affected_rois) {
    idx <- match(r, colnames(X_sim))
    X_sim[, idx] <- as.numeric(X_sim[, idx]) + delta_thickness * U
  }
}

# Create the oracle
set.seed(seed + 2L)
sub_w_oracle <- sub_w

# Shuffle the W_SIM values
sub_w_oracle$W_SIM <- sample(sub_w_oracle$W_SIM, size = nrow(sub_w_oracle), replace = FALSE)

# Replace the old W_SIM values with the new ones for oracle
cov_oracle <- cov_sim %>%
  select(-W_SIM, -scanner_group) %>%
  left_join(sub_w_oracle, by = "subid")

# Like before, inject the signal
X_oracle <- X_raw
if (length(affected_rois)) {
  Uo <- cov_oracle$W_SIM
  for (r in affected_rois) {
    idx <- match(r, colnames(X_oracle))
    X_oracle[, idx] <- as.numeric(X_oracle[, idx]) + delta_thickness * Uo
  }
}

# Get the candidate list for the simulation which involves adding on W_SIM
sim_candidates <- unique(c(candidates, "W_SIM"))

# Select greedy variables
greedy_vars_sim <- select_greedy_vars(
  df_for_cv = cov_sim %>% bind_cols(as.data.frame(X_sim)),
  feat_cols = feat_cols,
  candidates = sim_candidates, 
  seed_base = 42,
  keep_threshold = 0.5,
  min_gain = 0.5,
  max_k = 30
)

# Define the covariate sets
cov_sets_sim <- list(
  RAW    = NULL,
  CORE   = CORE_VARS,
  GREEDY = greedy_vars_sim,
  ALL    = unique(c(CORE_VARS, sim_candidates))
)

method_rows <- list()
method_rows_d <- list()

# Extract the targets (Siemens and W_SIM)
y_sim_num  <- as.numeric(cov_sim$W_SIM)
y_sim_site <- factor(ifelse(cov_sim$Siemens == "Siemens", "Siemens", "Other"), 
                     levels = c("Other", "Siemens"))

# Now determine the effect sizes
for (m in names(cov_sets_sim)) {
  
  # Harmonise
  if (m == "RAW") {
    X_h <- X_sim
  } else {
    cov_names <- cov_sets_sim[[m]]
    current_cov_df <- cov_sim[, cov_names, drop = FALSE]
    
    # Run ComBat
    X_h <- run_combat(X_sim, batch_vec, current_cov_df)
  }
  
  # Calculate slopes
  slopes <- apply(X_h[, affected_rois, drop = FALSE], 2, function(x) {
    compute_slope(x, y_sim_num)
  })
  
  # Calculate Cohen's d
  dvals <- apply(X_h, 2, function(x) {
    compute_cohen_d(x, y_sim_site, level0 = "Other", level1 = "Siemens")
  })
  
  method_rows[[m]]   <- tibble(method = m, roi = affected_rois, slope = as.numeric(slopes))
  method_rows_d[[m]] <- tibble(method = m, roi = feat_cols, d = as.numeric(dvals))
}

# Run the same pipeline for oracle
X_orc_core <- run_combat(X_oracle, batch_vec, cov_oracle[, CORE_VARS, drop = FALSE])
y_orc_num  <- as.numeric(cov_oracle$W_SIM)
y_orc_site <- factor(ifelse(cov_oracle$Siemens == "Siemens", "Siemens", "Other"), 
                     levels = c("Other", "Siemens"))

orc_slopes <- apply(X_orc_core[, affected_rois, drop = FALSE], 2, function(x) {
  compute_slope(x, y_orc_num)
})
orc_dvals <- apply(X_orc_core, 2, function(x) {
  compute_cohen_d(x, y_orc_site, level0 = "Other", level1 = "Siemens")
})

method_rows[["ORACLE"]]   <- tibble(method = "ORACLE", roi = affected_rois, slope = as.numeric(orc_slopes))
method_rows_d[["ORACLE"]] <- tibble(method = "ORACLE", roi = feat_cols, d = as.numeric(orc_dvals))

slope_df <- bind_rows(method_rows)
d_df     <- bind_rows(method_rows_d)

order_methods <- c("RAW", "CORE", "GREEDY", "ORACLE", "ALL")
slope_df$method <- factor(slope_df$method, levels = order_methods)
d_df$method     <- factor(d_df$method, levels = order_methods)

slope_summary <- slope_df %>%
  group_by(method) %>%
  summarise(mean_slope = mean(slope), sd_slope = sd(slope), .groups = "drop")

print(as.data.frame(slope_summary))

p_slope <- ggplot(slope_df, aes(x = method, y = slope, fill = method)) +
  geom_boxplot(outlier.size = 0.7, width = 0.65, color = "black", linewidth = 0.4) +
  labs(x = NULL, y = expression(W[SIM] ~ "slope (signal ROIs)"), title = NULL) +
  scale_fill_discrete(guide = "none") +
  theme_classic(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        strip.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

p_d <- ggplot(d_df, aes(x = method, y = d, fill = method)) +
  geom_boxplot(outlier.size = 0.7, width = 0.65, color = "black", linewidth = 0.4) +
  labs(x = NULL, y = "Siemens Cohen's d (all ROIs)", title = NULL) +
  scale_fill_discrete(guide = "none") +
  theme_classic(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        strip.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

p_slope <- p_slope + annotate("text", x = -Inf, y = Inf, label = "(A)", hjust = -0.2, vjust = 1.2, fontface = "bold")
p_d <- p_d + annotate("text", x = -Inf, y = Inf, label = "(B)", hjust = -0.2, vjust = 1.2, fontface = "bold")
p_pub <- patchwork::wrap_plots(p_slope, p_d, ncol = 2)
print(p_pub)

ggsave("Figure_6_Simulation_Results.pdf", p_pub, width = 11, height = 5.5, device = cairo_pdf)
ggsave("Figure_6_Simulation_Results.png", p_pub, width = 11, height = 5.5, dpi = 300)

message("Figure 6 saved successfully.")



# ==============================================================================
# FIGURE S1
# ==============================================================================

# Simulation parameters
seed <- 42
delta_thickness <- 0.0 # Change: No biological effect
prop_rois <- 0.5
sd_eps_W <- 0.4
Delta_scanner <- 1.0

# Get unique subjects and scanners
sub_primary <- df_harm %>%
  select(subid, batch = .data[[batch_col]]) %>%
  distinct(subid, .keep_all = TRUE)

# Assign scanners to -1 and +1 groups
scanners <- sort(unique(sub_primary$batch))
n_scan <- length(scanners)
Gk <- setNames(rep(c(-1, 1), length.out = n_scan), sample(scanners))

# Generate the simulated signal
set.seed(seed)
sub_w <- sub_primary %>%
  mutate(
    # Look up the group for the scanner
    scanner_sign = Gk[batch],
    # Calculate the scanner shift
    scanner_shift = Delta_scanner * scanner_sign,
    # Add noise
    W_SIM = scanner_shift + rnorm(n(), mean = 0, sd = sd_eps_W)
  ) %>%
  select(subid, W_SIM, scanner_group = scanner_sign)

# Drop ROIs, keep covariates and join it with W_SIM
cov_sim <- df_harm %>%
  select(-all_of(feat_cols)) %>%
  left_join(sub_w, by = "subid")

# Determine randomly the affected ROIs
set.seed(seed + 1L)
m_aff <- as.integer(round(prop_rois * length(feat_cols)))
affected_rois <- if (m_aff == 0L) character(0) else sample(feat_cols, size = m_aff, replace = FALSE)

# Inject the signal
X_sim <- X_raw
if (length(affected_rois)) {
  U <- cov_sim$W_SIM
  for (r in affected_rois) {
    idx <- match(r, colnames(X_sim))
    X_sim[, idx] <- as.numeric(X_sim[, idx]) + delta_thickness * U
  }
}

# Create the oracle
set.seed(seed + 2L)
sub_w_oracle <- sub_w

# Shuffle the W_SIM values
sub_w_oracle$W_SIM <- sample(sub_w_oracle$W_SIM, size = nrow(sub_w_oracle), replace = FALSE)

# Replace the old W_SIM values with the new ones for oracle
cov_oracle <- cov_sim %>%
  select(-W_SIM, -scanner_group) %>%
  left_join(sub_w_oracle, by = "subid")

# Like before, inject the signal
X_oracle <- X_raw
if (length(affected_rois)) {
  Uo <- cov_oracle$W_SIM
  for (r in affected_rois) {
    idx <- match(r, colnames(X_oracle))
    X_oracle[, idx] <- as.numeric(X_oracle[, idx]) + delta_thickness * Uo
  }
}

# Get the candidate list for the simulation which involves adding on W_SIM
sim_candidates <- unique(c(candidates, "W_SIM"))

# Select greedy variables
greedy_vars_sim <- select_greedy_vars(
  df_for_cv = cov_sim %>% bind_cols(as.data.frame(X_sim)),
  feat_cols = feat_cols,
  candidates = sim_candidates, 
  seed_base = 42,
  keep_threshold = 0.5,
  min_gain = 0.5,
  max_k = 30
)

# Define the covariate sets
cov_sets_sim <- list(
  RAW    = NULL,
  CORE   = CORE_VARS,
  GREEDY = greedy_vars_sim,
  ALL    = unique(c(CORE_VARS, sim_candidates))
)

method_rows <- list()
method_rows_d <- list()

# Extract the targets (Siemens and W_SIM)
y_sim_num  <- as.numeric(cov_sim$W_SIM)
y_sim_site <- factor(ifelse(cov_sim$Siemens == "Siemens", "Siemens", "Other"), 
                     levels = c("Other", "Siemens"))

# Now determine the effect sizes
for (m in names(cov_sets_sim)) {
  
  # Harmonise
  if (m == "RAW") {
    X_h <- X_sim
  } else {
    cov_names <- cov_sets_sim[[m]]
    current_cov_df <- cov_sim[, cov_names, drop = FALSE]
    
    # Run ComBat
    X_h <- run_combat(X_sim, batch_vec, current_cov_df)
  }
  
  # Calculate slopes
  slopes <- apply(X_h[, affected_rois, drop = FALSE], 2, function(x) {
    compute_slope(x, y_sim_num)
  })
  
  # Calculate Cohen's d
  dvals <- apply(X_h, 2, function(x) {
    compute_cohen_d(x, y_sim_site, level0 = "Other", level1 = "Siemens")
  })
  
  method_rows[[m]]   <- tibble(method = m, roi = affected_rois, slope = as.numeric(slopes))
  method_rows_d[[m]] <- tibble(method = m, roi = feat_cols, d = as.numeric(dvals))
}

# Run the same pipeline for oracle
X_orc_core <- run_combat(X_oracle, batch_vec, cov_oracle[, CORE_VARS, drop = FALSE])
y_orc_num  <- as.numeric(cov_oracle$W_SIM)
y_orc_site <- factor(ifelse(cov_oracle$Siemens == "Siemens", "Siemens", "Other"), 
                     levels = c("Other", "Siemens"))

orc_slopes <- apply(X_orc_core[, affected_rois, drop = FALSE], 2, function(x) {
  compute_slope(x, y_orc_num)
})
orc_dvals <- apply(X_orc_core, 2, function(x) {
  compute_cohen_d(x, y_orc_site, level0 = "Other", level1 = "Siemens")
})

method_rows[["ORACLE"]]   <- tibble(method = "ORACLE", roi = affected_rois, slope = as.numeric(orc_slopes))
method_rows_d[["ORACLE"]] <- tibble(method = "ORACLE", roi = feat_cols, d = as.numeric(orc_dvals))

slope_df <- bind_rows(method_rows)
d_df     <- bind_rows(method_rows_d)

order_methods <- c("RAW", "CORE", "GREEDY", "ORACLE", "ALL")
slope_df$method <- factor(slope_df$method, levels = order_methods)
d_df$method     <- factor(d_df$method, levels = order_methods)

slope_summary <- slope_df %>%
  group_by(method) %>%
  summarise(mean_slope = mean(slope), sd_slope = sd(slope), .groups = "drop")

print(as.data.frame(slope_summary))

p_slope <- ggplot(slope_df, aes(x = method, y = slope, fill = method)) +
  geom_boxplot(outlier.size = 0.7, width = 0.65, color = "black", linewidth = 0.4) +
  labs(x = NULL, y = expression(W[SIM] ~ "slope (signal ROIs)"), title = NULL) +
  scale_fill_discrete(guide = "none") +
  theme_classic(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        strip.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

p_d <- ggplot(d_df, aes(x = method, y = d, fill = method)) +
  geom_boxplot(outlier.size = 0.7, width = 0.65, color = "black", linewidth = 0.4) +
  labs(x = NULL, y = "Siemens Cohen's d (all ROIs)", title = NULL) +
  scale_fill_discrete(guide = "none") +
  theme_classic(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        strip.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

p_slope <- p_slope + annotate("text", x = -Inf, y = Inf, label = "(A)", hjust = -0.2, vjust = 1.2, fontface = "bold")
p_d <- p_d + annotate("text", x = -Inf, y = Inf, label = "(B)", hjust = -0.2, vjust = 1.2, fontface = "bold")
p_pub <- patchwork::wrap_plots(p_slope, p_d, ncol = 2)
print(p_pub)

ggsave("Figure_S1_Simulation_Results.pdf", p_pub, width = 11, height = 5.5, device = cairo_pdf)
ggsave("Figure_S1_Simulation_Results.png", p_pub, width = 11, height = 5.5, dpi = 300)

message("Figure S1 saved successfully.")

# ==============================================================================
# FIGURE S2
# ==============================================================================

# Simulation parameters
seed <- 42
delta_thickness <- 0.2
prop_rois <- 0.5
sd_eps_W <- 0.4
Delta_scanner <- 1.0

# Get unique subjects and scanners
sub_primary <- df_harm %>%
  select(subid, batch = .data[[batch_col]]) %>%
  distinct(subid, .keep_all = TRUE)

# Assign scanners to -1 and +1 groups
scanners <- sort(unique(sub_primary$batch))
n_scan <- length(scanners)
Gk <- setNames(rep(c(-1, 1), length.out = n_scan), sample(scanners))

# Generate the simulated signal
set.seed(seed)
sub_w <- sub_primary %>%
  mutate(
    # Look up the group for the scanner
    scanner_sign = Gk[batch],
    # Calculate the scanner shift
    scanner_shift = Delta_scanner * scanner_sign,
    # Add noise
    W_SIM = scanner_shift + rnorm(n(), mean = 0, sd = sd_eps_W)
  ) %>%
  select(subid, W_SIM, scanner_group = scanner_sign)

# Drop ROIs, keep covariates and join it with W_SIM
cov_sim <- df_harm %>%
  select(-all_of(feat_cols)) %>%
  left_join(sub_w, by = "subid")

# Determine randomly the affected ROIs
set.seed(seed + 1L)
m_aff <- as.integer(round(prop_rois * length(feat_cols)))
affected_rois <- if (m_aff == 0L) character(0) else sample(feat_cols, size = m_aff, replace = FALSE)

# Inject the signal
X_sim <- harm_datasets[["CORE"]] # Change: use the already harmonised dataset to simulate no site effects
if (length(affected_rois)) {
  U <- cov_sim$W_SIM
  for (r in affected_rois) {
    idx <- match(r, colnames(X_sim))
    X_sim[, idx] <- as.numeric(X_sim[, idx]) + delta_thickness * U
  }
}

# Create the oracle
set.seed(seed + 2L)
sub_w_oracle <- sub_w

# Shuffle the W_SIM values
sub_w_oracle$W_SIM <- sample(sub_w_oracle$W_SIM, size = nrow(sub_w_oracle), replace = FALSE)

# Replace the old W_SIM values with the new ones for oracle
cov_oracle <- cov_sim %>%
  select(-W_SIM, -scanner_group) %>%
  left_join(sub_w_oracle, by = "subid")

# Like before, inject the signal
X_oracle <- harm_datasets[["CORE"]] # Change: use the already harmonised dataset to simulate no site effects
if (length(affected_rois)) {
  Uo <- cov_oracle$W_SIM
  for (r in affected_rois) {
    idx <- match(r, colnames(X_oracle))
    X_oracle[, idx] <- as.numeric(X_oracle[, idx]) + delta_thickness * Uo
  }
}

# Get the candidate list for the simulation which involves adding on W_SIM
sim_candidates <- unique(c(candidates, "W_SIM"))

# Select greedy variables
greedy_vars_sim <- select_greedy_vars(
  df_for_cv = cov_sim %>% bind_cols(as.data.frame(X_sim)),
  feat_cols = feat_cols,
  candidates = sim_candidates, 
  seed_base = 42,
  keep_threshold = 0.5,
  min_gain = 0.5,
  max_k = 30
)

# Define the covariate sets
cov_sets_sim <- list(
  RAW    = NULL,
  CORE   = CORE_VARS,
  GREEDY = greedy_vars_sim,
  ALL    = unique(c(CORE_VARS, sim_candidates))
)

method_rows <- list()
method_rows_d <- list()

# Extract the targets (Siemens and W_SIM)
y_sim_num  <- as.numeric(cov_sim$W_SIM)
y_sim_site <- factor(ifelse(cov_sim$Siemens == "Siemens", "Siemens", "Other"), 
                     levels = c("Other", "Siemens"))

# Now determine the effect sizes
for (m in names(cov_sets_sim)) {
  
  # Harmonise
  if (m == "RAW") {
    X_h <- X_sim
  } else {
    cov_names <- cov_sets_sim[[m]]
    current_cov_df <- cov_sim[, cov_names, drop = FALSE]
    
    # Run ComBat
    X_h <- run_combat(X_sim, batch_vec, current_cov_df)
  }
  
  # Calculate slopes
  slopes <- apply(X_h[, affected_rois, drop = FALSE], 2, function(x) {
    compute_slope(x, y_sim_num)
  })
  
  # Calculate Cohen's d
  dvals <- apply(X_h, 2, function(x) {
    compute_cohen_d(x, y_sim_site, level0 = "Other", level1 = "Siemens")
  })
  
  method_rows[[m]]   <- tibble(method = m, roi = affected_rois, slope = as.numeric(slopes))
  method_rows_d[[m]] <- tibble(method = m, roi = feat_cols, d = as.numeric(dvals))
}

# Run the same pipeline for oracle
X_orc_core <- run_combat(X_oracle, batch_vec, cov_oracle[, CORE_VARS, drop = FALSE])
y_orc_num  <- as.numeric(cov_oracle$W_SIM)
y_orc_site <- factor(ifelse(cov_oracle$Siemens == "Siemens", "Siemens", "Other"), 
                     levels = c("Other", "Siemens"))

orc_slopes <- apply(X_orc_core[, affected_rois, drop = FALSE], 2, function(x) {
  compute_slope(x, y_orc_num)
})
orc_dvals <- apply(X_orc_core, 2, function(x) {
  compute_cohen_d(x, y_orc_site, level0 = "Other", level1 = "Siemens")
})

method_rows[["ORACLE"]]   <- tibble(method = "ORACLE", roi = affected_rois, slope = as.numeric(orc_slopes))
method_rows_d[["ORACLE"]] <- tibble(method = "ORACLE", roi = feat_cols, d = as.numeric(orc_dvals))

slope_df <- bind_rows(method_rows)
d_df     <- bind_rows(method_rows_d)

order_methods <- c("RAW", "CORE", "GREEDY", "ORACLE", "ALL")
slope_df$method <- factor(slope_df$method, levels = order_methods)
d_df$method     <- factor(d_df$method, levels = order_methods)

slope_summary <- slope_df %>%
  group_by(method) %>%
  summarise(mean_slope = mean(slope), sd_slope = sd(slope), .groups = "drop")

print(as.data.frame(slope_summary))

p_slope <- ggplot(slope_df, aes(x = method, y = slope, fill = method)) +
  geom_boxplot(outlier.size = 0.7, width = 0.65, color = "black", linewidth = 0.4) +
  labs(x = NULL, y = expression(W[SIM] ~ "slope (signal ROIs)"), title = NULL) +
  scale_fill_discrete(guide = "none") +
  theme_classic(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        strip.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

p_d <- ggplot(d_df, aes(x = method, y = d, fill = method)) +
  geom_boxplot(outlier.size = 0.7, width = 0.65, color = "black", linewidth = 0.4) +
  labs(x = NULL, y = "Siemens Cohen's d (all ROIs)", title = NULL) +
  scale_fill_discrete(guide = "none") +
  theme_classic(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        strip.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

p_slope <- p_slope + annotate("text", x = -Inf, y = Inf, label = "(A)", hjust = -0.2, vjust = 1.2, fontface = "bold")
p_d <- p_d + annotate("text", x = -Inf, y = Inf, label = "(B)", hjust = -0.2, vjust = 1.2, fontface = "bold")
p_pub <- patchwork::wrap_plots(p_slope, p_d, ncol = 2)
print(p_pub)

ggsave("Figure_S2_Simulation_Results.pdf", p_pub, width = 11, height = 5.5, device = cairo_pdf)
ggsave("Figure_S2_Simulation_Results.png", p_pub, width = 11, height = 5.5, dpi = 300)

message("Figure S2 saved successfully.")




# ==============================================================================
# FIGURE S3
# ==============================================================================

# Simulation parameters
seed <- 42
delta_thickness <- 0.2
prop_rois <- 0.5
sd_eps_W <- 1.08 # Increased to keep marginal variance
Delta_scanner <- 0.0 # Change: No site association

# Get unique subjects and scanners
sub_primary <- df_harm %>%
  select(subid, batch = .data[[batch_col]]) %>%
  distinct(subid, .keep_all = TRUE)

# Assign scanners to -1 and +1 groups
scanners <- sort(unique(sub_primary$batch))
n_scan <- length(scanners)
Gk <- setNames(rep(c(-1, 1), length.out = n_scan), sample(scanners))

# Generate the simulated signal
set.seed(seed)
sub_w <- sub_primary %>%
  mutate(
    # Look up the group for the scanner
    scanner_sign = Gk[batch],
    # Calculate the scanner shift
    scanner_shift = Delta_scanner * scanner_sign,
    # Add noise
    W_SIM = scanner_shift + rnorm(n(), mean = 0, sd = sd_eps_W)
  ) %>%
  select(subid, W_SIM, scanner_group = scanner_sign)

# Drop ROIs, keep covariates and join it with W_SIM
cov_sim <- df_harm %>%
  select(-all_of(feat_cols)) %>%
  left_join(sub_w, by = "subid")

# Determine randomly the affected ROIs
set.seed(seed + 1L)
m_aff <- as.integer(round(prop_rois * length(feat_cols)))
affected_rois <- if (m_aff == 0L) character(0) else sample(feat_cols, size = m_aff, replace = FALSE)

# Inject the signal
X_sim <- X_raw
if (length(affected_rois)) {
  U <- cov_sim$W_SIM
  for (r in affected_rois) {
    idx <- match(r, colnames(X_sim))
    X_sim[, idx] <- as.numeric(X_sim[, idx]) + delta_thickness * U
  }
}

# Create the oracle
set.seed(seed + 2L)
sub_w_oracle <- sub_w

# Shuffle the W_SIM values
sub_w_oracle$W_SIM <- sample(sub_w_oracle$W_SIM, size = nrow(sub_w_oracle), replace = FALSE)

# Replace the old W_SIM values with the new ones for oracle
cov_oracle <- cov_sim %>%
  select(-W_SIM, -scanner_group) %>%
  left_join(sub_w_oracle, by = "subid")

# Like before, inject the signal
X_oracle <- X_raw
if (length(affected_rois)) {
  Uo <- cov_oracle$W_SIM
  for (r in affected_rois) {
    idx <- match(r, colnames(X_oracle))
    X_oracle[, idx] <- as.numeric(X_oracle[, idx]) + delta_thickness * Uo
  }
}

# Get the candidate list for the simulation which involves adding on W_SIM
sim_candidates <- unique(c(candidates, "W_SIM"))

# Select greedy variables
greedy_vars_sim <- select_greedy_vars(
  df_for_cv = cov_sim %>% bind_cols(as.data.frame(X_sim)),
  feat_cols = feat_cols,
  candidates = sim_candidates, 
  seed_base = 42,
  keep_threshold = 0.5,
  min_gain = 0.5,
  max_k = 30
)

# Define the covariate sets
cov_sets_sim <- list(
  RAW    = NULL,
  CORE   = CORE_VARS,
  GREEDY = greedy_vars_sim,
  ALL    = unique(c(CORE_VARS, sim_candidates))
)

method_rows <- list()
method_rows_d <- list()

# Extract the targets (Siemens and W_SIM)
y_sim_num  <- as.numeric(cov_sim$W_SIM)
y_sim_site <- factor(ifelse(cov_sim$Siemens == "Siemens", "Siemens", "Other"), 
                     levels = c("Other", "Siemens"))

# Now determine the effect sizes
for (m in names(cov_sets_sim)) {
  
  # Harmonise
  if (m == "RAW") {
    X_h <- X_sim
  } else {
    cov_names <- cov_sets_sim[[m]]
    current_cov_df <- cov_sim[, cov_names, drop = FALSE]
    
    # Run ComBat
    X_h <- run_combat(X_sim, batch_vec, current_cov_df)
  }
  
  # Calculate slopes
  slopes <- apply(X_h[, affected_rois, drop = FALSE], 2, function(x) {
    compute_slope(x, y_sim_num)
  })
  
  # Calculate Cohen's d
  dvals <- apply(X_h, 2, function(x) {
    compute_cohen_d(x, y_sim_site, level0 = "Other", level1 = "Siemens")
  })
  
  method_rows[[m]]   <- tibble(method = m, roi = affected_rois, slope = as.numeric(slopes))
  method_rows_d[[m]] <- tibble(method = m, roi = feat_cols, d = as.numeric(dvals))
}

# Run the same pipeline for oracle
X_orc_core <- run_combat(X_oracle, batch_vec, cov_oracle[, CORE_VARS, drop = FALSE])
y_orc_num  <- as.numeric(cov_oracle$W_SIM)
y_orc_site <- factor(ifelse(cov_oracle$Siemens == "Siemens", "Siemens", "Other"), 
                     levels = c("Other", "Siemens"))

orc_slopes <- apply(X_orc_core[, affected_rois, drop = FALSE], 2, function(x) {
  compute_slope(x, y_orc_num)
})
orc_dvals <- apply(X_orc_core, 2, function(x) {
  compute_cohen_d(x, y_orc_site, level0 = "Other", level1 = "Siemens")
})

method_rows[["ORACLE"]]   <- tibble(method = "ORACLE", roi = affected_rois, slope = as.numeric(orc_slopes))
method_rows_d[["ORACLE"]] <- tibble(method = "ORACLE", roi = feat_cols, d = as.numeric(orc_dvals))

slope_df <- bind_rows(method_rows)
d_df     <- bind_rows(method_rows_d)

order_methods <- c("RAW", "CORE", "GREEDY", "ORACLE", "ALL")
slope_df$method <- factor(slope_df$method, levels = order_methods)
d_df$method     <- factor(d_df$method, levels = order_methods)

slope_summary <- slope_df %>%
  group_by(method) %>%
  summarise(mean_slope = mean(slope), sd_slope = sd(slope), .groups = "drop")

print(as.data.frame(slope_summary))

p_slope <- ggplot(slope_df, aes(x = method, y = slope, fill = method)) +
  geom_boxplot(outlier.size = 0.7, width = 0.65, color = "black", linewidth = 0.4) +
  labs(x = NULL, y = expression(W[SIM] ~ "slope (signal ROIs)"), title = NULL) +
  scale_fill_discrete(guide = "none") +
  theme_classic(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        strip.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

p_d <- ggplot(d_df, aes(x = method, y = d, fill = method)) +
  geom_boxplot(outlier.size = 0.7, width = 0.65, color = "black", linewidth = 0.4) +
  labs(x = NULL, y = "Siemens Cohen's d (all ROIs)", title = NULL) +
  scale_fill_discrete(guide = "none") +
  theme_classic(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        strip.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

p_slope <- p_slope + annotate("text", x = -Inf, y = Inf, label = "(A)", hjust = -0.2, vjust = 1.2, fontface = "bold")
p_d <- p_d + annotate("text", x = -Inf, y = Inf, label = "(B)", hjust = -0.2, vjust = 1.2, fontface = "bold")
p_pub <- patchwork::wrap_plots(p_slope, p_d, ncol = 2)
print(p_pub)

ggsave("Figure_S3_Simulation_Results.pdf", p_pub, width = 11, height = 5.5, device = cairo_pdf)
ggsave("Figure_S3_Simulation_Results.png", p_pub, width = 11, height = 5.5, dpi = 300)

message("Figure S3 saved successfully.")


