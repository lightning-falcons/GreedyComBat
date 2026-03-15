suppressPackageStartupMessages({
  library(tidyverse)
  library(vroom)
})

# ==============================================================================
# Globals
# ==============================================================================

# Minimum number of patients required to keep a medication or medical history level
MIN_RECC_ENTRIES <- 2

# File paths
tables_dir  <- "Tables_07Aug2025"
registry_fp <- "REGISTRY_11Mar2026.csv"
ants_fp     <- "ADNI_ANTsSST_protocol_2019_03_14.csv" 
schema_fp   <- "adni_covariates_curated_schema_v5.csv"
out_dir     <- "combat_outputs"
if (!dir.exists(out_dir)) dir.create(out_dir)

# ID cleaner
norm_id  <- function(x) str_replace_all(str_trim(toupper(x)), "[^A-Z0-9]+", "_") |> str_replace_all("^_|_$", "")

# VISCODE cleaner
norm_vis <- function(x) str_remove_all(toupper(x), "\\s+")

# RID cleaner to store as an integer (character type)
norm_num <- function(x) suppressWarnings(as.character(as.integer(x)))

# Fast read everything as char
read_fast <- function(fp) vroom(fp, col_types = cols(.default = "c"), show_col_types = FALSE, progress = FALSE)

# ==============================================================================
# ANTs and Registry Setup
# ==============================================================================
# Ingest and filter the ANTs (cortical thickness) file for those with a valid ID and scan date
ants_raw <- read_fast(ants_fp) %>%
  mutate(PTID_norm = norm_id(ID), ANTS_SCAN_DATE = as.Date(dateAcquired)) %>%
  filter(!is.na(PTID_norm), !is.na(ANTS_SCAN_DATE))

# Identify only the earliest scan for each person
ants_first <- ants_raw %>%
  group_by(PTID_norm) %>%
  summarise(first_ants_scan = min(ANTS_SCAN_DATE, na.rm = TRUE), .groups = "drop") %>%
  distinct(PTID_norm, .keep_all = TRUE)

# Get the list of subids of interest
ptids_of_interest <- ants_first$PTID_norm

# Extract the registry for our subids of interest, which is used to link the various files
reg <- read_fast(registry_fp) %>%
  mutate(
    PTID_norm = norm_id(PTID),
    RID_norm  = norm_num(RID),
    VIS2_norm = norm_vis(VISCODE2),
    VIS_norm  = norm_vis(VISCODE),
    REG_DATE  = as.Date(EXAMDATE)
  ) %>% 
  filter(PTID_norm %in% ptids_of_interest)

# Create and RID <-> PTID dictionary
crosswalk <- reg %>% filter(!is.na(PTID_norm), !is.na(RID_norm)) %>% distinct(RID_norm, PTID_norm)

# Create a lookup table of visit dates for each patient visit of interest
reg_vis2 <- reg %>% filter(!is.na(VIS2_norm), !is.na(REG_DATE)) %>% group_by(PTID_norm, VIS2_norm) %>% summarise(REG_DATE_v2 = min(REG_DATE, na.rm=T), .groups="drop")
reg_vis  <- reg %>% filter(!is.na(VIS_norm), !is.na(REG_DATE)) %>% group_by(PTID_norm, VIS_norm) %>% summarise(REG_DATE_v1 = min(REG_DATE, na.rm=T), .groups="drop")

# ==============================================================================
# Table Linkage
# ==============================================================================
link_table <- function(fp) {
  
  # Read the file
  df <- read_fast(fp)
  
  # Strip the path
  src <- tools::file_path_sans_ext(basename(fp))
  
  # Normalises the columns for one file with different schema
  if (src == "All_Subjects_Study_Entry_07Aug2025") {
    if ("subject_id" %in% names(df)) df <- rename(df, PTID = subject_id)
    if ("entry_date" %in% names(df)) df <- rename(df, EXAMDATE = entry_date)
  }
  
  # Safe column extraction (returns NA if column is missing)
  pull_col <- function(col) if (col %in% names(df)) df[[col]] else NA_character_
  
  # Pull the consistent data headers for easier linking later
  df_std <- tibble(
    PTID_norm = norm_id(pull_col("PTID")),
    RID_norm  = norm_num(pull_col("RID")),
    VIS2_norm = norm_vis(pull_col("VISCODE2")),
    VIS_norm  = norm_vis(pull_col("VISCODE")),
    TAB_DATE  = coalesce(as.Date(pull_col("EXAMDATE")), as.Date(pull_col("VISDATE"))),
    update_stamp = pull_col("UPDATE_STAMP"),
    USERDATE  = pull_col("USERDATE")
  ) %>%
    # Attach the raw data
    bind_cols(df) %>% 
    
    # Patch RID -> PTID if the PTID is missing
    rows_patch(crosswalk, by = "RID_norm", unmatched = "ignore") %>%
    
    # Filter to only those scans in ANTs
    filter(!is.na(PTID_norm), PTID_norm %in% ptids_of_interest) %>%
    
    # Filter to only baseline records
    filter(VIS2_norm %in% c("SC", "BL") | VIS_norm %in% c("SC", "BL") | (is.na(VIS2_norm) & is.na(VIS_norm))) %>%
    
    # Record the source file
    mutate(source_file = src)
  
  # Do not join if no data remains
  if (nrow(df_std) == 0) return(NULL)
  
  df_std %>%
    
    # Join with the ANTs file
    inner_join(ants_first, by = "PTID_norm") %>%
    
    # Join with the VISCODEs
    left_join(reg_vis2, by = c("PTID_norm", "VIS2_norm")) %>%
    left_join(reg_vis, by = c("PTID_norm", "VIS_norm")) %>%
    mutate(
      # Find the date of the visit
      REG_DATE = coalesce(REG_DATE_v2, REG_DATE_v1),
      
      # Find
      diff_days = as.numeric(REG_DATE - first_ants_scan)
    ) %>%
    
    # Delete the intermediary columns
    select(-REG_DATE_v2, -REG_DATE_v1)
}

# ==============================================================================
# Filtering Covariate Records and medhist/meds one-hot encoding
# ==============================================================================
# Get every .csv file (all our data files)
all_csv <- list.files(tables_dir, pattern = "\\.csv$", full.names = TRUE)
files <- all_csv[basename(all_csv) != basename(registry_fp)]

# Read each file according to the link_table function and calculate date differences
linked_raw <- map_dfr(files, link_table) %>%
  group_by(source_file, PTID_norm) %>%
  mutate(
    abs_diff = abs(diff_days),
    tab_reg_diff = abs(as.numeric(TAB_DATE - REG_DATE))
  )

# For standard files (one row per visit). Keep only the single best row.
df_single <- linked_raw %>%
  filter(!str_detect(source_file, "(?i)All_Subjects_RECMHIST|All_Subjects_RECCMEDS")) %>%
  arrange(
    is.na(abs_diff), abs_diff,
    is.na(tab_reg_diff), tab_reg_diff,
    desc(as.Date(update_stamp, tryFormats = c("%Y-%m-%d", "%m/%d/%Y", "%Y/%m/%d"))),
    desc(as.Date(USERDATE, tryFormats = c("%Y-%m-%d", "%m/%d/%Y", "%Y/%m/%d"))),
    .by_group = TRUE
  ) %>%
  slice(1) %>%
  ungroup()

# For multi-row files (Meds/History). Keep ALL rows on the best date.
df_multi <- linked_raw %>%
  filter(str_detect(source_file, "(?i)All_Subjects_RECMHIST|All_Subjects_RECCMEDS")) %>%
  filter(if (all(is.na(abs_diff))) TRUE else is.na(abs_diff) | abs_diff == min(abs_diff, na.rm = TRUE)) %>%
  filter(if (all(is.na(tab_reg_diff))) TRUE else is.na(tab_reg_diff) | tab_reg_diff == min(tab_reg_diff, na.rm = TRUE)) %>%
  ungroup()

# Recombine and replace sentinel values with NA
linked_cross_sectional <- bind_rows(df_single, df_multi) %>%
  mutate(across(everything(), ~ ifelse(str_trim(.x) %in% c("-1", "-4", "-7", "-8", "-9", ""), NA_character_, .x)))

# Name cleaner
clean_feature <- function(x) {
  x <- str_trim(toupper(as.character(x)))
  x[is.na(x) | x %in% c("-1", "-4", "-7", "-8", "-9", "")] <- "MISSING"
  str_replace_all(x, "[^A-Z0-9]+", "_")
}

recc_binary <- linked_cross_sectional %>%
  
  # Extract the medical history and medications files
  filter(str_detect(source_file, "(?i)All_Subjects_RECMHIST|All_Subjects_RECCMEDS")) %>%
  mutate(
    
    # Flag that indicates if the file is medical history
    is_hx = str_detect(source_file, "(?i)RECMHIST"),
    
    # Look at MHNUM or CMMED depending on the file and clean the name
    feature = clean_feature(ifelse(is_hx, pull(., MHNUM), pull(., CMMED))),
    
    # Build the final name
    original_var = paste(ifelse(is_hx, "HX", "MED"), feature, sep = "_"),
    
    # Set the value to 1 since the patients has this history / med
    value = "1"
  ) %>%
  
  # Remove rows without data
  filter(feature != "MISSING") %>%  
  
  # Remove repeats
  distinct(PTID_norm, source_file, original_var, value) %>%
  
  # Retain only rows corresponding to those which appear at least MIN_RECC_ENTRIES times
  group_by(original_var) %>% filter(n() >= MIN_RECC_ENTRIES) %>% ungroup()

# Remove empty strings (sentinels were nuked above)
resolve_valid <- function(x) {
  v <- x[!is.na(x) & str_trim(x) != ""]
  if(length(v) == 0) NA_character_ else v[1]
}

raw_matrix <- linked_cross_sectional %>%
  
  # Filter to all the non-medhist/meds files
  filter(!str_detect(source_file, "(?i)All_Subjects_RECCMEDS|All_Subjects_RECMHIST")) %>%
  
  # Remove all columns which look like dates/timestamps/not a covariate
  select(-any_of(c("RID_norm", "VIS2_norm", "VIS_norm", "TAB_DATE", "REG_DATE", "first_ants_scan", "diff_days", "abs_diff", "tab_reg_diff", "update_stamp", "USERDATE", "USERDATE2", "ID", "RECNO", "SITEID", "PHASE", "subject_id", "ACCNO", "PTDOB"))) %>%
  select(-matches("(?i)date|year|time|stamp")) %>%
  
  # Turn everything into a long format
  mutate(across(everything(), as.character)) %>%
  pivot_longer(-c(PTID_norm, source_file), names_to = "original_var", values_drop_na = TRUE) %>%
  filter(str_trim(value) != "") %>%
  
  # Add on the meds/medhist rows from earlier
  bind_rows(recc_binary) %>%
  
  # Create unique column names
  mutate(final_col_name = paste0(source_file, "__", original_var)) %>%
  
  # Make the table wide again
  pivot_wider(id_cols = PTID_norm, names_from = final_col_name, values_from = value, values_fn = resolve_valid) %>%
  
  # Join with ANTs file
  right_join(tibble(PTID_norm = ptids_of_interest), by = "PTID_norm") %>%
  
  # Set medications or med history that isn't applicable as 0
  mutate(across(matches("__(MED|HX)_"), ~ replace_na(.x, "0")))

# ==============================================================================
# Unit Conversion and Missingness Filter
# ==============================================================================
# Column targets
v_wt <- "All_Subjects_VITALS_07Aug2025__VSWEIGHT"
v_wu <- "All_Subjects_VITALS_07Aug2025__VSWTUNIT"
v_ht <- "All_Subjects_VITALS_07Aug2025__VSHEIGHT"
v_hu <- "All_Subjects_VITALS_07Aug2025__VSHTUNIT"
v_tp <- "All_Subjects_VITALS_07Aug2025__VSTEMP"
v_tu <- "All_Subjects_VITALS_07Aug2025__VSTMPUNT"

# Extract cols returning NA if not available
pull_safe <- function(col) if (col %in% names(raw_matrix)) raw_matrix[[col]] else NA

raw_matrix <- raw_matrix %>%
  mutate(
    
    # Weight conversion (to kg)
    raw_wt = as.numeric(pull_safe(v_wt)),
    wu     = pull_safe(v_wu),
    weight_kg = ifelse(wu == "1" | (is.na(wu) & raw_wt > 120), raw_wt * 0.453592, raw_wt),
    
    # Height conversion (to cm)
    raw_ht = as.numeric(pull_safe(v_ht)),
    hu     = pull_safe(v_hu),
    height_cm = ifelse(hu == "1" | (is.na(hu) & raw_ht < 90), raw_ht * 2.54, raw_ht),
    
    # Temp conversion (to C)
    raw_tp = as.numeric(pull_safe(v_tp)),
    tu     = pull_safe(v_tu),
    temp_c    = ifelse(tu == "1" | (is.na(tu) & raw_tp > 90), (raw_tp - 32) * 5/9, raw_tp)
  ) %>%
  
  # Drop raw converted vitals
  select(-any_of(c(v_wt, v_wu, v_ht, v_hu, v_tp, v_tu, "raw_wt", "wu", "raw_ht", "hu", "raw_tp", "tu")))

# Remove columns which have too many missing values
cols_to_keep <- names(raw_matrix)[sapply(names(raw_matrix), function(c) {
  
  # Always keep  our IDs, cleaned vitals, and med hist / meds
  if (c %in% c("PTID_norm", "weight_kg", "height_cm", "temp_c") || str_detect(c, "(?i)__MED_|__HX_")) return(TRUE)
  
  # For everything else, drop if too many missing
  mean(is.na(raw_matrix[[c]])) < 0.05
})]

# Restrict the matrix to only as above
final_covs <- raw_matrix %>% select(all_of(cols_to_keep))

# ==============================================================================
# Merge
# ==============================================================================
# Extract the cols of interest from ANTs
ants <- read_fast(ants_fp) %>%
  mutate(
    IMAGE_ID = as.character(IMAGE_ID),
    PTID_norm = norm_id(ID),
    subid = as.character(subid),
    baselineAge = as.numeric(baselineAge),
    DIAGNOSIS = factor(DIAGNOSIS),
    SEX = factor(SEX, levels = c("F", "M")),
    Siemens = factor(ifelse(str_detect(tolower(manufac.model.coil.strength.site), "siemens"), "Siemens", "Other"), levels = c("Other","Siemens"))
  )

# Clean any invisible whitespace from column names in both dataframes
names(final_covs) <- trimws(names(final_covs))

# Join the two tables to create the big table
combat_ready <- ants %>%
  select(IMAGE_ID, ID, subid, baselineAge, SEX, DIAGNOSIS, Siemens, PTID_norm, manufac.model.coil.strength.site) %>%
  left_join(final_covs, by = "PTID_norm") %>%
  select(-PTID_norm) %>%
  distinct(IMAGE_ID, .keep_all = TRUE) %>%
  
  # Remove technical columns, comments, dates, and IDs
  select(-matches("(?i)\\.row_in_file|^n_ants|__PTID|__RID|__VISCODE|__Phase|^raw_|^wu|^hu|^tu|__COMMENT|__KitID|__TestID|__BATCH|__SOURCE|__APRECEIVE|__APAMBTEMP|__PTNOTRT|__PTSOURCE|__BSXONSET|__APTESTDT|__CMT|__HMT|__UAT|__SLT|__ORT"))

# ==============================================================================
# Data Imputation
# ==============================================================================

# Read the schema to get the data types
if (!file.exists(schema_fp)) stop("ERROR: Schema file not found!")
schema <- read_csv(schema_fp, show_col_types = FALSE)

type_map <- schema %>% 
  mutate(variable=str_trim(variable), vic_type=str_trim(tolower(vic_type))) %>% 
  select(variable, vic_type) %>% 
  deframe()

meta_cols <- c("IMAGE_ID", "ID", "subid", "Siemens", "manufac.model.coil.strength.site", "baselineAge", "SEX", "DIAGNOSIS")

schema_allowed_cols <- c(meta_cols, names(type_map))

# These are columns that we would keep but are now dropping due to schema omission
dropped_by_schema <- setdiff(names(combat_ready), schema_allowed_cols)

if (length(dropped_by_schema) > 0) {
  message("The following ", length(dropped_by_schema), " variable(s) to be dropped")
  message("are NOT defined in your manual schema ('", schema_fp, "').")
  print(sort(dropped_by_schema)) 
}

# Remove columns that are not allowed
combat_ready <- combat_ready %>% select(any_of(schema_allowed_cols))

# Imputation based on type
final_imputed <- combat_ready %>%
  mutate(across(everything(), ~ {
    cname <- cur_column()
    
    # Do NOT impute core columns
    if (cname %in% meta_cols) return(.x)
    
    # Get the type of the column
    ctype <- type_map[[cname]]
    if (is.null(ctype) || is.na(ctype)) return(.x)
    
    # Mean imputation for numerics
    if (str_detect(ctype, "numeric|integer|double")) {
      val <- suppressWarnings(as.numeric(.x))
      m <- mean(val, na.rm = TRUE)
      val[is.na(val)] <- if(!is.nan(m)) m else 0
      return(val)
    } 
    
    # Factor imputation by the new Unknown category
    if (str_detect(ctype, "factor|character")) {
      val <- as.character(.x)
      val[is.na(val) | str_trim(val) == ""] <- "Unknown"
      return(val)
    }
    
    return(.x)
  }))

# Save the processed file
imputed_out_fp <- file.path(out_dir, "combat_covariates_v7_FINAL_IMPUTED.csv")
write_csv(final_imputed, imputed_out_fp)

message("Success! Pipeline Complete.")
message("- Applied curated schema: ", schema_fp)
message("- Final Imputed Matrix Saved to: ", imputed_out_fp)
