


#This document describes the process and code used to explore PFAS concentrations
#in the Mill River watershed near Amherst, MA using spatial stream network (SSN) models.  
#Specifically, we fit competing model structures using different combinations of fixed 
#effect covariates and determine which are best suited to explain concentration variability 
#of different PFAS family compounds at the sites of each sampling event 
#(S1 = 17, S3 = 33, S4 = 20 sites) across the Mill River network. 

# --------1.  Set Up Analysis --------------- 
library(doParallel)
library(foreach)
library(foreign)
library(ellipse)
library(car)
library(knitr)
library(pander)
library(magrittr)
library(lubridate)
library(tidyverse)
library(htmltools)
library(DT)
library(sf)
library(SSN2)
library(purrr)
library(doSNOW)
library(progress)
panderOptions('digits',8)

#I didn't remove outliers in my SSN models
# since I was thinking that being able to predict why certain locations are extreme
# with our limited covariates would be novel. 
#It may be best to remove the certain outliers, I am not sure
remove_outliers <- FALSE

# ---- Analysis settings ----
# A compound is only modeled for a sampling event if it has at least this many sites with a
# non-zero (above-MDL) result. S1 has 17 sites, S4 has 20, S3 has 33, so 10 drops compounds
# detected at fewer than ~50-60% of the sites in S1/S4.
min_nonzero <- 10

# Estimation method for the final retained models (section 5.4). Model SELECTION always uses
# "ml" (required to compare different fixed effects with AIC). "reml" gives less biased variance
# parameters for the models you report and predict from; use "ml" to reproduce earlier runs.
final_estmethod <- "reml"

opts_chunk$set(tidy.opts=list(width.cutoff=100),tidy=FALSE) # set tidy r code not to run more than 100 characters wide
knitr::opts_chunk$set(message = FALSE) # hide all warning messages in document

# Run from the folder that contains data/ and scripts/
if (basename(getwd()) != "data") setwd("data")
source("../scripts/00 PFAS naming.R")   # PFAS compounds, functional groups and column-naming rules

#------ 2.  Functions used in data analysis---------  

#These functions are used to standardize the continuous covariates for observation 
# (`stand()`) and prediction (`stdpreds()`) sites. The back transformation function converts
# the standardized data model coefficients from standardized form to the raw unit form based on the observation input data that the model was fit to. 
# stand() function to standardize fitting data
source(file = "../scripts/helperfnxs/ssn_standardize_variables_function.R")

# stdpreds() function to standardize prediction data based on the fitting data
source(file = "../scripts/helperfnxs/ssn_standardize_prediction_variables_function.R")

# backtransformation function for generating estimate table
source(file = "../scripts/helperfnxs/ssn_backtransformation_fnx.R")

# formula combination functions
source(file = "../scripts/helperfnxs/ssn_formula_all_combo_function.R")

#-------3.  Import the SSN object --------------- 
#read in the SSN object and calculate the distance matrices for the observation and prediction points

# Import every event's SSN object created in script 03 (../ssnobj/PFAS_NHDHR_<event>.ssn)
ssn_paths <- list.files("../ssnobj", pattern = "^PFAS_NHDHR_.*\\.ssn$", full.names = TRUE)
stopifnot("No .ssn folders found in ../ssnobj - run script 03 first" = length(ssn_paths) > 0)
ssn_list <- ssn_paths |>
  purrr::set_names(sub("^PFAS_NHDHR_(.*)\\.ssn$", "\\1", basename(ssn_paths))) |>
  purrr::map(\(p) SSN2::ssn_import(p, predpts = "preds", overwrite = TRUE))

ssn_list <- ssn_list |>
  purrr::imap(function(ssn_obj, event_name) {
    cat("── Event:", event_name,"──────────────\n")
    
    # Distance matrices (dismat)
    # Required before any ssn_lm() calls.
    # among_predpts = FALSE because we don't need pred-to-pred distances
    # only_predpts  = FALSE because compute obs-to-obs AND obs-to-pred
    
    SSN2::ssn_create_distmat( #building distance matrices
      ssn.object    = ssn_obj,
      predpts       = "preds",
      overwrite     = TRUE,
      among_predpts = FALSE,
      only_predpts  = FALSE)
    
    # Return the ssn object unchanged and distmat is written to the .ssn folder
    ssn_obj
  })


# ------ 4. Bring in covariate data----------  

#Reach-based covariate data to join to observation and predictions sites.  

# Every samplling event gets its own standardized obs + preds SSN object because
#  Site counts differ (S1=17, S3=33, S4=20)
#  Mean and SD of covariates differ per event
#  Standardization needs to be based on the fitting data for that specific event

# Outputs per event:
#   ssn_obj_list[[event]]$ssn_obj_obspreds_std  <- SSN object ready for modeling
#   ssn_obj_list[[event]]$obs_covariates_s      <- standardized obs data frame
#   ssn_obj_list[[event]]$preds_covariates_std  <- standardized preds data frame
#   ssn_obj_list[[event]]$continuous            <- raw continuous obs (for back-transform)
#   ssn_obj_list[[event]]$key_compounds         <-compounds passing variation filter
#   ssn_obj_list[[event]]$low_var_compounds     <- compounds dropped

#---- 4.0 Covariate data-------

# lentic_reach_data.csv is an external input; it is not produced by scripts 01-03.
reach_data <- readRDS("../data/RVsummary_basin_covariates.rds") |>
  left_join(readRDS("../data/RVsummary_rca_covariates.rds"), by = "reach_id") |>
  left_join(read_csv("../data/lentic_reach_data.csv") |>
              dplyr::select(reach_id, is_lentic), by = "reach_id") |>
  left_join(readRDS("../data/RVsummary_riparian_covariates.rds"), by = "reach_id")

#Load covariate name lists made in script 02
cov_lists <- readRDS("../data/covariate_name_lists.RDS")

covariate_groups               <- cov_lists$groups
all_modeling_covariates        <- cov_lists$all_modeling_covariates
continuous_modeling_covariates <- cov_lists$continuous_modeling_covariates
binary_modeling_covariates      <- cov_lists$binary_modeling_covariates

cat("Total continuous covariates:", length(continuous_modeling_covariates), "\n")
cat("Total binary covariates:    ", length(binary_modeling_covariates), "\n")

# PFAS compound definitions (pfas_compounds) come from 00_PFAS_naming.R
log_pfas_compounds <- paste0("log1p_", pfas_compounds)

# 4.1 Per-event covariate conditioning
ssn_obj_list <- ssn_list |>
  purrr::imap(function(ssn_obj, event_name) {
    
    cat("\n── Event:", event_name,"──────────────────\n")
    
    #  4.1.1 Observation covariates
    obs_covariates <- SSN2::ssn_get_data(ssn_obj, name = "obs") |>
      dplyr::rename(site_id = Site_Name) |>
      dplyr::select(
        rid, pid, ratio, snapdist, upDist, afvArea,
        locID, netID, site_id, reach_id, ssnid,
        dplyr::any_of(pfas_compounds),
        netgeom
      ) |>
      dplyr::left_join(reach_data, by = "reach_id")
    
    cat("  Obs covariates joined:",
        nrow(obs_covariates), "sites |",
        ncol(obs_covariates), "columns\n")
    
    #  4.1.2 Outlier removal
    if (isTRUE(remove_outliers) && event_name == "S3") {
      # Outlier removal only meaningful for S3 (full network)
      # ssnid == 6 identified in original script as high-leverage site
      outlier_compounds <- c(
        "PFAS40", "PFCA", "PFBA", "PFBS",
        "PFHpA", "PFHxA", "PFOA", "PFPeA")
      for (cmp in intersect(outlier_compounds, names(obs_covariates))) {
        obs_covariates[[cmp]][obs_covariates$ssnid == 6] <- NA
      }
      cat("  Outlier removal applied (ssnid == 6)\n")
    }
    
    #  4.1.3 Log1p response variables 
    # Applied after outlier removal so NAs propagate correctly
    compounds_present <- intersect(pfas_compounds, names(obs_covariates))
    
    obs_covariates <- obs_covariates |>
      dplyr::mutate(dplyr::across(
        .cols  = dplyr::all_of(compounds_present),
        .fns   = ~log1p(.),
        .names = "log1p_{.col}"
      ))
    
    # check there's no negative raw values
    neg_check <- obs_covariates |>
      sf::st_drop_geometry() |>
      dplyr::select(dplyr::all_of(compounds_present)) |>
      dplyr::summarise(
        dplyr::across(everything(), ~sum(. < 0, na.rm = TRUE)))
    
    neg_cols <- names(neg_check)[neg_check > 0]
    if (length(neg_cols) > 0) {
      warning(event_name, ": negative values in ",
              paste(neg_cols, collapse = ", "),
              " — log1p will produce NaN")
    }
    
    cat("  Response variables created:",
        length(compounds_present), "raw +",
        length(compounds_present), "log1p\n")
    
    #  4.1.4 Low-variation filter
    # Drops compounds before model fitting where ssn_lm will always fail:
    #    All values zero or constant (SD = 0)
    #    Fewer than 3 non-NA observations
    #    Fewer than min_nonzero non-zero (above-MDL) values (set at the top of this script)
    # Checked on RAW scale — log1p version dropped automatically if raw fails
    
    all_response_vars <- c(
      compounds_present,
      paste0("log1p_", compounds_present))
    
    variation_check <- obs_covariates |>
      sf::st_drop_geometry() |>
      dplyr::select(dplyr::all_of(all_response_vars)) |>
      dplyr::summarise(dplyr::across(
        .cols = everything(),
        .fns  = function(x) {
          n_non_na  <- sum(!is.na(x))
          n_nonzero <- sum(x > 0,  na.rm = TRUE)
          sd_val    <- sd(x,       na.rm = TRUE)
          n_non_na < 3 | n_nonzero < min_nonzero | is.na(sd_val) | sd_val == 0
        }
      )) |>
      tidyr::pivot_longer(
        cols      = everything(),
        names_to  = "compound",
        values_to = "low_variation")
    
    low_var_compounds <- variation_check |>
      dplyr::filter(low_variation) |>
      dplyr::pull(compound)
    
    key_compounds <- setdiff(all_response_vars, low_var_compounds)
    
    cat("  Variation filter —",
        length(low_var_compounds), "dropped,",
        length(key_compounds), "retained\n")
    
    if (length(low_var_compounds) > 0) {
      # Only print raw-scale drops. If raw is dropped, log scale is too
      raw_dropped <- low_var_compounds[
        !startsWith(low_var_compounds, "log1p_")
      ]
      if (length(raw_dropped) > 0) {
        cat("  Dropped (raw):",
            paste(raw_dropped, collapse = ", "), "\n")
      }
    }
    
    # 4.1.5 Standardize continuous covariates 
    # Standardization is per-event so mean and SD come from this event's obs data
    # and  coefficients are comparable within an event, not necessarily across events
    
    continuous_vars <- intersect(
      continuous_modeling_covariates,
      names(obs_covariates))
    
    binary_vars <- intersect(
      binary_modeling_covariates,
      names(obs_covariates))
    
    cat("  Continuous vars available:", length(continuous_vars), "/",
        length(continuous_modeling_covariates), "\n")
    
    # Flag any covariates missing from this event — informational only
    missing_covs <- setdiff(continuous_modeling_covariates, names(obs_covariates))
    if (length(missing_covs) > 0) {
      cat(" Not in obs data:", paste(missing_covs, collapse = ", "), "\n")
    }
    
    continuous <- obs_covariates[, continuous_vars]
    
    # Save raw continuous data for back-transformation later
    saveRDS(continuous,
            paste0("../data/", event_name,
                   "_PFAS_obs_continuous_df.RDS"))
    
    # Standardize
    cont_s <- continuous |>
      modify_at(continuous_vars, stand) |>
      rename_with(.fn  = ~paste0(.x, "_s"),
                  .cols = everything())
    # Extract binary variables
    binary_raw <- obs_covariates[, binary_vars, drop = FALSE]
    
    # Join raw + standardized into one object
    obs_covariates_s <- data.frame(obs_covariates, cont_s, check.names = FALSE) |>
      sf::st_as_sf()
    
    saveRDS(obs_covariates_s,
            paste0("../data/", event_name,
                   "_PFAS_obs_continuous_df_standardized.RDS"))
    
    cat(" Continuous vars standardized:", length(continuous_vars), "\n")
    
    #  4.1.6 Put standardized obs back into SSN object
    ssn_obj_obs_std <- SSN2::ssn_put_data(
      obs_covariates_s,
      ssn_obj,
      name        = "obs",
      resize_data = FALSE)
    
    # 4.2 Prediction site covariates
    preds_covariates <- SSN2::ssn_get_data(
      ssn_obj_obs_std,
      name = "preds")
    
    # First, select the explicit columns
    preds_covariates <- preds_covariates |>
      dplyr::select(
        rid, pid, ratio, snapdist, upDist, afvArea,
        locID, netID, reach_id, ssnid, netgeom)   # geometry column kept automatically

    # Then add continuous vars that exist
    if (length(continuous_vars) > 0) {
      preds_covariates <- preds_covariates |>
        dplyr::bind_cols(
          SSN2::ssn_get_data(ssn_obj_obs_std, name = "preds") |>
            sf::st_drop_geometry() |>
            dplyr::select(dplyr::any_of(continuous_vars)))
    }
    # Now join reach_data for binary vars
    preds_covariates <- preds_covariates |>
      dplyr::left_join(
        reach_data |> dplyr::select(reach_id, dplyr::all_of(binary_vars)),
        by = "reach_id")
    
    preds_continuous_vars <- intersect(continuous_vars, names(preds_covariates))
    preds_binary_vars     <- intersect(binary_vars, names(preds_covariates))
    
    # Standardize predictions using this events obs mean and SD
    cont_cov_data_preds <- preds_covariates |>
      sf::st_drop_geometry() |>
      dplyr::select(dplyr::all_of(preds_continuous_vars))
    
    cont_cov_data_preds_stdzd <- stdpreds(
      preds_cont_var_df = sf::st_drop_geometry(cont_cov_data_preds),
      obs_cont_var_df   = sf::st_drop_geometry(
        continuous |> dplyr::select(dplyr::all_of(preds_continuous_vars))
      )
    )
    
    colnames(cont_cov_data_preds_stdzd) <-
      paste0(preds_continuous_vars, "_s")
    
    preds_covariates_std <- dplyr::bind_cols(
      preds_covariates,
      cont_cov_data_preds_stdzd
    ) |>
      sf::st_as_sf()
    
    saveRDS(preds_covariates_std,
            paste0("../data/", event_name,
                   "_PFAS_preds_covariates_df_standardized.RDS"))
    
    # 4.2.1 Put standardized preds back into SSN object
    ssn_obj_obspreds_std <- SSN2::ssn_put_data(
      preds_covariates_std,
      ssn_obj_obs_std,
      name        = "preds",
      resize_data = FALSE
    )
    
    #save everything needed down the line
    list(
      ssn_obj_obspreds_std  = ssn_obj_obspreds_std,
      obs_covariates_s      = obs_covariates_s,
      preds_covariates_std  = preds_covariates_std,
      continuous            = continuous,
      continuous_vars       = continuous_vars,
      binary_vars           = binary_vars,
      preds_continuous_vars = preds_continuous_vars,
      preds_binary_vars     = preds_binary_vars,
      compounds_present     = compounds_present,
      key_compounds         = key_compounds,
      low_var_compounds     = low_var_compounds)
  })

# 4.3 Summary across all sampling events

purrr::imap(ssn_obj_list, function(ev, event_name) {
  cat(event_name, ":\n")
  cat("  Obs sites:            ", nrow(ev$obs_covariates_s), "\n")
  cat("  Pred sites:           ", nrow(ev$preds_covariates_std), "\n")
  cat("  Compounds retained:   ", length(ev$key_compounds), "\n")
  cat("  Compounds dropped:    ", length(ev$low_var_compounds), "\n")
  cat("  Continuous vars:      ", length(ev$continuous_vars), "\n")
})

# ---- 5. Model fitting and selection ----
##--- 5.1 Covariates Forumla list---------
# Covariate groups developed in script 02_SSN_Covariates 
npdes_cov       <- covariate_groups$npdes
impervious_cov  <- covariate_groups$impervious
agriculture_cov <- covariate_groups$agriculture
baseflow_cov    <- covariate_groups$baseflow
soil_cov        <- covariate_groups$soil

# Formula builder appends _s suffix because ssn_lm uses standardized columns.
# Each group contributes one slot in the formula (1 covariate per group max).

# Helper adds "_s" to all standardized covariate names for use in formulas
add_s <- function(x) paste0(x, "_s")

# Build formula layers with one covariate group per layer
# Each layer adds one covariate from that group to existing formulas
ssn_formula_list1 <- paste(" ~", formula_all_combo_fnx(add_s(npdes_cov), 1))

ssn_formula_list2 <- c(
  paste(ssn_formula_list1, "+",
        rep(formula_all_combo_fnx(add_s(impervious_cov), 1),
            each = length(ssn_formula_list1))),
  paste(" ~", formula_all_combo_fnx(add_s(impervious_cov), 1)))

ssn_formula_list3 <- c(
  paste(ssn_formula_list2, "+",
        rep(formula_all_combo_fnx(add_s(agriculture_cov), 1),
            each = length(ssn_formula_list2))),
  paste(" ~", formula_all_combo_fnx(add_s(agriculture_cov), 1)))

ssn_formula_list4 <- c(
  paste(ssn_formula_list3, "+",
        rep(formula_all_combo_fnx(add_s(baseflow_cov), 1),
            each = length(ssn_formula_list3))),
  paste(" ~", formula_all_combo_fnx(add_s(baseflow_cov), 1)))

ssn_formula_list5 <- c(
  paste(ssn_formula_list4, "+",
        rep(formula_all_combo_fnx(add_s(soil_cov), 1),
            each = length(ssn_formula_list4))),
  paste(" ~", formula_all_combo_fnx(add_s(soil_cov), 1)))

# is_lentic is binary, it's NOT standardized so there's no _s suffix
ssn_formula_list6 <- c(
  paste(ssn_formula_list5, "+",
        rep(formula_all_combo_fnx("is_lentic", 1),
            each = length(ssn_formula_list5))),
  paste(" ~", formula_all_combo_fnx("is_lentic", 1)))

# Clean up and restrict forumulas to 2-3 covariates
ssn_formula_list <- gsub(" + 1", "",  ssn_formula_list6, fixed = TRUE) |>
  gsub(" 1 +", "", x = _, fixed = TRUE) |>
  unique() |>
  (\(x) x[order(nchar(x))])() |>
  (\(x) x[stringr::str_count(x, "\\+") %in% c(1, 2)])()

pander(paste("Total formula structures to fit:", length(ssn_formula_list)))
print(head(ssn_formula_list, 10)) # check the first few formulas


##-------5.2.  Fit models  ------------
#We nest the model fit within a diagnostic statistic function to provide the fit and the diagnostics for quick model comparison via AIC.  

# run parallel loop and save results into the following data frame

# Each sampling event runs independently and create their own:
#   SSN object (because they have different obs sites)
#   key_compounds (sample variation differs per event and some don't pass the filtering step)
#   checkpoint directory (so they save seperately)
#   PFAS_ssn_fits output (a fitted model per sampling event)

#I do the cluster creation and model fitting in one go
# It also closes the clusters at the end. 
{if (exists("cl")) try(stopCluster(cl), silent = TRUE)
  
  cl <- makeCluster(4)
  registerDoSNOW(cl)
  parallel::clusterEvalQ(cl, {
    library(SSN2)
    library(dplyr)
    library(purrr)
    library(stringr)
  })
  
  ssn_fits_list <- ssn_obj_list |>
    purrr::imap(function(ev, event_name) {
      
      cat("  Fitting models — Event:", event_name, "\n")
      
      ssn_obj_obspreds_std <- ev$ssn_obj_obspreds_std
      key_compounds        <- ev$key_compounds
      
      n_compounds <- length(key_compounds)
      n_formulas  <- length(ssn_formula_list)
      total_fits  <- n_compounds * n_formulas
      
      cat("  Compounds:", n_compounds,
          "| Formulas:", n_formulas,
          "| Total fits:", total_fits, "\n")
      
      # Event-specific checkpoint directory. One file per compound NAME (not index),
      # so a changed compound list can never load another compound's results.
      results_dir <- file.path("../outputs", paste0("ssn_fits_progress_", event_name))
      dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
      
      # If the data, formula set or compound list changed since
      # the checkpoints were written, clear them so everything is refit.
      # (Only checkpoints from a crashed run with IDENTICAL inputs are reused.)
      fit_signature <- rlang::hash(list(
        ssn_formula_list, key_compounds,
        sf::st_drop_geometry(SSN2::ssn_get_data(ssn_obj_obspreds_std, name = "obs"))))
      sig_file <- file.path(results_dir, "_signature.RDS")
      if (!file.exists(sig_file) || !identical(readRDS(sig_file), fit_signature)) {
        unlink(list.files(results_dir, pattern = "\\.rds$", full.names = TRUE))
        saveRDS(fit_signature, sig_file)
        cat("  New run or inputs changed - checkpoints cleared\n")
      }
      
      existing <- file.exists(file.path(results_dir,
                                        paste0("compound_", key_compounds, ".rds")))
      cat("  Checkpoints already done:", sum(existing), "/", n_compounds, "\n")
      
      # The standardized SSN object lives only in memory (the *_s columns are not on
      # disk), so it is exported to the workers explicitly through foreach's .export.
      
      # Progress bar
      pb <- progress::progress_bar$new(
        format = paste0("  ", event_name,
                        " [:bar] :current/:total fits | :percent | ETA: :eta"),
        total  = total_fits, clear  = FALSE, width  = 75)
      pb$tick(0)
      
      prog <- function(n) pb$tick(n_formulas)
      opts <- list(progress = prog)
      
      #Parallel model fitting loop
      t0 <- Sys.time() #this is for tracking how long it takes to run the model
      
      PFAS_ssn_fits <-
        foreach(i = 1:n_compounds,
                .combine = dplyr::bind_rows,
                .packages = c("SSN2", "dplyr", "purrr", "stringr"),
                .export = c("ssn_obj_obspreds_std", "ssn_formula_list",
                            "key_compounds", "results_dir"),
                .options.snow = opts) %dopar% {
                  out_file <- file.path(results_dir,
                                        paste0("compound_", key_compounds[[i]], ".rds"))
                  
                  # Resume after crash and skips if already done
                  if (file.exists(out_file)) return(readRDS(out_file))
                  
                  formula_list <- lapply(
                    paste(key_compounds[[i]], ssn_formula_list),
                    as.formula
                  )
                  
                  result <- purrr::map(.x = formula_list,
                                       .f = function(x) {
                                         tryCatch({SSN2::glance(SSN2::ssn_lm(
                                           formula     = x,
                                           ssn.object  = ssn_obj_obspreds_std,
                                           tailup_type = "exponential",
                                           additive    = "afvArea",
                                           estmethod   = "ml")) |>
                                             dplyr::mutate(
                                               response_var  = stringr::str_split_i(
                                                 deparse1(x), " ~ ", 1),
                                               predictor_var = stringr::str_split_i(
                                                 deparse1(x), " ~ ", 2))
                                         },
                                         error = function(e) {
                                           data.frame(
                                             response_var  = stringr::str_split_i(
                                               deparse1(x), " ~ ", 1),
                                             predictor_var = stringr::str_split_i(
                                               deparse1(x), " ~ ", 2),
                                             error_msg     = conditionMessage(e))
                                         })
                                       }) |> purrr::list_rbind()
                  
                  saveRDS(result, out_file)
                  result}
      
      t1 <- Sys.time()
      e0 <- t1 - t0
      cat("\n ", event_name, "—", total_fits,
          "fits in", round(e0, 2), attr(e0, "units"), "\n")
      
      # Rebuild from checkpoints if the models already been run
      #this is helpful if it crashes halfway through a run
      PFAS_ssn_fits <- file.path(results_dir,
                                 paste0("compound_", key_compounds, ".rds")) |>
        purrr::keep(file.exists) |>
        purrr::map(readRDS) |>
        dplyr::bind_rows()
      
      saveRDS(PFAS_ssn_fits,
              paste0("../outputs/", event_name, "_PFAS_ssn_fits_all.RDS"))
      
      cat("  Rows in fits table:", nrow(PFAS_ssn_fits), "\n")
      PFAS_ssn_fits
    })
  
  stopCluster(cl)}

##-- 5.3. Select best models based on AIC for each compound-------- 

#We use the delta AIC value of 2 to retain models that have comparable fit to the observed data. At the same time we also calculate the AIC weight of each model among the retained models for a given response variable or PFAS compound/family.  

#According to Burnham and Anderson (2002), the AIC weights for a set of competing models can be calculated as:  
#  $$w_i = \frac{exp(-\frac{1}{2}\Delta_i)}{\sum_{r=1}^R exp(-\frac{1}{2}\Delta_r)}$$  
#   where $w_i$ is the AIC weight for the *i^th^* model relative to the model with the lowest AIC value, *R* is the total number of models in the set being considered for the AIC weighting, $exp(-\frac{1}{2}\Delta_i)$ is the Likelihood of model *i* and $\sum_{r=1}^R exp(-\frac{1}{2}\Delta_r)$ is the sum of all likelihoods in the set of models being considered.  

#select best models based on AIC
best_models_list <- ssn_fits_list |>
  purrr::imap(function(PFAS_ssn_fits, event_name) {
    
    best <- PFAS_ssn_fits |>
      dplyr::group_by(response_var) |>
      dplyr::group_split() |>
      purrr::map(function(x) {
        x |>
          dplyr::arrange(AIC) |>
          dplyr::filter(AIC < min(AIC, na.rm = TRUE) + 2) |>
          tibble::rownames_to_column("rank") |>
          dplyr::mutate(
            delta_aic    = AIC - min(AIC, na.rm = TRUE),
            likelihood_i = exp(-0.5 * delta_aic),
            AIC_weight_i = likelihood_i / sum(likelihood_i, na.rm = TRUE)
          ) |>
          dplyr::relocate(response_var)
      })
    
    names(best) <- PFAS_ssn_fits |>
      dplyr::group_by(response_var) |>
      dplyr::group_keys() |>
      dplyr::pull(response_var)
    
    # Save per-event output
    dir.create("../outputs", showWarnings = FALSE)
    saveRDS(
      best |> dplyr::bind_rows(),
      paste0("../outputs/", event_name,
             "_PFAS_best_model_fit_diagnostics_table.RDS"))
    cat(event_name, "— compounds with retained models:", length(best), "\n")
    best
  })

DT::datatable(
  best_models_list[["S3"]] |>
    dplyr::bind_rows() |>
    dplyr::mutate(dplyr::across(where(is.numeric), ~round(.x, 2))),
  extensions = "FixedColumns",
  options    = list(scrollX = TRUE,
                    fixedColumns = list(leftColumns = 1)),
  class      = "display nowrap compact")

#Each row in these tables represents a model retained based on its delta AIC value (<= 2 AIC units from minimum AIC).  

##----- 5.4.  Best fit models fit individually---------------   

#Fit only the models retained by AIC comparison from above for model evaluation and diagnostics.  
#Selection (5.2) used ML; these final fits use final_estmethod set at the top of this script (default REML).
#The AIC weights and pseudo-R2 used downstream always come from the ML selection table; AIC in the
#6.1 glance tables comes from these final fits and is only comparable within one estimation method.

#r fit suite of best models, cache = TRUE}
# create a named list of formula from best model suites

fit_model_ls_list <- best_models_list |>
  purrr::imap(function(best_PFAS_models, event_name) {
    
    # Pull a specific event's SSN object
    ssn_obj_obspreds_std <- ssn_obj_list[[event_name]]$ssn_obj_obspreds_std
    
    # Build named formula lists per compound
    best_model_formula_ls <- purrr::map(best_PFAS_models, function(x) {
      
      if (nrow(x) == 0) return(list())
      
      form_list <- x |>
        dplyr::mutate(formula = paste(response_var, "~", predictor_var)) |>
        dplyr::select(formula) |>
        as.list() |> unlist() |> as.list() |>
        lapply(FUN = as.formula)
      
      names(form_list) <- x |>
        dplyr::mutate(model_rank = paste0(response_var, "_", rank)) |>
        dplyr::select(model_rank) |>
        as.list() |> unlist()
      form_list})
    
    # Fit each retained model
    fit_model_ls <- purrr::map2(
      .x = best_model_formula_ls,
      .y = names(best_model_formula_ls),
      .f = function(x, y) {
        
        if (length(x) == 0) {
          message(event_name, " — Skipping ", y,
                  " — no valid candidate models.")
          return(NULL)}
        
        loop_model_ls        <- vector("list", length(x))
        names(loop_model_ls) <- names(x)
        
        for (i in seq_len(length(x))) {
          loop_model_ls[[i]] <- SSN2::ssn_lm(
            formula     = x[[i]],
            ssn.object  = ssn_obj_obspreds_std,  # ← per-event object
            tailup_type = "exponential",
            additive    = "afvArea",
            estmethod   = final_estmethod)  # selection above always used "ml"
        }
        loop_model_ls
      })
    
    # Record and drop the skipped compounds
    skipped <- names(fit_model_ls)[sapply(fit_model_ls, is.null)]
    fit_model_ls <- fit_model_ls[!sapply(fit_model_ls, is.null)]
    
    cat(event_name, "— models fit:", length(fit_model_ls),
        "| skipped:", length(skipped), "\n")
    
    if (length(skipped) > 0) {
      cat("  Skipped:", paste(skipped, collapse = ", "), "\n")
    }
    
    list(
      fit_model_ls      = fit_model_ls,
      skipped_compounds = skipped)
  })


#----- 6. Model evaluation and diagnostics---------

##--- 6.1. Individual model summaries ---------
# Calculate the model diagnostics for each individual model in each suite of
# models for a PFAS compound/family. Runs on ALL compounds (both raw and
# log1p scales) before scale comparison in Section 6.2. 

dir.create("../outputs/model_summ_tables", showWarnings = FALSE, recursive = TRUE)
indv_model_summary_ls_list <- fit_model_ls_list |>
  purrr::imap(function(ev, event_name) {
    
    fit_model_ls          <- ev$fit_model_ls
    best_PFAS_models      <- best_models_list[[event_name]]
    continuous            <- ssn_obj_list[[event_name]]$continuous
    
    best_PFAS_models_filtered <- best_PFAS_models[
      !names(best_PFAS_models) %in% ev$skipped_compounds]
    
    purrr::map2(
      .x = fit_model_ls,
      .y = best_PFAS_models_filtered,
      .f = function(x, y) {
        
        modsuite_glance_df   <- data.frame()
        modsuite_loocv_df    <- data.frame()
        modsuite_varcomp_df  <- data.frame()
        modsuite_esttable_df <- data.frame()
        
        for (i in seq_len(length(x))) {
          tryCatch({  # tryCatch around each individual model
            
            loop_glance_df <- SSN2::glance(x[[i]]) |>
              dplyr::mutate(modelrank    = i,
                            response_var = y$response_var[i],
                            modelformula = y$predictor_var[i])
            modsuite_glance_df <- dplyr::bind_rows(modsuite_glance_df,
                                                   loop_glance_df)
            
            loop_loocv_df <- SSN2::loocv(x[[i]]) |>
              dplyr::mutate(modelrank    = i,
                            response_var = y$response_var[i],
                            modelformula = y$predictor_var[i])
            modsuite_loocv_df <- dplyr::bind_rows(modsuite_loocv_df,
                                                  loop_loocv_df)
            
            loop_varcomp_df <- SSN2::varcomp(x[[i]]) |>
              dplyr::mutate(modelrank    = i,
                            response_var = y$response_var[i],
                            modelformula = y$predictor_var[i])
            modsuite_varcomp_df <- dplyr::bind_rows(modsuite_varcomp_df,
                                                    loop_varcomp_df)
            
            ssn_tidy_out      <- SSN2::tidy(x[[i]], conf.int = TRUE)
            continuous_nogeom <- sf::st_drop_geometry(continuous)
            esttable <- std_to_raw_estimate_table_fnx(
              continuous_nogeom, ssn_tidy_out, 5) |>
              dplyr::mutate(modelrank    = i,
                            response_var = y$response_var[i],
                            modelformula = y$predictor_var[i])
            modsuite_esttable_df <- dplyr::bind_rows(modsuite_esttable_df,
                                                     esttable)
            
          }, error = function(e) {
            message(event_name, " — skipping model ", i,
                    " (", y$response_var[i], "): ", conditionMessage(e))
          })
        }
        
        # Only save CSVs if we got at least one successful model
        if (nrow(modsuite_glance_df) > 0) {
          out_prefix <- paste0("../outputs/model_summ_tables/",
                               event_name, "_",
                               y$response_var[nrow(modsuite_glance_df)])
          readr::write_csv(modsuite_glance_df,
                           paste0(out_prefix, "_modsuite_glance.csv"))
          readr::write_csv(modsuite_loocv_df,
                           paste0(out_prefix, "_modsuite_loocv.csv"))
          readr::write_csv(modsuite_varcomp_df,
                           paste0(out_prefix, "_modsuite_varcomp.csv"))
          readr::write_csv(modsuite_esttable_df,
                           paste0(out_prefix, "_modsuite_esttable.csv"))
        }
        
        list(modsuite_glance_df   = modsuite_glance_df,
             modsuite_loocv_df    = modsuite_loocv_df,
             modsuite_varcomp_df  = modsuite_varcomp_df,
             modsuite_esttable_df = modsuite_esttable_df)
      })
  })

# The following note comes from the function after this step.
#    When the response variable is not transformed, 
#    scaled, ortransformed, the '(Intercept)' coefficient values
#    are thesame for both standardized and raw columns in the table.

##----6.2. Comparing Raw vs log1p scale models--------

# AIC cannot compare log vs linear models - different likelihoods, different
# scales. We use AIC-weighted LOOCV RMSPE from Section 6.1 (true out-of-sample
# performance) as the comparison metric.
#
# Raw RMSPE is in ng/L; log1p RMSPE is in log(ng/L+1). Direct comparison is
# not possible, so we normalize each RMSPE by the mean observed value on its own scale
# to produce a dimensionless CV-RMSPE, comparable across both scales and all
# 49 compounds.
#
# Compounds with less than 5% CV-RMSPE difference default to the raw scale for
# interpretability in ng/L
#
# final_key_compounds gives one winning variable per compound
#(one model per compound, either log or raw)

final_key_compounds_list <- best_models_list |>
  purrr::imap(function(best_PFAS_models, event_name) {
    
    indv_model_summary_ls <- indv_model_summary_ls_list[[event_name]]
    obs_covariates_s      <- ssn_obj_list[[event_name]]$obs_covariates_s
    fit_model_ls          <- fit_model_ls_list[[event_name]]$fit_model_ls
    
    # Get raw compounds present in this event
    event_pfas_compounds <- ssn_obj_list[[event_name]]$compounds_present
    
    scale_comparison_ls <- list()
    
    for (compound in event_pfas_compounds) {
      
      log_compound <- paste0("log1p_", compound)
      
      has_raw <- compound     %in% names(indv_model_summary_ls)
      has_log <- log_compound %in% names(indv_model_summary_ls)
      
      if (!has_raw || !has_log) {
        warning(event_name, ": missing summaries for ",
                compound, " — skipping.")
        next
      }
      
      
      raw_loocv <- indv_model_summary_ls[[compound]]$modsuite_loocv_df
      log_loocv <- indv_model_summary_ls[[log_compound]]$modsuite_loocv_df
      
      if (nrow(raw_loocv) == 0 || !"modelrank" %in% names(raw_loocv) ||
          nrow(log_loocv) == 0 || !"modelrank" %in% names(log_loocv)) {
        warning(event_name, ": empty LOOCV data for ",
                compound, " — skipping scale comparison.")
        next
      }
      
      loocv_raw <- indv_model_summary_ls[[compound]]$modsuite_loocv_df |>
        dplyr::left_join(
          best_PFAS_models[[compound]] |>
            dplyr::select(rank, AIC_weight_i) |>
            dplyr::mutate(rank = as.integer(rank)),
          by = c("modelrank" = "rank")) |>
        dplyr::summarise(
          rmspe_wtd = sum(RMSPE * AIC_weight_i, na.rm = TRUE),
          bias_wtd  = sum(bias  * AIC_weight_i, na.rm = TRUE))
      
      loocv_log <- indv_model_summary_ls[[log_compound]]$modsuite_loocv_df |>
        dplyr::left_join(
          best_PFAS_models[[log_compound]] |>
            dplyr::select(rank, AIC_weight_i) |>
            dplyr::mutate(rank = as.integer(rank)),
          by = c("modelrank" = "rank")) |>
        dplyr::summarise(
          rmspe_wtd = sum(RMSPE * AIC_weight_i, na.rm = TRUE),
          bias_wtd  = sum(bias  * AIC_weight_i, na.rm = TRUE))
      
      log_models  <- fit_model_ls[[log_compound]]
      log_weights <- best_PFAS_models[[log_compound]] |>
        dplyr::mutate(rank = as.integer(rank)) |>
        dplyr::arrange(rank) |>
        dplyr::pull(AIC_weight_i)
      
      rmspe_log_ng_per_L <- purrr::map2_dbl(
        log_models, log_weights,
        function(m, w) {
          aug     <- SSN2::augment(m) |> sf::st_drop_geometry()
          obs_ng  <- expm1(aug[[log_compound]])
          pred_ng <- expm1(aug$.fitted)
          sqrt(mean((obs_ng - pred_ng)^2, na.rm = TRUE)) * w
        }
      ) |> sum()
      
      obs_mean_raw <- obs_covariates_s |> sf::st_drop_geometry() |>
        dplyr::pull(dplyr::all_of(compound)) |> mean(na.rm = TRUE)
      obs_mean_log <- obs_covariates_s |> sf::st_drop_geometry() |>
        dplyr::pull(dplyr::all_of(log_compound)) |> mean(na.rm = TRUE)
      
      cv_rmspe_raw <- loocv_raw$rmspe_wtd / obs_mean_raw
      cv_rmspe_log <- loocv_log$rmspe_wtd / obs_mean_log
      pct_diff     <- round(100 * (cv_rmspe_raw - cv_rmspe_log) / cv_rmspe_raw, 1)
      better       <- ifelse(cv_rmspe_log < cv_rmspe_raw, "log1p", "raw")
      
      scale_comparison_ls[[compound]] <- data.frame(
        compound             = compound,
        rmspe_raw            = round(loocv_raw$rmspe_wtd,  3),
        rmspe_log1p          = round(loocv_log$rmspe_wtd,  3),
        rmspe_log1p_ng_per_L = round(rmspe_log_ng_per_L,   3),
        cv_rmspe_raw         = round(cv_rmspe_raw,          4),
        cv_rmspe_log1p       = round(cv_rmspe_log,          4),
        bias_raw             = round(loocv_raw$bias_wtd,    3),
        bias_log1p           = round(loocv_log$bias_wtd,    3),
        better_scale         = better,
        cv_rmspe_pct_diff    = pct_diff)
    }
    
    scale_comparison_df <- dplyr::bind_rows(scale_comparison_ls) |>
      dplyr::arrange(desc(abs(cv_rmspe_pct_diff)))
    
    saveRDS(scale_comparison_df,
            paste0("../outputs/", event_name,
                   "_PFAS_scale_comparison_raw_vs_log1p.RDS"))
    
    cat("\n── Event:", event_name,
        "— scale comparison ─────────────────────────\n")
    cat("  log1p better:  ",
        sum(scale_comparison_df$better_scale == "log1p"), "\n")
    cat("  raw better:    ",
        sum(scale_comparison_df$better_scale == "raw"), "\n")
    cat("  within 5% diff:",
        sum(abs(scale_comparison_df$cv_rmspe_pct_diff) < 5),
        "(→ raw for interpretability)\n")
    
    best_scale_lookup <- scale_comparison_df |>
      dplyr::mutate(
        better_scale = ifelse(abs(cv_rmspe_pct_diff) < 5, "raw", better_scale)) |>
      dplyr::transmute(
        compound,
        best_response_var = ifelse(better_scale == "log1p",
                                   paste0("log1p_", compound),
                                   compound),
        better_scale,
        cv_rmspe_pct_diff)
    
    list(
      scale_comparison_df = scale_comparison_df,
      best_scale_lookup   = best_scale_lookup,
      final_key_compounds = best_scale_lookup$best_response_var)
  })

##---- 6.3. Relative importance of covariates-------
# Burnham and Anderson (2002): relative importance is the sum of AIC weights
# across all retained models that include that variable.

rel_imp_ls_list <- final_key_compounds_list |>
  purrr::imap(function(ev, event_name) {
    
    final_key_compounds <- ev$final_key_compounds
    best_PFAS_models    <- best_models_list[[event_name]]  # ← defined here
    
    rel_imp_df_ls <- best_PFAS_models[final_key_compounds] |>
      purrr::map(function(x) {

        presence_mat <- x |> # Build presence matrix
          dplyr::group_by(rank) |>
          dplyr::group_split() |>
          purrr::map(function(row) {
            data.frame(
              covariate  = strsplit(row$predictor_var,
                                    split = " + ", fixed = TRUE)[[1]],
              presence   = 1,
              model_name = paste0(row$response_var, "_", row$rank))
          }) |>
          dplyr::bind_rows() |>
          tidyr::pivot_wider(names_from  = model_name,
                             values_from = presence,
                             values_fill = 0)
        
        #2: Transpose — rows = models, cols = covariates
        mat_t <- presence_mat |>
          tibble::column_to_rownames("covariate") |>
          t() |>
          data.frame() |>
          tibble::rownames_to_column("model_name")
        
        #3: Join AIC weights
        weights_df <- x |>
          dplyr::mutate(model_name = paste0(response_var, "_", rank)) |>
          dplyr::select(model_name, AIC_weight_i, response_var,
                        rank, predictor_var)
        
        mat_weighted <- mat_t |>
          dplyr::left_join(weights_df, by = "model_name")
        
        # Step 4: Multiply each covariate column by AIC_weight_i
        meta_cols      <- c("model_name", "AIC_weight_i", "response_var",
                            "rank", "predictor_var")
        covariate_cols <- setdiff(names(mat_weighted), meta_cols)
        
        mat_weighted[covariate_cols] <- mat_weighted[covariate_cols] *
          mat_weighted$AIC_weight_i
        
        # Step 5: Sum → relative importance
        mat_weighted |>
          dplyr::select(dplyr::all_of(covariate_cols)) |>
          dplyr::summarise(dplyr::across(everything(),
                                         ~sum(.x, na.rm = TRUE))) |>
          tidyr::pivot_longer(cols      = dplyr::everything(),
                              names_to  = "covariate",
                              values_to = "rel_imp") |>
          dplyr::arrange(desc(rel_imp)) |>
          dplyr::mutate(response_var = x$response_var[1]) |>
          dplyr::relocate(response_var, covariate, rel_imp)
      })
    
    # Relative importance plot
    rel_imp_df_plotting <- purrr::list_rbind(rel_imp_df_ls) |>
      tidyr::pivot_wider(id_cols     = "covariate",
                         values_from = rel_imp,
                         names_from  = response_var,
                         values_fill = 0) |>
      dplyr::mutate(
        overall_mean = rowMeans(dplyr::across(where(is.numeric)),
                                na.rm = TRUE),
        covariate    = forcats::fct_reorder(covariate, overall_mean,
                                            .desc = FALSE)) |>
      tidyr::pivot_longer(cols      = !covariate,
                          names_to  = "response_var",
                          values_to = "rel_imp")
    
    relimp_plot <- ggplot(
      data = rel_imp_df_plotting |>
        dplyr::filter(response_var != "overall_mean"),
      aes(x = covariate, y = rel_imp, color = response_var)) +
      geom_point(size = 2.5, alpha = 0.7) +
      coord_flip() +
      theme_bw() +
      labs(y     = "Relative importance",
           x     = "Covariate",
           color = "Response variable",
           title = paste0("Event ", event_name,
                          " — AIC-weighted covariate importance")) +
      theme(legend.position = "right",
            axis.text.y     = element_text(size = 9))
    
    dir.create("../figs", showWarnings = FALSE)
    ggsave(relimp_plot,
           filename = paste0("../figs/", event_name,
                             "_Relative_Importance_plot.png"),
           width = 8, height = 5, units = "in", dpi = 300)
    
    cat(event_name, "— rel importance plot saved\n")
    
    list(rel_imp_df_ls       = rel_imp_df_ls,
         rel_imp_df_plotting  = rel_imp_df_plotting,
         relimp_plot          = relimp_plot)
  })


##----- 6.4. Averaged model suite summaries-----------
# AIC-weighted mean diagnostics for each compound's best model suite.
# Both map2() inputs filtered to final_key_compounds from Section 6.2.

dir.create("../outputs/model_summ_tables", showWarnings = FALSE, recursive = TRUE)
mean_model_summary_ls_list <- final_key_compounds_list |>
  purrr::imap(function(ev, event_name) {
    
    final_key_compounds  <- ev$final_key_compounds
    fit_model_ls         <- fit_model_ls_list[[event_name]]$fit_model_ls
    indv_model_summary_ls <- indv_model_summary_ls_list[[event_name]]
    continuous           <- ssn_obj_list[[event_name]]$continuous
    
    purrr::map2(
      .x = fit_model_ls[final_key_compounds],
      .y = indv_model_summary_ls[final_key_compounds],
      .f = function(x, y) {
        
        avg_esttable_df <- purrr::map2(x, names(x), function(m, nm) {
          std_to_raw_estimate_table_fnx(
            continuous_vars_df = sf::st_drop_geometry(continuous),
            tidy_out_tibble    = SSN2::tidy(m, conf.int = TRUE),
            roundval           = 5) |>
            dplyr::mutate(model_name = nm)
        }) |>
          dplyr::bind_rows() |>
          # reframe() evaluates arguments in order, so the *_sd columns must be
          # computed before the mean overwrites std_est / raw_est
          dplyr::reframe(.by    = term,
                         model_n    = dplyr::n(),
                         n_sig_0.05 = sum(p_val < 0.05, na.rm = TRUE),
                         std_est_sd = sd(std_est,      na.rm = TRUE),
                         raw_est_sd = sd(raw_est,      na.rm = TRUE),
                         std_est    = mean(std_est,    na.rm = TRUE),
                         raw_est    = mean(raw_est,    na.rm = TRUE),
                         t_stat     = mean(t_stat,     na.rm = TRUE)) |>
          dplyr::relocate(std_est_sd, .after = std_est) |>
          dplyr::relocate(raw_est_sd, .after = raw_est) |>
          dplyr::arrange(desc(model_n)) |>
          dplyr::rename(covariate = term)
        
        avg_loocv_df <- dplyr::bind_rows(
          dplyr::summarise(y$modsuite_loocv_df,
                           dplyr::across(
                             -dplyr::starts_with(c("response_var",
                                                   "modelformula")),
                             ~mean(.x, na.rm = TRUE))) |>
            dplyr::mutate(statistic = "mean"),
          dplyr::summarise(y$modsuite_loocv_df,
                           dplyr::across(
                             -dplyr::starts_with(c("response_var", "model")),
                             ~sd(.x, na.rm = TRUE))) |>
            dplyr::mutate(statistic = "sd")) |>
          dplyr::select(statistic, dplyr::everything(), -modelrank) |>
          tidyr::pivot_longer(cols      = -statistic,
                              names_to  = "LOOCV Statistic",
                              values_to = "values") |>
          tidyr::pivot_wider(names_from  = statistic,
                             values_from = values)
        
        avg_varcomp_df <- y$modsuite_varcomp_df |>
          dplyr::reframe(.by       = varcomp,
                         prop_mean = mean(proportion, na.rm = TRUE),
                         prop_sd   = sd(proportion,   na.rm = TRUE),
                         prop_min  = min(proportion,   na.rm = TRUE),
                         prop_max  = max(proportion,   na.rm = TRUE)) |>
          dplyr::rename(variance_component = varcomp)
        
        out_prefix <- paste0("../outputs/model_summ_tables/", event_name, "_",
                             y$modsuite_glance_df$response_var[1])
        readr::write_csv(avg_esttable_df, paste0(out_prefix, "_modsuite_esttable_MEANS.csv"))
        readr::write_csv(avg_loocv_df, paste0(out_prefix, "_modsuite_loocv_MEANS.csv"))
        readr::write_csv(avg_varcomp_df, paste0(out_prefix, "_modsuite_varcomp_MEANS.csv"))
        
        list(avg_esttable_df  = avg_esttable_df,
             avg_loocv_df     = avg_loocv_df,
             avg_varcomp_df   = avg_varcomp_df)
      })
  })


##------ 6.5. Model diagnostics plots------------------
# Obs/pred, residuals vs fitted, Q-Q, and Cook's Distance plots.

diagnostics_plot_ls_list <- final_key_compounds_list |>
  purrr::imap(function(ev, event_name) {
    
    final_key_compounds <- ev$final_key_compounds
    fit_model_ls        <- fit_model_ls_list[[event_name]]$fit_model_ls
    best_PFAS_models    <- best_models_list[[event_name]]
    obs_covariates      <- ssn_obj_list[[event_name]]$obs_covariates_s
    
    purrr::map2(
      .x = fit_model_ls[final_key_compounds],
      .y = best_PFAS_models[final_key_compounds],
      .f = function(x, y) {
        
        all_model_obspreds_df <- purrr::map2(x, names(x), function(m, nm) {
          SSN2::augment(m) |>
            dplyr::select(pid, observations = 1, dplyr::starts_with(".")) |>
            dplyr::mutate(model_name = nm) |>
            sf::st_drop_geometry() |>
            data.frame()
        }) |>
          dplyr::bind_rows() |>
          dplyr::rename(predictions = .fitted) |>
          dplyr::mutate(pred_type = "Individual model predictions",
                        pid       = as.numeric(pid))
        
        avg_model_obspreds_df <- all_model_obspreds_df |>
          dplyr::left_join(
            y |> dplyr::mutate(
              model_name = paste0(response_var, "_", rank)),
            by = "model_name") |>
          dplyr::mutate(
            weighted_pred     = predictions * AIC_weight_i,
            weighted_resid    = .resid      * AIC_weight_i,
            weighted_hat      = .hat        * AIC_weight_i,
            weighted_cooksd   = .cooksd     * AIC_weight_i,
            weighted_stdresid = .std.resid  * AIC_weight_i) |>
          dplyr::reframe(
            .by          = pid,
            observations = mean(observations, na.rm = TRUE),
            predictions  = sum(weighted_pred),
            .resid       = sum(weighted_resid),
            .hat         = sum(weighted_hat),
            .cooksd      = sum(weighted_cooksd),
            .std.resid   = sum(weighted_stdresid),
            model_name   = "Model_Average") |>
          dplyr::mutate(pred_type = "AIC-weighted predictions",
                        pid       = as.numeric(pid))
        
        obspreds_df <- dplyr::bind_rows(
          all_model_obspreds_df, avg_model_obspreds_df) |>
          dplyr::left_join(
            dplyr::select(obs_covariates, pid, site_id),
            by = "pid")
        
        threshold <- 4 / nrow(avg_model_obspreds_df)
        
        obspred_plot <- ggplot(obspreds_df,
                               aes(x = predictions, y = observations)) +
          geom_point() + theme_minimal() +
          geom_abline(intercept = 0, slope = 1) +
          facet_wrap(~pred_type, ncol = 2) +
          labs(title = paste0(y$response_var[1], ": Observed vs Predicted"),
               x = "Predicted", y = "Observed")
        
        residfit_plot <- ggplot(obspreds_df,
                                aes(x = predictions, y = .resid)) +
          geom_point(alpha = 0.5) +
          geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
          geom_smooth(method = "loess", color = "blue", se = FALSE) +
          theme_minimal() + facet_wrap(~pred_type, ncol = 2) +
          labs(title = paste0(y$response_var[1], ": Residuals vs Fitted"),
               x = "Fitted", y = "Residuals")
        
        qq_plot <- ggplot(obspreds_df, aes(sample = .resid)) +
          stat_qq() + stat_qq_line(color = "red") +
          theme_minimal() + facet_wrap(~pred_type, ncol = 2) +
          labs(title = paste0(y$response_var[1], ": Normal Q-Q"))
        
        cooksd_plot <- ggplot(obspreds_df,
                              aes(x = site_id, y = .cooksd)) +
          geom_hline(yintercept = threshold, linetype = "dashed",
                     color = "red", linewidth = 0.8) +
          geom_segment(aes(xend = site_id, y = 0, yend = .cooksd),
                       color = "gray60", alpha = 0.7) +
          geom_point(aes(color = .cooksd > threshold),
                     size = 2.5, alpha = 0.5, show.legend = FALSE) +
          scale_color_manual(values = c("FALSE" = "black",
                                        "TRUE"  = "firebrick")) +
          theme_minimal(base_size = 13) + coord_flip() +
          facet_wrap(~pred_type, ncol = 2) +
          labs(title    = paste0(y$response_var[1], ": Cook's Distance"),
               subtitle = paste0("Threshold = 4/n = ",
                                 round(threshold, 4)), x = "Site", y = "Weighted Cook's D")
        
        plot_prefix <- paste0("../figs/diagnostics/",
                              event_name, "_", y$response_var[1])
        dir.create("../figs/diagnostics", showWarnings = FALSE, recursive = TRUE)
        
        ggsave(paste0(plot_prefix, "_obspred.png"), obspred_plot,  width = 8, height = 4, dpi = 200)
        ggsave(paste0(plot_prefix, "_residfit.png"), residfit_plot, width = 8, height = 4, dpi = 200)
        ggsave(paste0(plot_prefix, "_qq.png"), qq_plot, width = 8, height = 4, dpi = 200)
        ggsave(paste0(plot_prefix, "_cooksd.png"), cooksd_plot, width = 8, height = 8, dpi = 200)
        
        list(obspreds_df   = obspreds_df,
             obspred_plot  = obspred_plot,
             residfit_plot = residfit_plot,
             qq_plot       = qq_plot,
             cooksd_plot   = cooksd_plot) })
  })


##------ 7. Save objects needed by script 05 ------------------
# Script 05 loads these automatically if they are not already in the session, so it can
# be run on its own (e.g. after restarting R). Fitted models are the slow part - keep them.
dir.create("../outputs", showWarnings = FALSE)
for (nm in c("ssn_obj_list", "fit_model_ls_list", "best_models_list",
             "final_key_compounds_list")) {
  saveRDS(get(nm), file.path("../outputs", paste0(nm, ".RDS")))
}
# Script 05 only needs the obs/pred tables from the diagnostics, not the ggplot objects
# (which carry large environments), so save just those in the same list structure.
saveRDS(
  purrr::map(diagnostics_plot_ls_list,
             \(ev) purrr::map(ev, \(x) list(obspreds_df = x$obspreds_df))),
  "../outputs/diagnostics_plot_ls_list.RDS")

