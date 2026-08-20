
#This document describes the process and code used to explore PFAS concentrations in the Mill River watershed near Amherst, MA using spatial stream network (SSN) models.  
#Specifically, we fit competing model structures using different combinations of fixed effect covariates and determine which are best suited to explain concentration variability of different PFAS family compounds at 33 sites across the Mill River network. The results presented below offer some evidence for what landscape features may be related to PFAS concentrations.  
  
  #  1.  Set Up Analysis  
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
panderOptions('digits',8)

remove_outliers <- FALSE

# set tidy r code not to run more than 100 characters wide
opts_chunk$set(tidy.opts=list(width.cutoff=100),tidy=FALSE) 

# hide all warning messages in document
knitr::opts_chunk$set(message = FALSE)
setwd("C:/Users/Marston User/Documents/LWMR Isoscape/PFAS_MillRiver_SSN-main/PFAS_MillRiver_SSN-main/data")

#  2.  Functions used in data analysis  

#These functions are used to standardize the continuous covariates for observation (`stand()`) and prediction (`stdpreds()`) sites. The back transformation function converts the standardized data model coefficients from standardized form to the raw unit form based on the observation input data that the model was fit to.  

# stand() function to standardize fitting data
source(file = "../scripts/helperfnxs/ssn_standardize_variables_function.R")

# stdpreds() function to standardize prediction data based on the fitting data
source(file = "../scripts/helperfnxs/ssn_standardize_prediction_variables_function.R")

# backtransformation function for generating estimate table
source(file = "../scripts/helperfnxs/ssn_backtransformation_fnx.R")

# formula combination functions
source(file = "../scripts/helperfnxs/ssn_formula_all_combo_function.R")

#  3.  Import the SSN object  

#read in the SSN object and calculate the distance matrices for the observation and prediction points. Next plot the network and observation sites.  

#r load_ssn_object_and_plot_data, warning=FALSE, message=FALSE, fig.cap="Mill River network with 33 observation sites (dots)."}
# Load the ssn and all sets of prediction points
ssn_obj <- ssn_import("../ssnobj/PFAS_NHDHR.ssn", include_obs = TRUE,
                      predpts = c("preds"), overwrite = TRUE )

# create distance matrix for ssn object and "preds" prediction points
ssn_create_distmat(ssn_obj, predpts = "preds", overwrite = TRUE, 
                   among_predpts = FALSE, only_predpts = FALSE) # ~ 1 min

# some summary plot for initial data exploration
ggplot() +
  geom_sf(data = ssn_obj$edges) +
  geom_sf(data = ssn_obj$obs, color = "brown", size = 2) +
  theme_bw()

#  4.  Read in and condition covariate data  

#Reach-based covariate data to join to observation and predictions sites.  

reach_data <- read_csv("../data/RVsummary_basin_covariates.csv") |>
  left_join(read_csv("../data/RVsummary_rca_covariates.csv"),
            by = "reach_id") |>
  left_join(select(read_csv("../data/lentic_reach_data.csv"), reach_id, is_lentic),
            by = "reach_id")

# Check that new covariates are present
names(reach_data)

pfas_compounds <- c(
  "PFAS40","X11Cl.PF3OUdS_Results","X3.3_FTCA_Results","X4.2FTS_Results",
  "X5.3_FTCA_Results","X6.2FTS_Results","X7.3_FTCA_Results","X8.2FTS_Results",
  "X9Cl.PF3ONS_Results", "ADONA_Results","HFPO.DA_Results","NEtFOSA_Results",
  "NEtFOSAA_Results","NEtFOSE_Results","NFDHA_Results","NMeFOSA_Results",
  "NMeFOSAA_Results","NMeFOSE_Results","PFBA_Results","PFBS_Results",
  "PFDA_Results", "PFDoA_Results","PFDoS_Results","PFDS_Results", "PFEESA_Results",
  "PFHpA_Results","PFHpS_Results","PFHxA_Results","PFHxS_Results","PFMBA_Results",
  "PFMPA_Results","PFNA_Results","PFNS_Results", "PFOA_Results","PFOS_Results",
  "PFOSA_Results", "PFPeA_Results", "PFPeS_Results","PFTeDA_Results",
  "PFTrDA_Results","PFUnA_Results","PFCA", "PFSA","FTSA","PFOSA", "FOSAA",
  "PFECA","PFESA", "FTCA")

##  4.1.  Observation sites  

#Read in the observation data.   

#r read in observation covariate data}
obs_covariates <- ssn_get_data(ssn_obj, name = "obs") |> 
  dplyr::select(rid, pid, ratio, snapdist, upDist, afvArea, locID, netID,                   
                site_id, reach_id, ssnid, 
                all_of(pfas_compounds), netgeom, geometry) |>
  left_join(reach_data, by = "reach_id")

#Some observations are outside the typical range experienced across a majority of the system. These can cause model fitting issues and therefore we may consider removing those points for the model fitting process. We can predict those values to see if their extreme highs are addressed by the model and we can always add them back in and refit.  

if(isTRUE(remove_outliers)) {
  obs_covariates$PFAS40[obs_covariates$ssnid == 6] <- NA
  obs_covariates$PFCA[obs_covariates$ssnid == 6] <- NA
  obs_covariates$PFBA_Results[obs_covariates$ssnid == 6] <- NA
  obs_covariates$PFBS_Results[obs_covariates$ssnid == 6] <- NA
  obs_covariates$PFHpA_Results[obs_covariates$ssnid == 6] <- NA
  obs_covariates$PFHxA_Results[obs_covariates$ssnid == 6] <- NA
  obs_covariates$PFOA_Results[obs_covariates$ssnid == 6] <- NA
  obs_covariates$PFPeA_Results[obs_covariates$ssnid == 6] <- NA
}

# ---- Create log1p response variables --------------------------------
# log(x + 1) = log1p(x) in R. Back-transform with expm1().
# Applied AFTER remove_outliers so NA values propagate correctly.
# Applied to all pfas_compounds; creates parallel set of response columns.

log_pfas_compounds <- paste0("log1p_", pfas_compounds)

obs_covariates <- obs_covariates |>
  mutate(across(
    .cols = all_of(pfas_compounds),
    .fns  = ~log1p(.),
    .names = "log1p_{.col}"
  ))

# Sanity check — confirm no negative values snuck in
neg_check <- obs_covariates |>
  st_drop_geometry() |>
  select(all_of(pfas_compounds)) |>
  summarise(across(everything(), ~sum(. < 0, na.rm = TRUE)))

if (any(neg_check > 0)) {
  warning("Negative values found in raw PFAS columns — log1p will produce NaN. Check these compounds:")
  print(names(neg_check)[neg_check > 0])
}

cat("Raw PFAS columns:    ", length(pfas_compounds), "\n")
cat("Log1p PFAS columns:  ", length(log_pfas_compounds), "\n")
cat("Total response vars: ", length(pfas_compounds) + length(log_pfas_compounds), "\n")
#Standardize the continuous variables used in the model fit.  

#r standardize_continuous_variables}
# standardize continuous variables so their estimates can be compared directly
continuous_vars <- obs_covariates |> st_drop_geometry() |> 
  select(all_of(names(reach_data)), -reach_id, -is_lentic) |> names()

# grab only to continuous covariates we want to use in our modeling
continuous <- obs_covariates[ ,continuous_vars]

# save continuous variables to file for standardizing prediction data later
saveRDS(continuous, file = "../data/PFAS_obs_continuous_df.RDS")

# standardize the continuous variables
cont_s <- continuous |> modify_at(continuous_vars, stand) |>
  rename_with(.fn = ~ paste0(.x, "_s"), .cols = everything())

# join raw and standardize observation site data into one data frame
obs_covariates_s <- data.frame(obs_covariates, cont_s) |> st_as_sf()

# write out the raw & standardized observation input data to be used later
saveRDS(obs_covariates_s,
        file = "../data/PFAS_obs_continuous_df_standardized.RDS")

# input the observation data into a new SSN object
ssn_obj_obs_std <- ssn_put_data(obs_covariates_s, ssn_obj, 
                                name = "obs", resize_data = FALSE)


#View the observation data.  

DT::datatable(
  obs_covariates_s |> 
    mutate(across(where(is.numeric), \(x) round(x,3))) |> 
    column_to_rownames(var = "site_id"),
  extensions = "FixedColumns",
  options = list(
    scrollX = TRUE,
    fixedColumns = list(leftColumns = 1)
  ),
  class = "display nowrap compact"
)



## 4.2. Prediction sites  

#Read in prediction data contained in the SSN object and join to the reach-based covariate data set.  

#r read in prediction covariate data}
preds_covariates <- ssn_get_data(ssn_obj_obs_std, name = "preds") |> 
  dplyr::select(rid, pid, ratio, snapdist, upDist, afvArea, locID, netID,                   
                reach_id, ssnid, netgeom, geom) |>
  left_join(reach_data, by = "reach_id")
# names(preds_covariates)


#Standardize the prediction data based on the mean and SD of the observation data being used to fit the SSN model.  

# extract the continuous data columns
cont_cov_data_preds <- preds_covariates |> st_drop_geometry() |>
  dplyr::select(all_of(continuous_vars))

# use the "stdpreds()" function to standardize each prediction covariate
# based on the observation data for each covariate in the fitted model
cont_cov_data_preds_stdzd <- 
  stdpreds(preds_cont_var_df = st_drop_geometry(cont_cov_data_preds), 
           obs_cont_var_df = st_drop_geometry(continuous) )


# add prefix ("s_") to standardized column names
colnames(cont_cov_data_preds_stdzd) <- paste0(names(st_drop_geometry(cont_cov_data_preds)),"_s")

# join standardized data to original data frame
preds_covariates_std <- 
  bind_cols(preds_covariates, cont_cov_data_preds_stdzd) %>%
  st_as_sf()

saveRDS(preds_covariates_std, "../data/PFAS_preds_covariates_df_standardized.RDS")

# put the data frame with standardized preds back into the ssn object
ssn_obj_obspreds_std <- 
  ssn_put_data(preds_covariates_std, ssn_obj_obs_std, 
               name = "preds", resize_data = FALSE)
# summary(preds_covariates_std)


# ---- 5. Model fitting and selection ----

# Extract the continuous covariate names directly from reach_data
available_continuous <- names(reach_data)[names(reach_data) != "reach_id"]

# Build covariate groups from actual column names
other_cov <- c("basin_area_km2")

npdes_cov <- c("npdes_density_km2_bas", "npdes_count_bas", 
               "npdes_density_km2_rca", "npdes_count_rca")

impervious_cov <- c("impervious_2025_bas", "impervious_2025_rca")

agriculture_cov <- c("pct_crops_2025_bas", "pct_crops_2025_rca",
                     "pct_pasture_hay_2025_bas", "pct_pasture_hay_2025_rca",
                     "pct_crops_and_pasture_2025_bas", "pct_crops_and_pasture_2025_rca")

baseflow_cov <- c("bfi_bas", "bfi_rca")

# Soil properties 
soil_cov <- c("clay_pct_bas", "clay_pct_rca",
              "om_bas", "om_rca",
              "ksat_bas", "ksat_rca")

# Build the formula list by starting with the "other" covariates
#ssn_formula_list0 <- paste(" ~", formula_all_combo_fnx(other_cov, 1))

# ---- 5. Model fitting and selection ----

available_continuous <- names(reach_data)[names(reach_data) != "reach_id"]

# DELETED: other_cov <- c("basin_area_km2") 

npdes_cov <- c("npdes_density_km2_bas", "npdes_count_bas",
               "npdes_density_km2_rca", "npdes_count_rca")

impervious_cov <- c("impervious_2025_bas", "impervious_2025_rca")

agriculture_cov <- c("pct_crops_2025_bas", "pct_crops_2025_rca",
                     "pct_pasture_hay_2025_bas", "pct_pasture_hay_2025_rca",
                     "pct_crops_and_pasture_2025_bas", "pct_crops_and_pasture_2025_rca")

baseflow_cov <- c("bfi_bas", "bfi_rca")

soil_cov <- c("clay_pct_bas", "clay_pct_rca",
              "om_bas", "om_rca",
              "ksat_bas", "ksat_rca")

# [CHANGED] No longer built on top of ssn_formula_list0.
# npdes is now the new starting group — standalone formulas only.
ssn_formula_list1 <- paste(" ~", formula_all_combo_fnx(npdes_cov, 1))  # ← CHANGED

# Everything below is identical to before
ssn_formula_list2 <- c(
  paste(ssn_formula_list1, "+", rep(formula_all_combo_fnx(impervious_cov, 1),
                                    each = length(ssn_formula_list1))),
  paste(" ~", formula_all_combo_fnx(impervious_cov, 1))
)

ssn_formula_list3 <- c(
  paste(ssn_formula_list2, "+", rep(formula_all_combo_fnx(agriculture_cov, 1),
                                    each = length(ssn_formula_list2))),
  paste(" ~", formula_all_combo_fnx(agriculture_cov, 1))
)

ssn_formula_list4 <- c(
  paste(ssn_formula_list3, "+", rep(formula_all_combo_fnx(baseflow_cov, 1),
                                    each = length(ssn_formula_list3))),
  paste(" ~", formula_all_combo_fnx(baseflow_cov, 1))
)

ssn_formula_list5 <- c(
  paste(ssn_formula_list4, "+", rep(formula_all_combo_fnx(soil_cov, 1),
                                    each = length(ssn_formula_list4))),
  paste(" ~", formula_all_combo_fnx(soil_cov, 1))
)

ssn_formula_list6 <- c(
  paste(ssn_formula_list5, "+", rep(formula_all_combo_fnx(c("is_lentic"), 1),
                                    each = length(ssn_formula_list5))),
  paste(" ~", formula_all_combo_fnx(c("is_lentic"), 1))
)
# Clean up formula list and filter to 2-3 covariates only
ssn_formula_list <- gsub(" + 1", "", ssn_formula_list6, fixed = TRUE) |>
  gsub(" 1 +", "", x = _, fixed = TRUE) |>
  unique() |>
  (\(x) x[order(nchar(x))])() |>
  (\(x) x[stringr::str_count(x, "\\+") %in% c(1, 2)])()

pander(paste("total number of model structures to fit:", length(ssn_formula_list)))

##  5.2.  Fit competing models  

#Identify the key compounds to use as response variables from the various PFAS families/compounds.  

# These are the main compounds for response variables 
# After model fitting we compare on the SAME scale (back-transformed RMSE)
# AIC is NOT comparable between log and linear models — different likelihoods

all_response_vars <- c(pfas_compounds, log_pfas_compounds)
key_compounds     <- all_response_vars

cat("Total model sets to fit:", length(key_compounds), "\n")
cat("  Raw:   ", length(pfas_compounds), "compounds\n")
cat("  Log1p: ", length(log_pfas_compounds), "compounds\n")
cat("  Formulas per compound:", length(ssn_formula_list), "\n")
cat("  Total model fits:", length(key_compounds) * length(ssn_formula_list), "\n")


#We nest the model fit within a diagnostic statistic function to provide the fit and the diagnostics for quick model comparison via AIC.  


# If cl exists from a previous attempt, stop it first
if (exists("cl")) {
  try(stopCluster(cl), silent = TRUE)
}
# Also good practice: check for and kill any zombie Rscript workers
# at the OS level if this keeps happening (Task Manager on Windows)

cl <- makeCluster(4)
registerDoParallel(cl)

parallel::clusterExport(cl, 
                        varlist = c("ssn_formula_list", "ssn_obj_obspreds_std"),
                        envir = environment())
parallel::clusterEvalQ(cl, library("SSN2")) # make sure the SSN2 library is loaded on each cluster

# Recreate distance matrices with updated data
ssn_create_distmat(ssn_obj_obspreds_std, predpts = "preds", overwrite = TRUE)



# run parallel loop and save results into the following data frame
{
  t0 <- Sys.time()
  times_per_compound <- c()   # track actual pace as it runs
  
  PFAS_ssn_fits <-
    foreach(i          = 1:length(key_compounds),
            .combine   = "rbind",
            .packages  = c("SSN2", "dplyr", "purrr", "stringr")) %dopar% {
              
              formula_list <-
                lapply(paste(key_compounds[[i]], ssn_formula_list), as.formula)
              
              purrr::map(.x = formula_list,
                         .f = function(x) {
                           SSN2::glance(SSN2::ssn_lm(
                             formula     = x,
                             ssn.object  = ssn_obj_obspreds_std,
                             tailup_type = "exponential",
                             additive    = "afvArea",
                             estmethod   = "ml")) |>
                             dplyr::mutate(
                               response_var  = stringr::str_split_i(deparse1(x), " ~ ", 1),
                               predictor_var = stringr::str_split_i(deparse1(x), " ~ ", 2))
                         }) |> purrr::list_rbind()
            }
  
  stopCluster(cl)
  
  t1  <- Sys.time()
  e0  <- t1 - t0
  
  pander(paste(length(key_compounds) * length(ssn_formula_list),
               "models fit in:", round(e0, 2), attr(e0, "units")))
}

## 5.3. Select best models based on AIC for each compound {.tabset}  

#We use the delta AIC value of 2 to retain models that have comparable fit to the observed data. At the same time we also calculate the AIC weight of each model among the retained models for a given response variable or PFAS compound/family.  

#According to Burnham and Anderson (2002), the AIC weights for a set of competing models can be calculated as:  
#  $$w_i = \frac{exp(-\frac{1}{2}\Delta_i)}{\sum_{r=1}^R exp(-\frac{1}{2}\Delta_r)}$$  
#   where $w_i$ is the AIC weight for the *i^th^* model relative to the model with the lowest AIC value, *R* is the total number of models in the set being considered for the AIC weighting, $exp(-\frac{1}{2}\Delta_i)$ is the Likelihood of model *i* and $\sum_{r=1}^R exp(-\frac{1}{2}\Delta_r)$ is the sum of all likelihoods in the set of models being considered.  

#r select best models based on AIC}
best_PFAS_models <- PFAS_ssn_fits |>
  group_by(response_var) |>
  group_split() |>
  map(.x = _,
      .f = function(x) {
        x |>
          arrange(AIC) |>
          filter(AIC < min(AIC) + 2) |>
          rownames_to_column("rank") |> 
          dplyr::mutate(delta_aic = AIC - min(AIC, na.rm = TRUE)) |>
          dplyr::mutate(likelihood_i = exp(-0.5*delta_aic)) |>
          dplyr::mutate(AIC_weight_i = likelihood_i/sum(likelihood_i, na.rm = TRUE)) |>
          dplyr::relocate(response_var)
      }) 
names(best_PFAS_models) <- PFAS_ssn_fits |>
  group_by(response_var) |>
  group_keys() |> as.list() |> unlist()

DT::datatable(best_PFAS_models |> bind_rows() |> mutate(across(where(is.numeric),~round(.x, 2))))

saveRDS(best_PFAS_models |> bind_rows(), file = "../outputs/PFAS_best_model_fit_diagnostics_table.RDS")


#Each row in these tables represents a model retained based on its delta AIC value (<= 2 AIC units from minimum AIC).  

#r, results='asis', echo=FALSE}
for (resp_var in names(best_PFAS_models)) {
  
  model <- best_PFAS_models[[resp_var]]
  
  # Tab header
  cat("\n### ", resp_var, "{.tabset}  \n\n", sep = "")
  
  tab_content <- tagList(
    datatable(model |> mutate(across(where(is.numeric), \(x) round(x, 2))),
              extensions = "FixedColumns",
              options = list(
                scrollX = TRUE,
                fixedColumns = list(leftColumns = 1)
              ),
              class = "display nowrap compact")
  )
  
  print(tab_content)
  cat("\n")
}

## 5.4.  Best fit models fit individually   

#Fit only the models retained by AIC comparison from above for model evaluation and diagnostics.  

#r fit suite of best models, cache = TRUE}
# create a named list of formula from best model suites
best_model_formla_ls <- 
  map(.x = best_PFAS_models,
      .f = function(x,y){
        form_list <- x |>
          mutate(formula = paste(response_var, "~", predictor_var)) |>
          select(formula) |> as.list() |> # creates a list object, but only 1 element
          unlist() |> as.list() |> # breaks list element into individual pieces
          lapply(FUN = as.formula) # converts to formula for each element
        
        names(form_list) <- x  |>
          mutate(model_rank = paste0(response_var, "_", rank)) |>
          select(model_rank) |> as.list() |> unlist()
        
        form_list
      })

fit_model_ls <- 
  map2(.x = best_model_formla_ls,
       .y = names(best_model_formla_ls),
       .f = function(x,y) {
         # create an empty list to hold each SSN model fit
         loop_model_ls <- as.vector(rep(NA,length(x)),mode="list")
         
         names(loop_model_ls) <- names(x)
         
         for (i in 1:length(x)) {
           loop_model_ls[[i]] <- ssn_lm(x[[i]],
                                        ssn.object = ssn_obj_obspreds_std,
                                        tailup_type = "exponential",
                                        additive = "afvArea",
                                        estmethod = "ml")
         }
         loop_model_ls
       })

# 6. Model evaluation and diagnostics


## 6.0. Individual model summaries {.tabset}
# Calculate the model diagnostics for each individual model in each suite of
# models for a PFAS compound/family. Runs on ALL compounds (both raw and
# log1p scales) before scale comparison in Section 6.1. Do not filter here.

library(htmltools)

indv_model_summary_ls <-
  map2(.x = fit_model_ls,
       .y = best_PFAS_models,
       .f = function(x, y) {
         modsuite_glance_df   <- data.frame()
         modsuite_loocv_df    <- data.frame()
         modsuite_varcomp_df  <- data.frame()
         modsuite_esttable_df <- data.frame()
         
         for (i in 1:length(x)) {
           loop_glance_df <- glance(x[[i]]) |>
             mutate(modelrank    = i,
                    response_var = y$response_var[i],
                    modelformula = y$predictor_var[i])
           modsuite_glance_df <- bind_rows(modsuite_glance_df, loop_glance_df)
           
           loop_loocv_df <- loocv(x[[i]]) |>
             mutate(modelrank    = i,
                    response_var = y$response_var[i],
                    modelformula = y$predictor_var[i])
           modsuite_loocv_df <- bind_rows(modsuite_loocv_df, loop_loocv_df)
           
           loop_varcomp_df <- SSN2::varcomp(x[[i]]) |>
             mutate(modelrank    = i,
                    response_var = y$response_var[i],
                    modelformula = y$predictor_var[i])
           modsuite_varcomp_df <- bind_rows(modsuite_varcomp_df, loop_varcomp_df)
           
           ssn_tidy_out      <- tidy(x[[i]], conf.int = TRUE)
           continuous_nogeom <- st_drop_geometry(continuous)
           esttable <- std_to_raw_estimate_table_fnx(continuous_nogeom,
                                                     ssn_tidy_out, 5) |>
             mutate(modelrank    = i,
                    response_var = y$response_var[i],
                    modelformula = y$predictor_var[i])
           modsuite_esttable_df <- bind_rows(modsuite_esttable_df, esttable)
         }
         
         write_csv(modsuite_glance_df,
                   file = paste0("../outputs/model_summ_tables/",
                                 y$response_var[i], "_modsuite_glance.csv"))
         write_csv(modsuite_loocv_df,
                   file = paste0("../outputs/model_summ_tables/",
                                 y$response_var[i], "_modsuite_loocv.csv"))
         write_csv(modsuite_varcomp_df,
                   file = paste0("../outputs/model_summ_tables/",
                                 y$response_var[i], "_modsuite_varcomp.csv"))
         write_csv(modsuite_esttable_df,
                   file = paste0("../outputs/model_summ_tables/",
                                 y$response_var[i], "_modsuite_esttable.csv"))
         
         list(modsuite_glance_df    = modsuite_glance_df,
              modsuite_loocv_df     = modsuite_loocv_df,
              modsuite_varcomp_df   = modsuite_varcomp_df,
              modsuite_esttable_df  = modsuite_esttable_df)
       })


## 6.1. Raw vs log1p scale comparison {.tabset}
# AIC cannot compare log vs linear models - different likelihoods, different
# scales. We use AIC-weighted LOOCV RMSPE from Section 6.0 (true out-of-sample
# performance) as the comparison metric.
#
# Raw RMSPE is in ng/L; log1p RMSPE is in log(ng/L+1). Direct comparison is
# invalid. We normalize each RMSPE by the mean observed value on its own scale
# to produce a dimensionless CV-RMSPE, comparable across both scales and all
# 49 compounds.
#
# Compounds with less than 5% CV-RMSPE difference default to raw scale for
# interpretability since ng/L units are directly reportable.
#
# Output: final_key_compounds - one winning response variable per compound.
# All downstream sections 6.2 onward are filtered to this list.

scale_comparison_ls <- list()

for (compound in pfas_compounds) {
  
  log_compound <- paste0("log1p_", compound)
  
  has_raw <- compound     %in% names(indv_model_summary_ls)
  has_log <- log_compound %in% names(indv_model_summary_ls)
  
  if (!has_raw || !has_log) {
    warning("Missing model summaries for ", compound, " - skipping.")
    next
  }
  
  # ---- AIC-weighted LOOCV RMSPE — raw model (ng/L) ----------------------
  loocv_raw <- indv_model_summary_ls[[compound]]$modsuite_loocv_df |>
    left_join(
      best_PFAS_models[[compound]] |>
        select(rank, AIC_weight_i) |>
        mutate(rank = as.integer(rank)),
      by = c("modelrank" = "rank")
    ) |>
    summarise(
      rmspe_wtd = sum(RMSPE * AIC_weight_i, na.rm = TRUE),
      bias_wtd  = sum(bias  * AIC_weight_i, na.rm = TRUE)
    )
  
  # ---- AIC-weighted LOOCV RMSPE — log1p model (log(ng/L+1)) ------------
  loocv_log <- indv_model_summary_ls[[log_compound]]$modsuite_loocv_df |>
    left_join(
      best_PFAS_models[[log_compound]] |>
        select(rank, AIC_weight_i) |>
        mutate(rank = as.integer(rank)),
      by = c("modelrank" = "rank")
    ) |>
    summarise(
      rmspe_wtd = sum(RMSPE * AIC_weight_i, na.rm = TRUE),
      bias_wtd  = sum(bias  * AIC_weight_i, na.rm = TRUE)
    )
  
  # ---- Back-transformed RMSE for log1p model in ng/L -------------------
  # loocv() only returns aggregate stats — individual LOO predictions are
  # not exposed by SSN2. Instead we use augment() (in-sample) predictions,
  # back-transform to ng/L, then compute AIC-weighted RMSE.
  # Labeled as "RMSE" (not LOOCV RMSPE) to be transparent about this.
  
  log_models  <- fit_model_ls[[log_compound]]
  log_weights <- best_PFAS_models[[log_compound]] |>
    mutate(rank = as.integer(rank)) |>
    arrange(rank) |>
    pull(AIC_weight_i)
  
  rmspe_log_ng_per_L <- map2_dbl(
    log_models,
    log_weights,
    function(m, w) {
      aug     <- augment(m) |> st_drop_geometry()
      obs_ng  <- expm1(aug[[log_compound]])  # back-transform observed
      pred_ng <- expm1(aug$.fitted)           # back-transform predicted
      rmse_i  <- sqrt(mean((obs_ng - pred_ng)^2, na.rm = TRUE))
      rmse_i * w                              # weight by AIC weight
    }
  ) |> sum()
  
  # ---- CV-RMSPE (dimensionless, for scale comparison) ------------------
  obs_mean_raw <- obs_covariates_s |> st_drop_geometry() |>
    pull(all_of(compound)) |> mean(na.rm = TRUE)
  
  obs_mean_log <- obs_covariates_s |> st_drop_geometry() |>
    pull(all_of(log_compound)) |> mean(na.rm = TRUE)
  
  cv_rmspe_raw <- loocv_raw$rmspe_wtd / obs_mean_raw
  cv_rmspe_log <- loocv_log$rmspe_wtd / obs_mean_log
  
  pct_diff <- round(100 * (cv_rmspe_raw - cv_rmspe_log) / cv_rmspe_raw, 1)
  better   <- ifelse(cv_rmspe_log < cv_rmspe_raw, "log1p", "raw")
  
  scale_comparison_ls[[compound]] <- data.frame(
    compound               = compound,
    rmspe_raw              = round(loocv_raw$rmspe_wtd,  3),  # LOOCV, ng/L
    rmspe_log1p            = round(loocv_log$rmspe_wtd,  3),  # LOOCV, log(ng/L+1)
    rmspe_log1p_ng_per_L   = round(rmspe_log_ng_per_L,   3),  # in-sample RMSE, ng/L
    cv_rmspe_raw           = round(cv_rmspe_raw,          4),
    cv_rmspe_log1p         = round(cv_rmspe_log,          4),
    bias_raw               = round(loocv_raw$bias_wtd,    3),
    bias_log1p             = round(loocv_log$bias_wtd,    3),
    better_scale           = better,
    cv_rmspe_pct_diff      = pct_diff
  )
}

scale_comparison_df <- bind_rows(scale_comparison_ls) |>
  arrange(desc(abs(cv_rmspe_pct_diff)))

saveRDS(scale_comparison_df, "../outputs/PFAS_scale_comparison_raw_vs_log1p.RDS")

DT::datatable(scale_comparison_df |>
                mutate(across(where(is.numeric), ~round(.x, 3))))

cat("\n=== Scale comparison summary ===\n")
cat("  Compounds better as log1p:    ",
    sum(scale_comparison_df$better_scale == "log1p"), "\n")
cat("  Compounds better as raw:      ",
    sum(scale_comparison_df$better_scale == "raw"), "\n")
cat("  Compounds within 5% CV-RMSPE:",
    sum(abs(scale_comparison_df$cv_rmspe_pct_diff) < 5),
    " (defaulting to raw for interpretability)\n")

# Build final lookup: one winning response variable per compound.
# Compounds within 5% CV-RMSPE difference are assigned raw scale.
best_scale_lookup <- scale_comparison_df |>
  mutate(
    better_scale = ifelse(abs(cv_rmspe_pct_diff) < 5, "raw", better_scale)
  ) |>
  transmute(
    compound,
    best_response_var = ifelse(
      better_scale == "log1p",
      paste0("log1p_", compound),
      compound
    ),
    better_scale,
    cv_rmspe_pct_diff
  )

final_key_compounds <- best_scale_lookup$best_response_var

cat("\nFinal response variables selected per compound:\n")
print(best_scale_lookup)


## 6.2. Relative importance of covariates {.tabset}
# Burnham and Anderson (2002): relative importance is the sum of AIC weights
# across all retained models that include that variable.
# [CHANGED] best_PFAS_models filtered to final_key_compounds from Section 6.1.

rel_imp_df_ls <- best_PFAS_models[final_key_compounds] |>   # [CHANGED]
  map(.x = _,
      .f = function(x) {
        x |> group_by(rank) |> group_split() |>
          map(.x = _,
              .f = function(x){
                model_df <-
                  data.frame(covariate = strsplit(x = x$predictor_var,
                                                  split = " + ",
                                                  fixed = TRUE)[[1]],
                             presence = 1)
                colnames(model_df) <- c("covariate",
                                        paste0(x$response_var,"_",x$rank))
                model_df
              }) |>
          purrr::reduce(full_join) |>
          column_to_rownames("covariate") |>
          t()  %>%
          replace(is.na(.), 0) |>
          data.frame() |>
          rownames_to_column("model_name") %>% data.frame() |>
          left_join(x |> mutate(model_name = paste0(x$response_var,"_",x$rank)),
                    by = c("model_name") ) |>
          select(-any_of(c("AIC", "AICc", "BIC", "delta_aic", "likelihood_i",
                           "n", "p","npar", "value", "logLik", "deviance",
                           "pseudo.r.squared"))) %>%
          mutate_at(vars(-model_name,-AIC_weight_i,-response_var,
                         -rank, -predictor_var),
                    function(x){x*.$AIC_weight_i}) |>
          select(-model_name,-AIC_weight_i,-response_var, -rank, -predictor_var) |>
          summarise_all(sum, na.rm = TRUE) |>
          pivot_longer(cols = everything()) |>
          rename(covariate = name, rel_imp = value) |>
          arrange(desc(rel_imp)) |>
          data.frame() |> mutate(response_var = x$response_var[1]) |>
          relocate(response_var, covariate, rel_imp)
      })

# Relative importance plot. Covariates ordered top-to-bottom by mean rel. importance.
rel_imp_df_plotting <- list_rbind(rel_imp_df_ls) |>
  pivot_wider(id_cols = "covariate", values_from = rel_imp,
              names_from = response_var, values_fill = 0) |>
  mutate(overall_mean = rowMeans(across(where(is.numeric)), na.rm = TRUE)) |>
  mutate(covariate = forcats::fct_reorder(covariate, overall_mean, .desc = FALSE)) |>
  pivot_longer(cols = !covariate, names_to = "response_var", values_to = "rel_imp")

color_vector <- c("#4daf4a","#e41a1c","#984ea3","#ff7f00","#377eb8",
                  "blue", "red", "orange", "maroon", "purple")

(relimp_point_plot <-
    ggplot(data = rel_imp_df_plotting,
           aes(x = covariate,
               fill = response_var, shape = response_var, size = response_var,
               y = rel_imp)) +
    geom_point() +
    theme_bw() +
    labs(y     = "Relative importance\n\n\n\n",
         fill  = "Response variable",
         shape = "Response variable",
         size  = "Response variable",
         x     = "Covariate") +
    scale_fill_manual(values  = c("black",color_vector)) +
    scale_shape_manual(values = c(3,21,24,22,25,23)) +
    scale_size_manual(values  = c(4,3.5,2.25,2.25,2.25,2.25)) +
    coord_cartesian(xlim = c(0,1)) +
    theme(legend.position   = c(0.05,-0.2),
          legend.title      = element_text(size = 12),
          legend.text       = element_text(size = 10),
          legend.background = element_blank(),
          axis.text.y       = ggtext::element_markdown(),
          legend.key        = element_blank()) +
    guides(fill=guide_legend(nrow=2, byrow=FALSE)) +
    coord_flip()
)

ggsave(relimp_point_plot,
       filename="../figs/Figure00_Relative_Importance_plot.png",
       width = 6.5, height = 4.5, units = "in", dpi = 600)

for (resp_var in names(rel_imp_df_ls)) {
  model <- rel_imp_df_ls[[resp_var]]
  cat("\n### ", resp_var, " \n\n", sep = "")
  tab_content <- tagList(
    datatable(model |> mutate(across(where(is.numeric), \(x) round(x, 2))),
              extensions = "FixedColumns",
              options = list(scrollX = TRUE,
                             fixedColumns = list(leftColumns = 1)),
              class = "display nowrap compact")
  )
  print(tab_content)
  cat("\n")
}


## 6.3. Averaged model suite summaries {.tabset}
# AIC-weighted mean diagnostics for each compound's best model suite.
# [CHANGED] Both map2() inputs filtered to final_key_compounds from Section 6.1.

mean_model_summary_ls <-
  map2(.x = fit_model_ls[final_key_compounds],          # [CHANGED]
       .y = indv_model_summary_ls[final_key_compounds], # [CHANGED]
       .f = function(x, y) {
         
         avg_esttable_df <-
           map2(.x = x,
                .y = names(x),
                .f = ~ {
                  std_to_raw_estimate_table_fnx(
                    continuous_vars_df = st_drop_geometry(continuous),
                    tidy_out_tibble    = tidy(.x, conf.int = TRUE),
                    roundval = 5,
                    message = FALSE) |>
                    mutate(model_name = .y)
                }) |>
           bind_rows() |>
           reframe(.by = term,
                   model_n    = n(),
                   std_est    = mean(std_est,    na.rm = TRUE),
                   std_est_sd = sd(std_est,      na.rm = TRUE),
                   raw_est    = mean(raw_est,    na.rm = TRUE),
                   raw_est_sd = sd(raw_est,      na.rm = TRUE),
                   t_stat     = mean(t_stat,     na.rm = TRUE),
                   p_val      = mean(p_val,      na.rm = TRUE)) |>
           arrange(desc(model_n)) |>
           rename(covariate = "term")
         
         write_csv(avg_esttable_df,
                   file = paste0("../outputs/model_summ_tables/",
                                 y$modsuite_glance_df$response_var[1],
                                 "_modsuite_esttable_MEANS.csv"))
         
         avg_modsuite_loocv_df <-
           bind_rows(
             summarize(y$modsuite_loocv_df,
                       across(-starts_with(c("response_var", "modelformula")),
                              .fns = function(q){ mean(q, na.rm = TRUE) })) |>
               mutate(statistic = "mean"),
             summarize(y$modsuite_loocv_df,
                       across(-starts_with(c("response_var", "model")),
                              .fns = function(p) {sd(p, na.rm = TRUE) })) |>
               mutate(statistic = "sd")) |>
           select(statistic, everything(), -modelrank) |>
           pivot_longer(cols = -statistic,
                        names_to = "LOOCV Statistic", values_to = "values") |>
           pivot_wider(names_from = statistic, values_from = values)
         
         write_csv(avg_modsuite_loocv_df,
                   file = paste0("../outputs/model_summ_tables/",
                                 y$modsuite_glance_df$response_var[1],
                                 "_modsuite_loocv_MEANS.csv"))
         
         avg_modsuite_varcomp_df <- y$modsuite_varcomp_df |>
           reframe(.by = varcomp,
                   prop_mean = mean(proportion, na.rm = TRUE),
                   prop_sd   = sd(proportion,   na.rm = TRUE),
                   prop_min  = min(proportion,   na.rm = TRUE),
                   prop_max  = max(proportion,   na.rm = TRUE)) |>
           rename(variance_component = "varcomp")
         
         write_csv(avg_modsuite_varcomp_df,
                   file = paste0("../outputs/model_summ_tables/",
                                 y$modsuite_glance_df$response_var[1],
                                 "_modsuite_varcomp_MEANS.csv"))
         
         list(avg_esttable_df         = avg_esttable_df,
              avg_modsuite_loocv_df   = avg_modsuite_loocv_df,
              avg_modsuite_varcomp_df = avg_modsuite_varcomp_df)
       })

for (resp_var in names(mean_model_summary_ls)) {
  model <- mean_model_summary_ls[[resp_var]]
  cat("\n### ", resp_var, "{.tabset}  \n\n", sep = "")
  tab_content <- tagList(
    tags$h4("LOOCV"),
    datatable(model$avg_modsuite_loocv_df |>
                mutate(across(where(is.numeric), \(x) round(x, 3))),
              extensions = "FixedColumns",
              options = list(scrollX = TRUE,
                             fixedColumns = list(leftColumns = 1)),
              class = "display nowrap compact"),
    tags$h4("Variance components"),
    datatable(model$avg_modsuite_varcomp_df |>
                mutate(across(where(is.numeric), \(x) round(x, 3))),
              extensions = "FixedColumns",
              options = list(scrollX = TRUE,
                             fixedColumns = list(leftColumns = 1)),
              class = "display nowrap compact"),
    tags$h4("Estimates"),
    datatable(model$avg_esttable_df |>
                mutate(across(where(is.numeric), \(x) round(x, 3))),
              extensions = "FixedColumns",
              options = list(scrollX = TRUE,
                             fixedColumns = list(leftColumns = 1)),
              class = "display nowrap compact")
  )
  print(tab_content)
  cat("\n")
}


## 6.4. Model diagnostics plots {.tabset}
# Obs/pred, residuals vs fitted, Q-Q, and Cook's Distance plots.
# [CHANGED] Both map2() inputs filtered to final_key_compounds from Section 6.1.

diagnostics_plot_ls <-
  map2(.x = fit_model_ls[final_key_compounds],       # [CHANGED]
       .y = best_PFAS_models[final_key_compounds],   # [CHANGED]
       .f = function(x, y) {
         
         obs_df <- augment(x[[1]]) |>
           select(pid, observations = 1) |>
           st_drop_geometry()
         
         all_model_obspreds_df <-
           map2(.x = x,
                .y = names(x),
                .f = function (x,y) {
                  SSN2::augment(x = x) |>
                    dplyr::select(pid, observations = 1, starts_with(".")) %>%
                    dplyr::mutate(model_name = y) %>%
                    st_drop_geometry() %>% data.frame()
                }) |>
           bind_rows() |>
           rename(predictions = .fitted) |>
           mutate(pred_type = "Individual model predictions",
                  pid = as.numeric(pid))
         
         avg_model_obspreds_df <-
           left_join(x = all_model_obspreds_df,
                     y = y |> mutate(model_name = paste0(response_var,"_",rank)),
                     by = "model_name") %>%
           dplyr::mutate(weighted_pred     = predictions * AIC_weight_i,
                         weighted_resid    = .resid      * AIC_weight_i,
                         weighted_hat      = .hat        * AIC_weight_i,
                         weighted_cooksd   = .cooksd     * AIC_weight_i,
                         weighted_stdresid = .std.resid  * AIC_weight_i) %>%
           reframe(.by = c(pid),
                   observations = mean(observations, na.rm = TRUE),
                   predictions  = sum(weighted_pred),
                   .resid       = sum(weighted_resid),
                   .hat         = sum(weighted_hat),
                   .cooksd      = sum(weighted_cooksd),
                   .std.resid   = sum(weighted_stdresid),
                   model_name   = "Model_Average") |>
           mutate(pred_type = "AIC-weighted predictions",
                  pid = as.numeric(pid))
         
         obspreds_df <- bind_rows(all_model_obspreds_df, avg_model_obspreds_df) |>
           left_join(select(obs_covariates, pid, site_id), by = "pid")
         
         obspred_plot <-
           ggplot(data = obspreds_df, aes(x = predictions, y = observations)) +
           geom_point() +
           theme_minimal() +
           geom_abline(intercept = 0, slope = 1) +
           labs(y     = "Observed values",
                title = paste0(y$response_var,": Observed vs. Predicted"),
                x     = "Predicted values") +
           facet_wrap(~pred_type, ncol = 2)
         
         residfit_plot <-
           ggplot(data = obspreds_df, aes(x = predictions, y = .resid)) +
           geom_point(alpha = 0.5) +
           geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
           geom_smooth(method = "loess", color = "blue", se = FALSE) +
           theme_minimal() +
           labs(x     = "Fitted",
                title = paste0(y$response_var,": Residuals vs. Fitted"),
                y     = "Residuals") +
           facet_wrap(~pred_type, ncol = 2)
         
         qq_plot <-
           ggplot(data = obspreds_df, aes(sample = .resid)) +
           stat_qq() +
           stat_qq_line(color = "red") +
           theme_minimal() +
           labs(title = paste0(y$response_var,": Normal Q-Q Plot")) +
           facet_wrap(~pred_type, ncol = 2)
         
         threshold <- 4 / nrow(avg_model_obspreds_df)
         cooksd_plot <-
           ggplot(obspreds_df, aes(x = site_id, y = .cooksd)) +
           geom_hline(yintercept = threshold, linetype = "dashed",
                      color = "red", linewidth = 0.8) +
           geom_segment(aes(x = site_id, xend = site_id, y = 0, yend = .cooksd),
                        color = "gray60", alpha = 0.7) +
           geom_point(aes(color = .cooksd > threshold), size = 2.5, alpha = 0.5,
                      show.legend = FALSE) +
           scale_color_manual(values = c("FALSE" = "black", "TRUE" = "firebrick")) +
           theme_minimal(base_size = 13) +
           labs(title    = paste0(y$response_var,": Cook's Distance"),
                subtitle = paste0("Red line indicates traditional threshold (4/n = ",
                                  round(threshold, 4), ")"),
                x = "Observation Site",
                y = "Weighted Cook's Distance (D)") +
           theme(panel.grid.minor = element_blank(),
                 plot.title = element_text(face = "bold")) +
           coord_flip() +
           facet_wrap(~pred_type, ncol = 2)
         
         list(obspreds_df   = obspreds_df,
              obspred_plot  = obspred_plot,
              residfit_plot = residfit_plot,
              qq_plot       = qq_plot,
              cooksd_plot   = cooksd_plot)
       })

# Helper: convert ggplot to base64 HTML image tag for inline rendering
plot_to_html <- function(plot_obj, width_in, height_in) {
  tmp <- tempfile(fileext = ".png")
  png(tmp, width = width_in * 150, height = height_in * 150, res = 150)
  print(plot_obj)
  dev.off()
  img_txt <- base64enc::dataURI(file = tmp, mime = "image/png")
  unlink(tmp)
  tags$img(src = img_txt,
           style = "width: 100%; max-width: 650px; display: block; margin-bottom: 20px;")
}

for (resp_var in names(diagnostics_plot_ls)) {
  cat("\n### ", resp_var, "\n\n", sep = "")
  current_dat <- diagnostics_plot_ls[[resp_var]]
  tab_content <- tagList(
    plot_to_html(current_dat$obspred_plot,  width_in = 6.5, height_in = 4),
    plot_to_html(current_dat$residfit_plot, width_in = 6.5, height_in = 4),
    plot_to_html(current_dat$qq_plot,       width_in = 6.5, height_in = 4),
    plot_to_html(current_dat$cooksd_plot,   width_in = 6.5, height_in = 8)
  )
  print(tab_content)
  cat("\n")
}

# ----7.0 Build prediction maps from MF's Mill River HR-NHD SSN models--------
# ---- 7.0 Pre-build diagnostics for ALL compounds (both raw + log1p) -----
# Section 6.4 diagnostics_plot_ls was filtered to final_key_compounds only.
# For per-compound comparison plots we need both scales, so rebuild here
# without the filter. Stored separately to not overwrite Section 6 output.
library(sf)
library(tidyverse)
library(SSN2)
library(patchwork)
library(wesanderson)
library(scales)

out_dir <- "PFAS_maps_HR"
if (!dir.exists(out_dir)) dir.create(out_dir)



diagnostics_plot_all_ls <-
  map2(.x = fit_model_ls,                  # ALL compounds, no filter
       .y = best_PFAS_models[names(fit_model_ls)],
       .f = function(x, y) {
         
         all_model_obspreds_df <-
           map2(.x = x, .y = names(x),
                .f = function(x, y) {
                  SSN2::augment(x = x) |>
                    dplyr::select(pid, observations = 1, starts_with(".")) %>%
                    dplyr::mutate(model_name = y) %>%
                    st_drop_geometry() %>% data.frame()
                }) |>
           bind_rows() |>
           rename(predictions = .fitted) |>
           mutate(pred_type = "Individual model predictions",
                  pid = as.numeric(pid))
         
         avg_model_obspreds_df <-
           left_join(x = all_model_obspreds_df,
                     y = y |> mutate(model_name = paste0(response_var, "_", rank)),
                     by = "model_name") %>%
           dplyr::mutate(weighted_pred = predictions * AIC_weight_i) %>%
           reframe(.by    = c(pid),
                   observations = mean(observations, na.rm = TRUE),
                   predictions  = sum(weighted_pred),
                   model_name   = "Model_Average") |>
           mutate(pred_type = "AIC-weighted predictions",
                  pid = as.numeric(pid))
         
         bind_rows(all_model_obspreds_df, avg_model_obspreds_df) |>
           left_join(select(obs_covariates, pid, site_id), by = "pid")
         # returns obspreds_df directly — simpler than full list
         # since maps section only needs obspreds_df
       })


# ---- 7.0. CHECK YOUR JOIN KEY BEFORE RUNNING FURTHER -----------------
# Your original script joined predictions to edges by `comid`. This
# Rmd never explicitly attaches reach_id (or any id) onto
# ssn_obj_obspreds_std$edges itself -- only onto obs/preds. Confirm:
print(names(ssn_obj_obspreds_std$edges))
# If you see "reach_id", set join_key <- "reach_id" below.
# If you only see "rid", set join_key <- "rid" and confirm
# preds_covariates_std$rid lines up with edges$rid (it should, since
# both come from the same lsn build).
join_key <- "reach_id"   # <-- change to "rid" if reach_id isn't present

# ---- 7.1. Compute AIC-weighted average predictions per compound -----
# For each compound: predict at each "preds" point from every
# retained model in its suite, weight each model's fitted value by
# its AIC_weight_i (already summing to 1 within a compound's set),
# and sum to get the model-averaged prediction.

get_weighted_preds <- function(compound) {
  mods    <- fit_model_ls[[compound]]
  weights <- best_PFAS_models[[compound]]$AIC_weight_i  # same row order as mods
  
  weighted_long <- map2_dfr(mods, weights, function(m, w) {
    SSN2::augment(m, newdata = "preds", pred.type = "preds") |>
      sf::st_drop_geometry() |>
      dplyr::transmute(pid = as.numeric(pid), weighted_pred = .fitted * w)
  })
  
  weighted_long |>
    dplyr::group_by(pid) |>
    dplyr::summarise(pred_val = sum(weighted_pred, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(compound = compound)
}

key_compounds <- names(fit_model_ls)  # e.g. PFAS40, PFSA, PFCA, PFOS_Results, PFOA_Results

weighted_preds_all <- purrr::map_dfr(key_compounds, get_weighted_preds)

# attach the join key (reach_id or rid) from the prediction points
pred_id_lookup <- preds_covariates_std |>
  sf::st_drop_geometry() |>
  dplyr::select(pid, dplyr::all_of(join_key)) |>
  dplyr::mutate(pid = as.numeric(pid))

weighted_preds_all <- weighted_preds_all |>
  dplyr::left_join(pred_id_lookup, by = "pid")

# wide format: one column per compound, one row per join_key
weighted_preds_wide <- weighted_preds_all |>
  dplyr::select(-pid) |>
  tidyr::pivot_wider(names_from = compound, values_from = pred_val,
                     names_glue = "{compound}_pred") |>
  dplyr::distinct(dplyr::across(dplyr::all_of(join_key)), .keep_all = TRUE)

# ---- 7.2. Join weighted predictions onto edges -----------------------
edges_with_preds <- ssn_obj_obspreds_std$edges |>
  dplyr::left_join(weighted_preds_wide, by = join_key)

pred_cols_present <- grep("_pred$", names(edges_with_preds), value = TRUE)
cat("Pred columns joined onto edges:", length(pred_cols_present), "\n")

# ---- 7.21 Pre-compute shared scale_max per base_compound (always in ng/L) ----
# Takes the maximum 0.98 quantile across both raw and log1p versions so
# RAW_PFOA and LOG_PFOA maps are on an identical color scale.

shared_scale_max <- list()

for (compound in names(fit_model_ls)) {
  is_log1p      <- startsWith(compound, "log1p_")
  base_compound <- ifelse(is_log1p, sub("^log1p_", "", compound), compound)
  pred_col      <- paste0(compound, "_pred")
  
  if (!pred_col %in% names(edges_with_preds)) next
  
  preds <- edges_with_preds[[pred_col]]
  if (is_log1p) preds <- expm1(preds)   # back-transform to ng/L
  
  vals <- c(preds, obs_sf[[base_compound]])
  qval <- quantile(vals, 0.99, na.rm = TRUE)
  
  # Keep the larger of the two (raw vs log1p) so both use the same ceiling
  if (is.null(shared_scale_max[[base_compound]])) {
    shared_scale_max[[base_compound]] <- qval
  } else {
    shared_scale_max[[base_compound]] <- max(shared_scale_max[[base_compound]], qval)
  }
}

cat("Shared scale_max per compound:\n")
print(unlist(shared_scale_max))

# ---- 7.3. Per-compound maps -------------------------------------------
pal <- colorRampPalette(c("#012A4A", "#014F86",
                          wes_palette("Zissou1", type = "continuous"),
                          "#8B0000"))(256)

obs_sf <- obs_covariates_s  # has geometry + raw PFAS columns + site_id

for (compound in key_compounds) {
  
  pred_col <- paste0(compound, "_pred")
  
  # ---- Scale detection and display setup -------------------------------
  is_log1p      <- startsWith(compound, "log1p_")
  base_compound <- ifelse(is_log1p, sub("^log1p_", "", compound), compound)
  obs_col       <- base_compound
  
  #Prefix for filename and plot title
  prefix        <- ifelse(is_log1p, "LOG", "RAW")
  
  if (!pred_col %in% names(edges_with_preds)) {
    message("Skipping ", compound, " - no prediction column found")
    next
  }
  
  # Back-transform log1p predictions to ng/L for display
  plot_edges <- edges_with_preds
  if (is_log1p) {
    plot_edges[[pred_col]] <- expm1(plot_edges[[pred_col]])
  }
  
  vals      <- c(plot_edges[[pred_col]], obs_sf[[obs_col]])
  scale_min <- 0
  scale_max <- shared_scale_max[[base_compound]]
  
  best_model_info <- best_PFAS_models[[compound]] |> slice(1)
  pseudo_r2  <- best_model_info$pseudo.r.squared
  predictors <- best_model_info$predictor_var
  
  scale_label <- ifelse(is_log1p, "ng/L\n(back-transformed)", "ng/L")
  
  subtitle <- paste0(
    "[", prefix, " model]  ",                          # scale method in subtitle
    "R\u00b2 = ", round(pseudo_r2, 3),
    "\nPredictors: ", predictors
  )
  
  p <- ggplot() +
    geom_sf(data = plot_edges,
            aes(color = .data[[pred_col]]), linewidth = 1.5) +
    geom_sf(data = obs_sf,
            aes(color = .data[[obs_col]]),
            shape = 21, fill = "white", size = 2, stroke = 2.5,
            show.legend = FALSE) +
    scale_color_gradientn(
      colors   = pal,
      limits   = c(scale_min, scale_max),
      oob      = scales::squish,
      na.value = "grey80",
      name     = scale_label) +
    guides(color = guide_colorbar(
      barwidth       = grid::unit(0.6, "cm"),
      barheight      = grid::unit(6,   "cm"),
      label.theme    = element_text(size = 12),
      title.theme    = element_text(size = 13, face = "bold"),
      title.position = "top")) +
    coord_sf(datum = sf::st_crs(4326)) +
    labs(title    = paste0("[", prefix, "] Predicted ", base_compound),  # prefix in title
         subtitle = subtitle) +
    theme_classic() +
    theme(
      plot.title      = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle   = element_text(size = 10, hjust = 0.5, face = "italic",
                                     margin = margin(b = 10)),
      legend.title    = element_text(size = 13, face = "bold"),
      legend.text     = element_text(size = 12),
      legend.position = "right",
      axis.text.x     = element_text(size = 9, angle = 45, hjust = 1),
      axis.text.y     = element_text(size = 9)
    )
  
  # RAW_ or LOG_ prefix in filename
  full_path <- file.path(out_dir,
                         paste0(prefix, "_", base_compound, "_HR_map.png"))
  
  ggsave(filename = full_path, plot = p, width = 8, height = 6, dpi = 300)
  
  if (file.exists(full_path)) {
    message("Saved: ", basename(full_path))
  } else {
    warning("File NOT found after save: ", full_path)
  }
}

# ---- Obs vs Predicted scatter plots with R2 and AIC ----

for (compound in key_compounds) {
  
  is_log1p      <- startsWith(compound, "log1p_")
  base_compound <- ifelse(is_log1p, sub("^log1p_", "", compound), compound)
  prefix        <- ifelse(is_log1p, "LOG", "RAW")
  
  obspreds_raw <- diagnostics_plot_all_ls[[compound]]
  if (is.null(obspreds_raw)) {
    warning("No diagnostics found for ", compound, " — skipping.")
    next
  }
  
  # ---- Pull AIC for BOTH scales of this compound ----------------------
  # These are NOT comparable to each other (different likelihoods/scales)
  # but shown side by side for reference
  aic_raw <- tryCatch(
    best_PFAS_models[[base_compound]] |>
      slice(1) |> pull(AIC),
    error = \(e) NA_real_
  )
  
  aic_log <- tryCatch(
    best_PFAS_models[[paste0("log1p_", base_compound)]] |>
      slice(1) |> pull(AIC),
    error = \(e) NA_real_
  )
  
  # ---- Pull AIC-weighted LOOCV RMSPE for current model from 6.1 -------
  # rmspe_raw is in ng/L; rmspe_log1p is in log(ng/L+1) — label accordingly
  scale_row <- scale_comparison_df |>
    filter(compound == base_compound)
  
  rmspe_display <- if (is_log1p) {
    scale_row$rmspe_log1p_ng_per_L   # in-sample RMSE, back-transformed
  } else {
    scale_row$rmspe_raw               # LOOCV RMSPE, already in ng/L
  }
  
  rmspe_label <- if (is_log1p) {
    "RMSE = "          # in-sample — honest label
  } else {
    "LOOCV RMSPE = "   # true out-of-sample
  }
  
  # pseudo R2 for current model
  pseudo_r2 <- best_PFAS_models[[compound]] |> slice(1) |> pull(pseudo.r.squared)
  
  # ---- Annotation label -----------------------------------------------
  annot_label <- paste0(
    "R\u00b2 = ",           round(pseudo_r2,    3),                   "\n",
    "AIC\u2081\u209a = ",   round(aic_log,      1), "  (log\u2081p)", "\n",
    "AIC\u1d63\u2090\u1d67 = ", round(aic_raw,  1), "  (raw)",        "\n",
    rmspe_label,             round(rmspe_display, 2), " ng/L"
  )
  
  model_avg_df <- obspreds_raw |>
    filter(pred_type == "AIC-weighted predictions")
  
  if (is_log1p) {
    model_avg_df <- model_avg_df |>
      mutate(observations = expm1(observations),
             predictions  = expm1(predictions))
  }
  
  p <- ggplot(model_avg_df, aes(x = predictions, y = observations)) +
    geom_point(size = 3, alpha = 0.6, color = "steelblue") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                color = "red", linewidth = 1) +
    annotate("text",
             x = Inf, y = -Inf,
             label  = annot_label,
             hjust  = 1.1, vjust = -0.3,
             size   = 3.5, fontface = "bold",
             color  = "black",
             lineheight = 1.4) +
    labs(title    = paste0("[", prefix, "] Observed vs Predicted: ", base_compound),
         subtitle = ifelse(is_log1p,
                           "log\u2081p model \u2014 axes back-transformed to ng/L",
                           "Raw (untransformed) model"),
         x = "Predicted (ng/L)",
         y = "Observed (ng/L)") +
    theme_minimal() +
    theme(
      plot.title    = element_text(size = 14, face = "bold",  hjust = 0.5),
      plot.subtitle = element_text(size = 10, face = "italic", hjust = 0.5),
      axis.title    = element_text(size = 12, face = "bold"),
      axis.text     = element_text(size = 11)
    ) +
    coord_fixed(ratio = 1)
  
  full_path <- file.path(out_dir,
                         paste0(prefix, "_", base_compound, "_obspred_scatter.png"))
  ggsave(filename = full_path, plot = p, width = 6, height = 6, dpi = 300)
  
  if (file.exists(full_path)) {
    message("Saved: ", basename(full_path))
  } else {
    warning("File NOT found: ", full_path)
  }
}

# ---- Assess leverage and influence of extreme sites ----

for (compound in key_compounds) {
  
  obspreds_df <- diagnostics_plot_ls[[compound]]$obspreds_df |>
    filter(pred_type == "AIC-weighted predictions")
  
  # Calculate Cook's Distance threshold
  n <- nrow(obspreds_df)
  cooksd_threshold <- 4 / n
  
  # Identify high-leverage points (.hat > 0.3)
  high_leverage <- obspreds_df |>
    filter(.hat > 0.3) |>
    arrange(desc(.hat)) |>
    mutate(compound = compound)
  
  if (nrow(high_leverage) > 0) {
    cat("\n", compound, " - High leverage points (.hat > 0.3):\n", sep = "")
    print(high_leverage |> select(pid, observations, predictions, .hat, .cooksd))
  }
}

# ---- Refit PFOA without the extreme site to compare R2 ----

# Get the best model formula for PFOA
pfoa_formula <- as.formula(paste(
  best_PFAS_models[["PFOA_Results"]]$response_var[1], "~",
  best_PFAS_models[["PFOA_Results"]]$predictor_var[1]
))

# Fit with all data
fit_pfoa_all <- fit_model_ls[["PFOA_Results"]][[1]]
r2_all <- glance(fit_pfoa_all)$pseudo.r.squared
aic_all <- glance(fit_pfoa_all)$AIC

# Remove pid 8 from obs data and rebuild SSN object
obs_reduced <- obs_covariates_s |> filter(pid != 8)
ssn_reduced <- ssn_put_data(obs_reduced, ssn_obj_obspreds_std, 
                            name = "obs", resize_data = TRUE)

# CRITICAL: Recreate distance matrices for the reduced dataset
ssn_create_distmat(ssn_reduced, overwrite = TRUE)

# Refit without site 8
fit_pfoa_no_site8 <- ssn_lm(pfoa_formula, 
                            ssn.object = ssn_reduced,
                            tailup_type = "exponential",
                            additive = "afvArea",
                            estmethod = "ml")
r2_no_site8 <- glance(fit_pfoa_no_site8)$pseudo.r.squared
aic_no_site8 <- glance(fit_pfoa_no_site8)$AIC

# Compare
comparison <- data.frame(
  scenario = c("All sites (including pid=8 @ 217 ng/L)", "Excluding pid=8"),
  R2 = c(r2_all, r2_no_site8),
  AIC = c(aic_all, aic_no_site8),
  R2_drop = c(0, r2_all - r2_no_site8)
)

print(comparison)

# ---- 4. Multipanel comparison figure ---------------------------------
wanted <- c("PFAS40", "PFSA", "PFOS_Results", "PFCA", "PFOA_Results")
wanted <- intersect(wanted, key_compounds)  # only ones actually fit

title_map <- c(
  PFAS40 = "PFAS40", PFCA = "PFCA", PFSA = "PFSA",
  PFOS_Results = "PFOS", PFOA_Results = "PFOA"
)

# shared scale across the multipanel so compounds are comparable
shared_vals <- purrr::map_dbl(wanted, function(cp) {
  stats::quantile(c(edges_with_preds[[paste0(cp, "_pred")]], obs_sf[[cp]]),
                  0.98, na.rm = TRUE)
})
scale_max_multiplot <- max(shared_vals, na.rm = TRUE)

plot_list <- list()

for (compound in wanted) {
  pred_col    <- paste0(compound, "_pred")
  obs_col     <- compound
  clean_title <- ifelse(compound %in% names(title_map), title_map[[compound]], compound)
  
  # Extract model info for best model of this compound
  best_model_info <- best_PFAS_models[[compound]] |> 
    slice(1)  # first row is best model (lowest AIC)
  
  pseudo_r2 <- best_model_info$pseudo.r.squared
  predictors <- best_model_info$predictor_var
  
  # Create subtitle with R² and predictors on separate lines
  subtitle <- paste0(
    "R² = ", round(pseudo_r2, 3), 
    "\nPredictors: ", predictors
  )
  
  plot_list[[compound]] <- ggplot() +
    geom_sf(data = edges_with_preds,
            aes(color = .data[[pred_col]]), linewidth = 1.5) +
    geom_sf(data = obs_sf,
            aes(color = .data[[obs_col]]),
            shape = 21, fill = "white", size = 2, stroke = 2.5,
            show.legend = FALSE) +
    scale_color_gradientn(
      colors = pal,
      limits = c(0, scale_max_multiplot),
      oob = scales::squish,
      na.value = "grey80",
      name = "PFAS (ng/L)") +
    guides(color = guide_colorbar(
      barwidth  = grid::unit(0.4, "cm"),
      barheight = grid::unit(4, "cm"))) +
    labs(title = clean_title,
         subtitle = subtitle) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 6, hjust = 0.5, face = "italic"),
      axis.title = element_blank(),
      axis.text  = element_text(size = 7),
      legend.position = "right")
}

final_plot <- plot_list[[1]]
for (i in 2:length(plot_list)) final_plot <- final_plot + plot_list[[i]]
final_plot <- final_plot +
  plot_layout(ncol = min(3, length(plot_list)), guides = "collect") +
  plot_annotation(
    title = "Predicted PFAS concentrations — Mill River Watershed (HR-NHD, AIC-weighted)",
    caption = "Circles = observed values | Stream color = AIC-weighted predicted concentration",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))

final_plot

ggsave(filename = file.path(out_dir, "HR_multipanel_SSN_map.png"),
       plot = final_plot, width = 14, height = 8, dpi = 400)

# ---- 5. StreamCat-style covariate maps across the watershed -----------
# Uses reach_data (basin, RCA, NPDES, AgTile, lentic) joined onto edges.
# No catchment polygon in this Rmd's environment -- edges (lines) are
# colored directly by covariate value. Add a catchment boundary layer
# below if you have one from the LSN build step.

covariate_edges <- ssn_obj_obspreds_std$edges |>
  dplyr::left_join(reach_data, by = "reach_id")

supp_predictors <- c(
  "basin_area_km2", "pct_nlcd_imp_bas", "pct_nlcd_imp_rca",
  "npdes_density_km2_bas", "npdes_density_km2_rca",
  "pct_agtile_bas", "pct_agtile_rca", "is_lentic"
)
supp_predictors <- intersect(supp_predictors, names(covariate_edges))

supp_labels <- c(
  basin_area_km2         = "Drainage area (km2)",
  pct_nlcd_imp_bas       = "% impervious (basin)",
  pct_nlcd_imp_rca       = "% impervious (RCA)",
  npdes_density_km2_bas  = "NPDES density (basin)",
  npdes_density_km2_rca  = "NPDES density (RCA)",
  is_lentic              = "Lentic vs lotic"
)

for (pred in supp_predictors) {
  
  p <- ggplot() +
    geom_sf(data = covariate_edges,
            aes(color = .data[[pred]]), linewidth = 1) +
    # geom_sf(data = catchment, fill = NA, color = "black", linewidth = 0.8) +
    scale_color_viridis_c(
      name = supp_labels[[pred]],
      option = "magma",
      na.value = "gray90") +
    labs(title = supp_labels[[pred]]) +
    theme_void() +
    theme(
      plot.title   = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text  = element_text(size = 9))
  
  ggsave(filename = file.path(out_dir, paste0(pred, "_covariate_map.png")),
         plot = p, width = 8, height = 6, dpi = 300)
  
  message("Saved covariate map: ", pred)
}

cat("\nAll maps saved to:", normalizePath(out_dir), "\n")

#-------------------------------------------------------------------
# Sensitivity analysis: keep vs remove the high-leverage site (ssnid == 6)
# instead of always dropping it via remove_outliers <- TRUE.
#
# RUN THIS AFTER your Rmd has executed all the way through section 5.3
# (so best_PFAS_models, continuous, continuous_vars, reach_data, and
# ssn_obj already exist in your session -- run it as its own chunk at
# the end, don't insert it earlier in the pipeline).
#-------------------------------------------------------------------

library(SSN2)
library(tidyverse)
library(broom)

# ---- 1. Re-fetch genuinely untouched obs data -----------------------
# obs_covariates in your session already has the NAs applied (since the
# main pipeline ran with remove_outliers <- TRUE), so pull a fresh copy
# straight from the ssn object rather than reusing obs_covariates.
obs_covariates_raw <- ssn_get_data(ssn_obj, name = "obs") |>
  dplyr::select(rid, pid, ratio, snapdist, upDist, afvArea, locID, netID,
                site_id, reach_id, ssnid,
                dplyr::all_of(pfas_compounds), netgeom, geometry) |>
  dplyr::left_join(reach_data, by = "reach_id")

outlier_site_id    <- 6
affected_compounds <- c("PFAS40", "PFCA", "PFBA_Results", "PFBS_Results",
                        "PFHpA_Results", "PFHxA_Results", "PFOA_Results",
                        "PFPeA_Results")

make_obs_scenario <- function(remove_outlier) {
  obs <- obs_covariates_raw
  if (remove_outlier) {
    obs[obs$ssnid == outlier_site_id, affected_compounds] <- NA
  }
  obs
}

obs_scenarios_raw <- list(
  all_data        = make_obs_scenario(FALSE),
  outlier_removed = make_obs_scenario(TRUE)
)

# ---- 2. Standardize + build an ssn object for each scenario ---------
# Predictor covariates (continuous/continuous_vars) don't change between
# scenarios -- only the PFAS response columns get NA'd -- so this reuses
# the same standardization logic, just fed different obs data.

build_ssn_scenario <- function(obs_variant) {
  cont_s <- continuous |>
    dplyr::mutate(dplyr::across(dplyr::all_of(continuous_vars), stand)) |>
    dplyr::rename_with(.fn = ~ paste0(.x, "_s"), .cols = dplyr::everything())
  
  obs_s <- data.frame(obs_variant, cont_s) |> sf::st_as_sf()
  
  # Only obs data is needed here -- ssn_lm() fits from obs, and Cook's
  # distance for site 6 comes from augment() on the fitted obs, so preds
  # (built later in section 4.2) isn't required for this comparison.
  ssn_put_data(obs_s, ssn_obj, name = "obs", resize_data = FALSE)
}

ssn_scenarios <- purrr::map(obs_scenarios_raw, build_ssn_scenario)

# ---- 3. Refit only the AFFECTED compounds' already-selected best model ----
# Reuses the rank-1 formula from best_PFAS_models -- this checks whether
# the model you already picked is sensitive to this one site, without
# re-running the full AIC search twice.

affected_key_compounds <- intersect(affected_compounds, names(best_PFAS_models))

compare_scenarios <- function(compound) {
  form <- as.formula(paste(
    best_PFAS_models[[compound]]$response_var[1], "~",
    best_PFAS_models[[compound]]$predictor_var[1]
  ))
  
  fits <- purrr::map(ssn_scenarios, function(ssn_obj_scenario) {
    ssn_lm(form, ssn.object = ssn_obj_scenario,
           tailup_type = "exponential", additive = "afvArea", estmethod = "ml")
  })
  
  coefs_df  <- purrr::map_dfr(fits, ~broom::tidy(.x, conf.int = TRUE), .id = "scenario")
  glance_df <- purrr::map_dfr(fits, ~broom::glance(.x), .id = "scenario")
  
  # Cook's distance for THIS site specifically, from the all-data fit
  site6_diag <- SSN2::augment(fits$all_data) |>
    sf::st_drop_geometry() |>
    dplyr::mutate(pid = as.character(pid)) |>
    dplyr::left_join(
      dplyr::select(sf::st_drop_geometry(obs_covariates_raw), pid, ssnid) |>
        dplyr::mutate(pid = as.character(pid)),
      by = "pid") |>
    dplyr::filter(ssnid == outlier_site_id)
  
  list(compound = compound, fits = fits, coefs = coefs_df,
       glance = glance_df, site6_diag = site6_diag)
}

comparison_results <- purrr::map(affected_key_compounds, compare_scenarios)
names(comparison_results) <- affected_key_compounds

# ---- 4. Summary table: AIC / pseudo-R2 shift when the site is removed ----
comparison_summary <- purrr::map_dfr(comparison_results, function(x) {
  x$glance |>
    dplyr::mutate(compound = x$compound) |>
    dplyr::select(compound, scenario, AIC, AICc, pseudo.r.squared)
})

print(comparison_summary)

# ---- 5. Coefficient stability check ---------------------------------------
# If a coefficient flips sign or its magnitude changes drastically between
# scenarios, that's real evidence the site is driving the relationship
# rather than the model just underpredicting one point.
coef_comparison <- purrr::map_dfr(comparison_results, function(x) {
  x$coefs |>
    dplyr::mutate(compound = x$compound) |>
    dplyr::select(compound, scenario, term, estimate, std.error, p.value)
})

coef_comparison |>
  tidyr::pivot_wider(names_from = scenario, values_from = c(estimate, std.error, p.value)) |>
  dplyr::mutate(pct_shift = 100 * (estimate_all_data - estimate_outlier_removed) /
                  abs(estimate_outlier_removed)) |>
  print(n = Inf)

# ---- 6. Objective threshold: is site 6 actually influential? --------------
# Standard rule of thumb: Cook's D > 4/n flags a high-leverage point.
purrr::walk(comparison_results, function(x) {
  n <- nrow(sf::st_drop_geometry(ssn_scenarios$all_data$obs))
  threshold <- 4 / n
  cat(sprintf("%-15s Cook's D at site 6: %6.4f | threshold: %6.4f | flagged: %s\n",
              x$compound, x$site6_diag$.cooksd, threshold,
              x$site6_diag$.cooksd > threshold))
})

# ---- Interpretation guide ---------------------------------------------------
# - If Cook's D for site 6 is well above 4/n AND coefficients/AIC shift a lot
#   between scenarios -> exclusion is defensible; document why (e.g. known
#   local source) rather than a blanket TRUE/FALSE flag.
# - If Cook's D is near/below threshold OR coefficients barely move -> keep
#   the site in; the model just isn't capturing that local signal, which is
#   worth noting as a limitation rather than grounds for removal.
# - Either way, report BOTH scenarios in your writeup so the choice is
#   transparent and reviewers can see what changes.

#-------------------------------------------------------------------
# Extend the Cook's D sensitivity check to ALL compounds site 6 was
# originally NA'd for (not just PFAS40/PFCA/PFOA_Results), and build
# a per-compound exclusions table -- not a single site-level flag --
# since influence is model/response-specific.
#
# RUN THIS AFTER outlier_sensitivity_comparison.R has already run
# (reuses ssn_scenarios, obs_covariates_raw, outlier_site_id from that
# script).
#-------------------------------------------------------------------

library(SSN2)
library(tidyverse)
library(broom)

# The five compounds never actually tested -- they don't have a
# pre-selected "best" formula from best_PFAS_models since only
# key_compounds went through the AIC search. Use the same predictor
# set as the PFCA family model (their parent family) as a reasonable
# default -- swap in a formula of your choice if you'd rather run each
# through its own AIC search.
untested_compounds <- c("PFBA_Results", "PFBS_Results", "PFHpA_Results",
                        "PFHxA_Results", "PFPeA_Results")

# Use the RIGHT family's formula per compound -- PFCA family (carboxylates)
# vs PFSA family (sulfonates) can have different covariate relationships.
pfca_family <- c("PFBA_Results", "PFHpA_Results", "PFHxA_Results", "PFPeA_Results")
pfsa_family <- c("PFBS_Results")

fallback_formula_lookup <- c(
  setNames(rep(best_PFAS_models[["PFCA"]]$predictor_var[1], length(pfca_family)), pfca_family),
  setNames(rep(best_PFAS_models[["PFSA"]]$predictor_var[1], length(pfsa_family)), pfsa_family)
)

test_compound_cooksd <- function(compound, formula_rhs) {
  form <- as.formula(paste(compound, "~", formula_rhs))
  
  fit_all <- ssn_lm(form, ssn.object = ssn_scenarios$all_data,
                    tailup_type = "exponential", additive = "afvArea",
                    estmethod = "ml")
  
  n <- nrow(sf::st_drop_geometry(ssn_scenarios$all_data$obs))
  threshold <- 4 / n
  
  site6_diag <- SSN2::augment(fit_all) |>
    sf::st_drop_geometry() |>
    dplyr::mutate(pid = as.character(pid)) |>
    dplyr::left_join(
      dplyr::select(sf::st_drop_geometry(obs_covariates_raw), pid, ssnid) |>
        dplyr::mutate(pid = as.character(pid)),
      by = "pid") |>
    dplyr::filter(ssnid == outlier_site_id)
  
  tibble(compound = compound, cooksd = site6_diag$.cooksd,
         threshold = threshold, flagged = site6_diag$.cooksd > threshold)
}

untested_results <- purrr::map_dfr(untested_compounds,
                                   test_compound_cooksd,
                                   formula_rhs = fallback_formula_rhs)

# combine with the three already tested in outlier_sensitivity_comparison.R
tested_results <- purrr::map_dfr(comparison_results, function(x) {
  n <- nrow(sf::st_drop_geometry(ssn_scenarios$all_data$obs))
  tibble(compound = x$compound, cooksd = x$site6_diag$.cooksd,
         threshold = 4 / n, flagged = x$site6_diag$.cooksd > (4 / n))
})

all_cooksd_results <- dplyr::bind_rows(tested_results, untested_results)
print(all_cooksd_results)

# ---- Build the exclusions table -------------------------------------
# One row per (site, compound) actually flagged -- not a blanket
# site-level exclusion. This is what drives NA-ing downstream.
site_exclusions <- all_cooksd_results |>
  dplyr::filter(flagged) |>
  dplyr::mutate(
    site_id  = outlier_site_id,
    reason   = "Known point-source contamination; Cook's D exceeds 4/n threshold",
    excluded_on = Sys.Date()
  ) |>
  dplyr::select(site_id, compound, cooksd, threshold, reason, excluded_on)

print(site_exclusions)

# save it so it's a persistent, documented record rather than a
# hardcoded ssnid == 6 check buried in a script
write_csv(site_exclusions, "../data/ssn_site_exclusions.csv")

# ---- How to apply it when building modeling input --------------------
# Generic version of the old remove_outliers block: NA out only the
# specific (site, compound) pairs recorded in the exclusions table,
# instead of a fixed compound list.
apply_exclusions <- function(obs_data, exclusions) {
  for (i in seq_len(nrow(exclusions))) {
    row <- exclusions[i, ]
    obs_data[obs_data$ssnid == row$site_id, row$compound] <- NA
  }
  obs_data
}

# example usage in place of the old if(isTRUE(remove_outliers)) block:
# obs_covariates <- apply_exclusions(obs_covariates_raw, site_exclusions)
```