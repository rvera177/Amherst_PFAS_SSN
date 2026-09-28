
#   title: "Mill River PFAS: SSN object development"
# author: "Matthew Fuller and Raul Vera"
# date: "September 28, 2026"

# This document describes the process and code used to explore PFAS concentrations
# in the Mill River watershed near Amherst, MA using spatial stream network (SSN) models.  
# Specifically, we build the SSN object in this script and it takes less than
# 5 minutes to run all the code.  


# ---- 1.  Set Up Analysis ------------ 

library(viridis)
library(stringr)
library(knitr)
library(pander)
library(tidyverse)
library(sf)
library(SSNbler)
library(SSN2)
library(readr)
panderOptions('digits',8)

# Run from here, all paths are relative to data/
if (basename(getwd()) != "data") setwd("data")
source("../scripts/00 PFAS naming.R")  # PFAS compounds, functional groups and column-naming rules

# set tidy r code not to run more than 100 characters wide
opts_chunk$set(tidy.opts=list(width.cutoff=100),tidy=FALSE) 

knitr::opts_chunk$set(message = FALSE) # hide all warning messages in document

# ----- 2.  Read in data  ------------

# Network data  
covariate_data <- readRDS("../data/RVsummary_basin_covariates.rds") |>
  left_join(readRDS("../data/RVsummary_rca_covariates.rds"), by = "reach_id")  |>
  left_join(readRDS("../data/RVsummary_riparian_covariates.rds"),
            by = "reach_id")

flines  <- st_read(dsn = "../data/shp/MillRiver_NHDPlusHR_EPSG5070_trace.shp") |>
  left_join(covariate_data, by = "reach_id")

flowline <- st_transform(flines, crs =5070) |> st_zm() 

##  Observation sites  

#bringing in the results from my PFAS analysis
PFAS_Spatial_Results_All <- 
  read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/Spatial_PFASResults.csv")

# clean_names() is defined in 00_PFAS_naming.R (keeps the leading "X" so names like X4.2FTS are valid in model formulas)
colnames(PFAS_Spatial_Results_All) <- clean_names(colnames(PFAS_Spatial_Results_All))
#names(PFAS_Spatial_Results_All)

#Remove results below Method Detection Limit
# Concentrations below this limit are not reportable.
mdl_table <- 
  read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/EPA1633_MDL_UMass_Amherst_Engineering.csv")
names(mdl_table) <- c("compound", "mdl_ng_l", "mrl_ng_l")
mdl_table <- mdl_table |> mutate(compound_clean = clean_names(compound))

pfas_result_cols <- grep("_Results$", names(PFAS_Spatial_Results_All), value = TRUE)

mdl_lookup <- tibble(result_col = pfas_result_cols) |>
  mutate(compound_clean = sub("_Results$", "", result_col)) |>
  left_join(mdl_table |> select(compound_clean, mdl_ng_l), by = "compound_clean")

# Every result column must have an MDL, and the compounds in the data must match 00_config.R
stopifnot(!anyNA(mdl_lookup$mdl_ng_l),
          setequal(mdl_lookup$compound_clean, pfas_individual))

S_clean <- PFAS_Spatial_Results_All
for (i in seq_len(nrow(mdl_lookup))) {
  col   <- mdl_lookup$result_col[i]
  mdl_v <- mdl_lookup$mdl_ng_l[i]
  if (is.na(mdl_v)) next
  S_clean[[col]] <- ifelse(S_clean[[col]] < mdl_v, 0, S_clean[[col]])
}

# Functional head groups (pfas_groups) are defined in 00_PFAS_naming.R

# Create PFAS group sums columns for each sampling event,
# and removes the "_Results" suffix, and adds OBSPRED_ID

prepare_event <- function(df) {
  df |>
    # Add family group sum columns
    bind_cols(lapply(names(pfas_groups), function(fam_name) {
      cols          <- paste0(pfas_groups[[fam_name]], "_Results")
      cols_existing <- intersect(cols, names(df))
      tibble(!!fam_name := rowSums(df[cols_existing], na.rm = TRUE))
    })
    ) |>
    # Compute PFAS40 total across all individual compounds
    mutate(PFAS40 = rowSums(
      across(all_of(intersect(pfas_result_cols, names(df)))),
      na.rm = TRUE)
    ) |>
    # Strip _Results suffix from compound columns
    rename_with(~ gsub("_Results$", "", .x)) |>
    # Add log1p transformed columns for all numeric PFAS columns
    mutate(
      across(
        .cols = all_of(
          gsub("_Results$", "",
               c(pfas_result_cols,
                 paste0(names(pfas_groups), "_Results")))
        ),
        .fns  = list(log1p = ~ log1p(.x)),
        .names = "log1p_{.col}"
      )
    ) |>
    mutate(OBSPRED_ID = row_number())
}

sampling_events <- S_clean |>
  dplyr::distinct(Sampling_Event) |>
  dplyr::pull(Sampling_Event) |>
  sort()
cat("Sampling events found:", paste(sampling_events, collapse = ", "), "\n")

S_list <- sampling_events |>
  purrr::set_names() |>
  purrr::map(function(event) {
    S_clean |>
      dplyr::filter(Sampling_Event == event) |>
      prepare_event()
  })

#Make the spatial objects 
make_obs_sf <- function(S_event) {
  # build points from CSV coordinates
  obs_pts <- S_event |>
    sf::st_as_sf(
      coords = c("Long", "Lat"),
      crs    = 4326
    ) |>
    sf::st_transform(crs = 5070) |>
    sf::st_zm()
  
  nearest_idx      <- sf::st_nearest_feature(obs_pts, flowline) # snap each point to nearest flowline
  obs_pts$reach_id <- flowline$reach_id[nearest_idx] #assign the reach_id from the closest flowline
  obs_pts$snap_dist_m <- sf::st_distance( # compute snap distance for QC
    obs_pts, 
    flowline[nearest_idx, ],
    by_element = TRUE
  ) |> as.numeric()
  #add covariates
  obs_pts |>
    dplyr::left_join(covariate_data, by = "reach_id") |>
    dplyr::mutate(ssnid = row_number())
}

obs_sf_list <- S_list |> purrr::map(make_obs_sf)

## Predictions sites  

#These are placed at the center of each flowline, but will be represented by data sets
#aggregated/summarized for the downstream pour point of the reach/flowline. Since the NHDPlus 
#high-resolution network has relatively short reach segments (mean length of ##m) this shouldn't be 
#problematic or poorly representative of the reach.

pred <- flines |> 
  st_zm() |>                  # drop Z and M dimensions if they exist
  st_centroid() |>            # centroid works for LINESTRING objects
  st_cast("POINT") |>         # makes sure geometry is POINT
  mutate(ssnid = reach_id + 10000) |>
  select(ssnid, reach_id, NHDPlusIDt, 
         ends_with(c("_km2","_bas","_rca","_rip","_rip_bas")),
         lengthkm, areasqkm, totdasqkm, geometry) |>
  st_transform(crs = 5070) |> st_zm() 

#--  3.  Build SSN object  ---------

# Directory to hold the files during SSN object development.  
ssnobj_devdir <- "../ssnobj_dev/"
dir.create(ssnobj_devdir, showWarnings = FALSE, recursive = TRUE)

##  3.1. Develop the Landscape network (LSN) 
flowlines_2 <- lines_to_lsn(streams = flowline, 
                            lsn_path = ssnobj_devdir,
                            check_topology = TRUE,
                            snap_tolerance = 0,
                            topo_tolerance = 0,
                            remove_ZM = TRUE,
                            overwrite = TRUE)

##---3.2. Bring in sites ------- 

#Observation sites  
obs_lsn_list <- obs_sf_list |>
  purrr::imap(function(obs_event, event_name) {
    sites_to_lsn(
      sites          = obs_event,
      edges          = flowlines_2,
      lsn_path       = ssnobj_devdir,
      file_name      = paste0("obs_", event_name), 
      snap_tolerance = 100,
      save_local     = TRUE,
      overwrite      = TRUE)
  })

#Prediction sites  
preds <- sites_to_lsn(sites = pred,
                      edges = flowlines_2,
                      lsn_path = ssnobj_devdir,
                      file_name = "preds",
                      snap_tolerance = 250,
                      save_local = TRUE,
                      overwrite = TRUE)


##----  3.3.  Calculate upstream distances for edges and sites  -------
edges <- updist_edges(edges = flowlines_2,
                      save_local = TRUE,
                      lsn_path = ssnobj_devdir,
                      calc_length = TRUE)

all_sites_for_updist <- c(
  obs_lsn_list,
  list(preds = preds))

site.list <- updist_sites(
  sites       = all_sites_for_updist,
  edges       = edges,
  length_col  = "Length",
  save_local  = TRUE,
  lsn_path    = ssnobj_devdir)


##---  3.4. Calculate additive function values  ---------
edges <- afv_edges(edges = edges,
                   infl_col = "basin_area_km2",
                   segpi_col = "areaPI",
                   afv_col = "afvArea",
                   lsn_path = ssnobj_devdir)

site.list <- afv_sites(sites = site.list,
                       edges = edges,
                       afv_col = "afvArea",
                       save_local = TRUE,
                       lsn_path = ssnobj_devdir)


# Visual to see how far each site had to snap onto the flowline, 
# and if it snapped to the correct location.

# Samples collected at the edge of a pond/reservoir/lake snapped to the center, 
# resulting in snap distances over 50 meters, which is expected.

ggplot() +
  geom_sf(data = flowline, color = "steelblue", linewidth = 0.4) +
  geom_sf(data = obs_sf_list[["S3"]], #this is for any specific spatial event. Here, I map S3.
          aes(color = snap_dist_m), size = 3) +
  scale_color_viridis_c(name = "Snap\ndistance (m)",
                        option = "plasma") +
  labs(title    = "S3 observation sites snapped to nearest flowline",
       subtitle = "Color = distance from original GPS point to nearest reach") +
  theme_minimal()

#  --- 4. Assemble SSN object  ----
dir.create("../ssnobj", showWarnings = FALSE, recursive = TRUE)
ssn_list <- names(obs_lsn_list) |>
  purrr::set_names() |>
  purrr::map(function(event_name) {
    obs_key <- paste0("obs_", event_name)
    ssn_assemble(
      edges      = edges,
      lsn_path   = ssnobj_devdir,
      obs_sites  = site.list[[event_name]],
      preds_list = site.list["preds"],
      ssn_path   = paste0("../ssnobj/PFAS_NHDHR_", event_name, ".ssn"),
      import     = TRUE,
      check      = TRUE,
      afv_col    = "afvArea",
      overwrite  = TRUE,
      verbose    = TRUE)
  })

