
#  title: "Mill River PFAS: Covariate development (generalized pipeline)"
#author: "adapted from Matt Fuller's original script"
#date: "8/13/26"
# output: html_document

  
  # Purpose
  # This code generalizes Matt's basin/RCA covariate workflow so new predictor rasters
  # (and re-derived layers can be added by editing
  # ONE registry list below, instead of hand-writing a new code block per raster.
  #
  # Two computation paths are preserved on purpose:
  #   - BASIN path: basins are nested/overlapping -> loop over basin_stack layers,
  #     mask() each covariate raster to each layer. Required because overlapping
  #     polygons cannot be represented as a single categorical zone raster.
  #   - RCA path: RCAs are mutually exclusive -> single categorical zone raster,
  #     use zonal() once across all reaches. Much faster; don't loop here.
  
library(whitebox); wbt_init()
library(sf)
library(terra)
library(tidyverse)
library(tidyterra)
getwd()
list.files("../outputs/basins/", pattern = "basin_", full.names = TRUE)
setwd("C:/Users/Marston User/Documents/LWMR Isoscape/PFAS_MillRiver_SSN-main/PFAS_MillRiver_SSN-main/data")

# ---- 1. Read in basin stack, RCA zone raster, and pour points --------------

# Basin stack: rebuild from individual files (or swap to readRDS if you trust the cache)
reach_basins <-
  list.files(path = "../outputs/basins/", pattern = "basin_", full.names = TRUE) |>
  (\(x) x[basename(x) != "basin_.tif"])() |>   # keep consistent w/ hydrography exclusion
  lapply(terra::rast)
basin_stack <- rast(reach_basins)
template <- basin_stack[[1]]

# RCA zone raster: this is the reclassified output from Step 6b, not final_pour_point_RCAs.tif
reach_rca_r <- rast("../data/raster/reach_rca_multipoint.tif")
cell_area_m2 <- prod(res(reach_rca_r))

# Pour points — same source shapefile used to build the zone_lookup crosswalk
final_pour_points <-
  st_read("../data/shp/final_pour_points_manual_corrections_20260518.shp") |>
  st_transform(crs(basin_stack))


# ---- 1b. Validate reach_id correspondence across basin/RCA/pour points -----

# Pull reach_id out of basin_stack layer names and set them explicitly
basin_ids <- as.numeric(gsub("basin_", "", names(basin_stack)))
names(basin_stack) <- basin_ids   # now layer name == reach_id, not filename

rca_ids   <- sort(unique(values(reach_rca_r)))
rca_ids   <- rca_ids[!is.na(rca_ids) & rca_ids != 0]   # drop background/NA

pp_ids    <- sort(unique(final_pour_points$reach_id))

# Three-way comparison
#this is checking that the basin and RCA reach_id's correspond 
# to each other correctly. you want to see numeric(0) results for the list
# and n basin layers: 479 | n RCA zones: 479 | n pour points: 479 
list(
  basin_not_in_rca   = setdiff(basin_ids, rca_ids),
  rca_not_in_basin   = setdiff(rca_ids, basin_ids),
  basin_not_in_pp    = setdiff(basin_ids, pp_ids),
  pp_not_in_basin    = setdiff(pp_ids, basin_ids),
  rca_not_in_pp      = setdiff(rca_ids, pp_ids),
  pp_not_in_rca      = setdiff(pp_ids, rca_ids)
) |> print()

cat("n basin layers:", length(basin_ids),
    "| n RCA zones:", length(rca_ids),
    "| n pour points:", length(pp_ids), "\n")

# ---- 2. Alignment helper (crops in native CRS BEFORE reprojecting) ---------
#
# For CONUS-wide (or otherwise huge) source rasters, reprojecting the FULL
# extent first is very slow and unnecessary. Instead: reproject just the
# template's bounding box into the source raster's native CRS, crop to that
# , THEN reproject the small clipped piece. A buffer
# is added before the native-CRS crop because reprojecting a rectangle can
# distort it into a non-rectangular shape - the buffer guards against
# clipping off real edge pixels before the final resample.

get_native_crop_extent <- function(template, target_crs, buffer_m = 1000) {
  template_poly      <- as.polygons(ext(template), crs = crs(template))
  template_poly_proj <- project(template_poly, target_crs)
  ext(template_poly_proj) + buffer_m
}
align_to_template <- function(path, template, method = c("near", "bilinear"),
                              pre_clip_buffer_m = 1000, cache_dir = NULL) {
  method <- match.arg(method)
  
  cache_path <- NULL
  if (!is.null(cache_dir)) {
    cache_path <- file.path(
      cache_dir,
      paste0(tools::file_path_sans_ext(basename(path)), "_aligned.tif")
    )
    if (file.exists(cache_path)) {
      message("  Loading from cache: ", basename(cache_path))
      return(rast(cache_path))
    }
  }
  
  r <- rast(path)
  
  if (!same.crs(r, template)) {
    native_ext <- get_native_crop_extent(template, crs(r),
                                         buffer_m = pre_clip_buffer_m)
    r <- crop(r, native_ext, snap = "out")
  } else {
    r <- crop(r, ext(template) + pre_clip_buffer_m, snap = "out")
  }
  
  # ✅ Single call: reprojects CRS + snaps to template grid + resamples
  # Replaces the broken crop() |> extend() |> resample() chain
  # extend() was silently failing for coarse rasters (e.g. 1km BFI)
  # because it cannot add fractional cells to match a finer template
  r_aligned <- project(r, template, method = method)
  
  if (!is.null(cache_path)) {
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    writeRaster(r_aligned, cache_path, overwrite = TRUE)
    message("  Cached: ", basename(cache_path))
  }
  
  r_aligned
}


# ---- Crop WorldPop to study area (before running main covariate script) ----
library(terra)

pop_full <- rast("../data/raster/usa_pop_2025_CN_100m_R2025A_v1.tif")

# Amherst, MA is roughly 42-43°N, 72-73°W
# Add buffer for safety
pop_cropped <- crop(pop_full, ext(-73.5, -71.5, 41.5, 43.5))

writeRaster(pop_cropped, "../data/raster/pop_2025_cropped.tif", overwrite = TRUE)


# ---- 3. Covariate registry --------------------------------------------------
# This is the only thing you edit to add a new predictor.
#   type:
#     "percent_cover" - binary/categorical raster (e.g. AgTile, wetlands mask);
#                        reports % of cells that are non-NA / "present"
#     "continuous"     - continuous raster (e.g. NLCD imperviousness, slope,
#                        precip); reports the mean value
#     "point_density"  - point vector layer (e.g. NPDES, dams); reports count
#                        and density per km2
#   resample_method: "near" for categorical/binary, "bilinear" for continuous
#
# For point layers, `path` should point to the vector file, not a raster.
# For raster layers, `path` should point to the GeoTIFF.
# reclassify_class: for binary land cover masks, specify the NLCD class(es) to keep as 1
# clipping population density first. too high resolution

covariate_registry <- list(
  #2025 impervious data source. https://www.mrlc.gov/data/type/fractional-impervious-surface
  nlcd_imp2025 = list( 
    path = "../data/raster/Annual_NLCD_FractionalImperviousSurface_2025_CU_C1V2/Annual_NLCD_FctImp_2025_CU_C1V2.tif",
    type = "continuous",
    resample_method = "bilinear",
    out_name = "impervious_2025"),
  # 2025 NLCD Cultivated Crops (class 81): Areas used for production of annual crops 
  # Data source: https://www.usgs.gov/centers/eros/science/annual-nlcd-land-cover-classification 
  # (corn, soybeans) and perennial woody crops (orchards). Crop vegetation > 20%.
  nlcd_cultivated_crops_2025 = list(
    path = "../data/raster/Annual_NLCD_LndCov_2025_CU_C1V2/Annual_NLCD_LndCov_2025_CU_C1V2.tif",
    type = "percent_cover",
    resample_method = "near",
    out_name = "crops_2025",
    reclassify_class = 81),
  # 2025 NLCD Pasture/Hay (class 82): Grasses, legumes, or grass/legume mixtures
  # planted for livestock grazing or hay/seed production on perennial cycle. 
  # Pasture/Hay vegetation > 20%.
  nlcd_pasture_hay_2025 = list(
    path = "../data/raster/Annual_NLCD_LndCov_2025_CU_C1V2/Annual_NLCD_LndCov_2025_CU_C1V2.tif",
    type = "percent_cover",
    resample_method = "near",
    out_name = "pasture_hay_2025",
    reclassify_class = 82),
  # 2025 NLCD Planted/Cultivated (classes 81 + 82): Combined cultivated crops and 
  # pasture/hay. Both fall under Planted/Cultivated classification and are sources 
  # of non-point contamination from pesticides and biosolid applications.
  nlcd_planted_cultivated_2025 = list(
    path = "../data/raster/Annual_NLCD_LndCov_2025_CU_C1V2/Annual_NLCD_LndCov_2025_CU_C1V2.tif",
    type = "percent_cover",
    resample_method = "near",
    out_name = "crops_and_pasture_2025",
    reclassify_class = c(81, 82)),
  npdes = list(
    path = "../data/shp/USEPA_NPDES_pts.shp",
    type = "point_density",
    out_name = "npdes"),
  # Baseflow Index (BFI): 1km resolution, interpolated from USGS streamgages
  # Represents base flow component of streamflow from ground-water discharge
  # link to source: https://www.sciencebase.gov/catalog/item/631405c5d34e36012efa3192
  # actual raster data is in the pfi48grd.zip file
  bfi = list(
    path              = "../data/raster/BaseFlowIndex_USGS_48grd/bfi48grd/w001001.adf",
    type              = "continuous",
    resample_method   = "bilinear",
    out_name          = "bfi",
    pre_clip_buffer_m = 5000 #BFI is 1 km cells and want 5× cell size as buffer
  ),
  # POLARIS Soil Properties: https://www.isric.org/evaluate/regional-soil-property-datasets/polaris-30m-americas
  # Organic Matter (OM): Affects PFAS sorption and degradation
  organic_matter = list(
    path = "../data/raster/OrganicMatter_lat4243_lon-73-72.tif",
    type = "continuous",
    resample_method = "bilinear",
    out_name = "om"),
  # Clay content: Sorption potential, affects permeability
  clay = list(
    path = "../data/raster/Clay_lat4243_lon-73-72.tif",
    type = "continuous",
    resample_method = "bilinear",
    out_name = "clay_pct"),
  # Saturated hydraulic conductivity (Ksat): Water/contaminant transport rate
  ksat = list(
    path = "../data/raster/hydraulicConductivity_lat4243_lon-73-72.tif",
    type = "continuous",
    resample_method = "bilinear",
    out_name = "ksat"),
  # populationdensity of US at 30m res for 2025
  # source: WorldPop https://hub.worldpop.org/geodata/summary?id=75983
  pop2025 = list(
    path = "../data/raster/pop_2025_cropped.tif",
    type = "continuous",
    resample_method = "bilinear",
    out_name = "pop_density_2025")
)


# ---- 4. Pre-align all registered rasters / load all point layers -----------

aligned_rasters <- list()
point_layers    <- list()

cache_dir <- "../data/raster_cache/"
dir.create(cache_dir, showWarnings = FALSE)

for (nm in names(covariate_registry)) {
  cov <- covariate_registry[[nm]]
  
  if (cov$type == "point_density") {
    point_layers[[nm]] <- st_read(cov$path, quiet = TRUE) |>
      vect() |> project(crs(template))
    
  } else {
    # Use per-covariate buffer if specified, otherwise default 1000 m
    buf <- if (!is.null(cov$pre_clip_buffer_m)) cov$pre_clip_buffer_m else 1000
    
    r_aligned <- align_to_template(
      path              = cov$path,
      template          = template,
      method            = cov$resample_method,
      pre_clip_buffer_m = buf,          # ← now passed through correctly
      cache_dir         = cache_dir
    )
    
    if (!is.null(cov$reclassify_class)) {
      r_aligned <- ifel(r_aligned %in% cov$reclassify_class, 1, NA)
    }
    aligned_rasters[[nm]] <- r_aligned
  }
}

# Sanity check alignment (fix: don't unlist SpatRasters)
compareGeom(template, aligned_rasters[[1]])

# ---- 5. Basin-side summary functions ----------------------------------------

summarize_percent_cover_basin <- function(basin_mask, data_raster) {
  data_masked <- mask(data_raster, basin_mask)
  n_cells <- sum(!is.na(values(basin_mask)))
  100 * sum(!is.na(values(data_masked))) / n_cells
}

summarize_continuous_mean_basin <- function(basin_mask, data_raster) {
  data_masked <- mask(data_raster, basin_mask)
  vals <- values(data_masked)
  mean(vals[!is.na(vals)])
}

summarize_point_density_basin <- function(basin_mask, pts) {
  n_cells <- sum(!is.na(values(basin_mask)))
  basin_area_km2 <- (n_cells * cell_area_m2) / 1e6
  
  pts_subset <- pts[ext(basin_mask), ]
  if (nrow(pts_subset) == 0) {
    n_pts <- 0
  } else {
    vals <- terra::extract(basin_mask, pts_subset)
    n_pts <- sum(!is.na(vals[, 2]))
  }
  list(n_points = n_pts, density_km2 = n_pts / basin_area_km2)
}


# ---- 6. Basin-side loop (generic over the registry) -------------------------
# so this loop is checking if i already have any covariates developed at the basin scale
# first it loads the data and checks what basin covariates i have. any new ones that are 
# in the registry but not in the results_cache will be run. Covariates that have already been
# calculated won't be re-done. This is done because it takes minutes for a covariate to load
# and i am constantly adding/iterating new covariates to the list and don't want to wait each time


# Create a valid watershed boundary (union of all basins with data)
watershed_mask <- !is.na(basin_stack)
watershed_mask <- any(watershed_mask)  # TRUE where any basin has data

results_cache <- "../data/RVsummary_basin_covariates.csv"

# Load existing results if they exist
if (file.exists(results_cache)) {
  results_old <- read_csv(results_cache)
  existing_cols <- names(results_old)
  cat("Loaded", nrow(results_old), "basins with", ncol(results_old), "columns.\n")
} else {
  results_old <- NULL
  existing_cols <- c("reach_id", "basin_area_km2")  # always have these
}

# Identify NEW covariates (not in existing results)
covs_to_calc <- names(covariate_registry)
new_covs <- c()

for (nm in covs_to_calc) {
  cov <- covariate_registry[[nm]]
  if (cov$type == "percent_cover") {
    col_name <- paste0("pct_", cov$out_name, "_bas")
  } else if (cov$type == "continuous") {
    col_name <- paste0(cov$out_name, "_bas")
  } else if (cov$type == "point_density") {
    col_name <- paste0(cov$out_name, "_count_bas")
  }
  
  if (!col_name %in% existing_cols) {
    new_covs <- c(new_covs, nm)
  }
}

# Calculate only NEW covariates
if (length(new_covs) > 0) {
  cat("Calculating", length(new_covs), "NEW covariates...\n")
  
  results_list <- lapply(seq_len(nlyr(basin_stack)), function(i) {
    
    basin_mask <- basin_stack[[i]]
    basin_id <- names(basin_mask)
    
    # MASK to valid watershed area (removes edge artifacts)
    basin_mask <- mask(basin_mask, watershed_mask)
    
    n_cells <- sum(!is.na(values(basin_mask)))
    basin_area_km2 <- (n_cells * cell_area_m2) / 1e6
    
    row <- data.frame(reach_id = basin_id, basin_area_km2 = basin_area_km2)
    
    # Only loop over NEW covariates
    for (nm in new_covs) {
      cov <- covariate_registry[[nm]]
      
      if (cov$type == "percent_cover") {
        row[[paste0("pct_", cov$out_name, "_bas")]] <-
          summarize_percent_cover_basin(basin_mask, aligned_rasters[[nm]])
        
      } else if (cov$type == "continuous") {
        row[[paste0(cov$out_name, "_bas")]] <-
          summarize_continuous_mean_basin(basin_mask, aligned_rasters[[nm]])
        
      } else if (cov$type == "point_density") {
        pd <- summarize_point_density_basin(basin_mask, point_layers[[nm]])
        row[[paste0(cov$out_name, "_count_bas")]] <- pd$n_points
        row[[paste0(cov$out_name, "_density_km2_bas")]] <- pd$density_km2
      }
    }
    row
  })
  
  results_new <- do.call(rbind, results_list) |>
    mutate(reach_id = as.numeric(gsub("basin_", "", reach_id))) |>
    filter(!is.na(reach_id)) |>
    select(reach_id, everything())
  
  # Merge old + new
  if (!is.null(results_old)) {
    results <- left_join(results_old, results_new, by = "reach_id")
  } else {
    results <- results_new
  }
  
  cat("✓ Calculated NEW covariates.\n")
  
} else {
  cat("No NEW covariates. Using existing results.\n")
  results <- results_old
}

summary(results)

write_csv(results, results_cache)
saveRDS(results, "../data/RVsummary_basin_covariates.rds")

# ---- 7. RCA-side: shared cell count -----------------------------------------

# Mask reach_rca_r to valid watershed area
reach_rca_masked <- mask(reach_rca_r, watershed_mask)

cell_count <- zonal(!is.na(reach_rca_masked), reach_rca_masked, fun = "sum") |>
  rename(reach_id = 1, cell_count = 2) |>
  mutate(rca_area_km2 = (cell_count * cell_area_m2) / 1e6)

# ---- 8. RCA-side loop -------
# Point layers still need rasterizing to a count-per-cell raster before zonal()
# can sum them, matching the original npdes_r_aligned approach.

rca_results <- cell_count

for (nm in names(covariate_registry)) {
  cov <- covariate_registry[[nm]]
  
  # Skip if this covariate is not a raster (e.g., npdes is point_density)
  if (is.null(aligned_rasters[[nm]])) {
    next
  }
  
  # Mask data raster to watershed
  data_masked <- mask(aligned_rasters[[nm]], watershed_mask)
  
  if (cov$type == "continuous") {
    stat <- zonal(data_masked, reach_rca_masked, fun = "mean", na.rm = TRUE) |>
      rename(reach_id = 1, value = 2) |>
      rename(!!paste0(cov$out_name, "_rca") := value)
    rca_results <- left_join(rca_results, stat, by = "reach_id")
    
  } else if (cov$type == "percent_cover") {
    covered_count <- zonal(data_masked, reach_rca_masked, fun = "sum", na.rm = TRUE) |>
      rename(reach_id = 1, covered_count = 2) |>
      mutate(covered_count = ifelse(is.na(covered_count), 0, covered_count))
    stat <- left_join(covered_count, cell_count, by = "reach_id") |>
      mutate(!!paste0("pct_", cov$out_name, "_rca") := (covered_count / cell_count) * 100) |>
      select(reach_id, !!paste0("pct_", cov$out_name, "_rca"))
    rca_results <- left_join(rca_results, stat, by = "reach_id")
  }
}

# Handle point_density separately (npdes)
for (nm in names(point_layers)) {
  cov <- covariate_registry[[nm]]
  pts_r <- rasterize(point_layers[[nm]], template, fun = "length", background = 0)
  pts_r_masked <- mask(pts_r, watershed_mask)
  pt_count <- zonal(pts_r_masked, reach_rca_masked, fun = "sum", na.rm = TRUE) |>
    rename(reach_id = 1, n_pts = 2) |>
    mutate(n_pts = ifelse(is.na(n_pts), 0, n_pts))
  stat <- left_join(pt_count, cell_count, by = "reach_id") |>
    mutate(
      !!paste0(cov$out_name, "_count_rca") := n_pts,
      !!paste0(cov$out_name, "_density_km2_rca") := (n_pts / (cell_count * cell_area_m2)) * 1e6
    ) |>
    select(reach_id, contains(cov$out_name))
  rca_results <- left_join(rca_results, stat, by = "reach_id")
}

summary(rca_results)

write_csv(rca_results, "../data/RVsummary_RCA_covariates.csv")
saveRDS(rca_results, "../data/RVsummary_RCA_covariates.rds")


# Get reach_id from basin_stack layer name
basin_223_layer <- which(grepl("basin_223", names(basin_stack)))
reach_id_from_basename <- as.numeric(gsub("basin_", "", names(basin_stack)[basin_223_layer]))

cat("Basin layer 223 has reach_id:", reach_id_from_basename, "\n")

# Extract all unique values in reach_rca_r
all_rca_ids <- sort(unique(values(reach_rca_r), na.rm = TRUE))
cat("Total unique:", length(all_rca_ids), "\n")

# Which zone has the most cells?
cell_counts <- table(values(reach_rca_r))
cell_counts_sorted <- sort(cell_counts, decreasing = TRUE)
cat("\nTop 10 zones by cell count:\n")
print(head(cell_counts_sorted, 10))

# ---- 9. Visualize -----------------------------------------------------------

results |>
  pivot_longer(cols = -reach_id, names_to = "covariate", values_to = "value") |>
  ggplot(aes(x = covariate, y = value)) +
  geom_violin() +
  scale_y_log10() +
  facet_wrap(~covariate, scales = "free", nrow = 1)

rca_results |>
  pivot_longer(cols = -reach_id, names_to = "covariate", values_to = "value") |>
  ggplot(aes(x = covariate, y = value)) +
  geom_violin() +
  scale_y_log10() +
  facet_wrap(~covariate, scales = "free", nrow = 1)


# ---- 10. Visualization: Spatial distribution of predictors ----------------
library(tmap)
library(purrr)

# ---- Setup spatial layers ----
reach_rca_sf <- as.polygons(reach_rca_r) |> 
  st_as_sf() |>
  rename(reach_id = reach_rca_multipoint) |>
  left_join(rca_results, by = "reach_id")

basin_summary_sf <- map_df(seq_len(nlyr(basin_stack)), function(i) {
  layer_name <- names(basin_stack)[i]
  reach_id_val <- as.numeric(gsub("basin_", "", layer_name))
  
  as.polygons(basin_stack[[i]]) |> 
    st_as_sf() |> 
    select(geometry) |>
    mutate(reach_id = reach_id_val)
}) |>
  filter(!is.na(reach_id)) |>
  left_join(results, by = "reach_id") |>
  st_as_sf() |>
  arrange(desc(basin_area_km2))  # <- Large basins drawn first, small ones on top

# ---- Define covariates  ----
covariates_to_plot <- tibble(
  raster_nm = c("nlcd_imp2025", "nlcd_cultivated_crops_2025", "nlcd_pasture_hay_2025", 
                "nlcd_planted_cultivated_2025", "bfi", "organic_matter", "clay", "ksat", "pop2025"),
  rca_col = c("impervious_2025_rca", "pct_crops_2025_rca", "pct_pasture_hay_2025_rca",
              "pct_crops_and_pasture_2025_rca", "bfi_rca", "om_rca", "clay_pct_rca", 
              "ksat_rca", "pop_density_2025_rca"),
  bas_col = c("impervious_2025_bas", "pct_crops_2025_bas", "pct_pasture_hay_2025_bas",
              "pct_crops_and_pasture_2025_bas", "bfi_bas", "om_bas", "clay_pct_bas", 
              "ksat_bas", "pop_density_2025_bas"),
  title = c("Impervious Cover (%)", "Cultivated Crops (%)", "Pasture/Hay (%)",
            "Planted/Cultivated (%)", "Baseflow Index", "Organic Matter (%)", 
            "Clay Content (%)", "Hydraulic Conductivity (Ksat)", "Population Density"),
  palette = c("Purples", "YlOrBr", "YlGn", "Oranges", "Blues", "Greens", "YlGn", "RdYlBu", "Reds")
)

# ---- Create 3-panel figures ----

make_3panel_map <- function(cov_row) {
  
  rca_vals <- na.omit(reach_rca_sf[[cov_row$rca_col]])
  bas_vals <- na.omit(basin_summary_sf[[cov_row$bas_col]])
  all_vals <- c(rca_vals, bas_vals)
  val_min <- min(all_vals, na.rm = TRUE)
  val_max <- max(all_vals, na.rm = TRUE)
  
  # Fewer breaks (5 instead of 7) and rounded labels
  breaks <- round(seq(val_min, val_max, length.out = 5), 2)
  labels <- format(breaks, big.mark = ",", scientific = FALSE, trim = TRUE)
  
  # Shared fill scale with clean labels and short legend title
  shared_scale <- scale_fill_distiller(
    palette = cov_row$palette,
    limits = c(val_min, val_max),
    breaks = breaks,
    labels = labels,
    name = cov_row$title,   # Clean legend title
    na.value = "grey80",
    direction = 1
  )
  
  # Panel 1: Raw raster
  p1 <- ggplot() +
    geom_spatraster(data = aligned_rasters[[cov_row$raster_nm]]) +
    geom_sf(data = reach_rca_sf, fill = NA, color = "black", linewidth = 0.3) +
    scale_fill_distiller(
      palette = cov_row$palette,
      na.value = "transparent",
      direction = 1,
      name = "Raw Value",
      labels = function(x) format(round(x, 2), scientific = FALSE, trim = TRUE)
    ) +
    labs(title = "A. Raw Raster") +
    theme_void() +
    theme(legend.position = "bottom",
          legend.key.width = unit(1.0, "cm"),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 8))
  
  # Shared theme for legends on p2 and p3
  shared_legend_theme <- theme(
    legend.position = "bottom",
    legend.key.width = unit(0.7, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 7, angle = 60, hjust = 1)
  )
  #reach scale
  p2 <- ggplot() +
    geom_sf(data = reach_rca_sf, 
            aes(fill = .data[[cov_row$rca_col]]),
            color = "black", linewidth = 0.3) +
    shared_scale +
    labs(title = "B. Reach-scale Summary") +
    theme_void() +
    shared_legend_theme
  #basin scale plot
  p3 <- ggplot() +
    geom_sf(data = basin_summary_sf, 
            aes(fill = .data[[cov_row$bas_col]]),
            color = "black", linewidth = 0.3) +
    shared_scale +
    labs(title = "C. Basin-scale Summary") +
    theme_void() +
    shared_legend_theme
  
  combined <- p1 + p2 + p3 + 
    plot_layout(ncol = 3) +
    plot_annotation(title = cov_row$title)
  
  ggsave(
    paste0("../figs/Supplementary_", cov_row$raster_nm, "_3panel.png"),
    combined, width = 16, height = 5, dpi = 300)
  
  cat("✓ Saved:", cov_row$raster_nm, "\n")
  invisible(NULL)
}

for (i in 1:nrow(covariates_to_plot)) {
  make_3panel_map(covariates_to_plot[i, ])
}

