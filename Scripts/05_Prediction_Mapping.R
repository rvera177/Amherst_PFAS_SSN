#ok we finished fitting the models in Script 04_PFAS_SSN_model_selection!!
#Now, moving on to plotting all the model fits onto the NHDV2-HR flowlines

# ---- 1. Setup and load objects from script 04 ----------------------------------
library(sf)
library(tidyverse)
library(SSN2)
library(patchwork)
library(wesanderson)
library(scales)
# Run from the project root (the folder that contains data/ and scripts/); all paths are relative to data/
if (basename(getwd()) != "data") setwd("data")
source("../scripts/00 PFAS naming.R")  # PFAS compounds, functional groups and column-naming rules

for (p in c("ggridges", "ggdist", "ggrepel")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Please install.packages('", p, "')")
}

# Use the objects already in the session, otherwise read what script 04 saved
for (nm in c("ssn_obj_list", "fit_model_ls_list", "best_models_list",
             "final_key_compounds_list", "diagnostics_plot_ls_list")) {
  if (!exists(nm)) assign(nm, readRDS(file.path("../outputs", paste0(nm, ".RDS"))))
}

out_dir <- "../figs/PFAS_maps_HR"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- 2.1. Compute AIC-weighted average predictions per compound -----
# For each compound: predict at each "preds" point from every
# retained model in its suite, weight each model's fitted value by
# its AIC_weight_i (already summing to 1 within a compound's set),
# and sum to get the model-averaged prediction.
join_key <- "reach_id"
# What variables were used in the models?
cat("\nVariables in first compound's first model:\n")
print(all.vars(fit_model_ls_list[[1]]$fit_model_ls[[1]][[1]]$formula))

# What variables are in preds?
cat("\nVariables in preds:\n")
print(names(ssn_obj_list[[1]]$preds_covariates_std))


weighted_preds_all_list <- fit_model_ls_list |>
  purrr::imap(function(ev, event_name) {
    
    fit_model_ls         <- ev$fit_model_ls
    best_PFAS_models     <- best_models_list[[event_name]]
    preds_covariates_std <- ssn_obj_list[[event_name]]$preds_covariates_std
    
    # Helper function to get weighted predictions for one compound
    get_weighted_preds <- function(compound) {
      mods    <- fit_model_ls[[compound]]
      weights <- best_PFAS_models[[compound]]$AIC_weight_i  # same row order as mods
      
      weighted_long <- purrr::map2_dfr(mods, weights, function(m, w) {
        SSN2::augment(m, newdata = "preds", pred.type = "preds") |>
          sf::st_drop_geometry() |>
          dplyr::transmute(pid = as.numeric(pid), weighted_pred = .fitted * w)
      })
      
      weighted_long |>
        dplyr::group_by(pid) |>
        dplyr::summarise(pred_val = sum(weighted_pred, na.rm = TRUE), .groups = "drop") |>
        dplyr::mutate(compound = compound)
    }
    
    key_compounds <- names(fit_model_ls)  # all compounds for this event
    
    weighted_preds_event <- purrr::map_dfr(key_compounds, get_weighted_preds)
    
    # Attach the join key (reach_id) from the prediction points
    pred_id_lookup <- preds_covariates_std |>
      sf::st_drop_geometry() |>
      dplyr::select(pid, dplyr::all_of(join_key)) |>
      dplyr::mutate(pid = as.numeric(pid))
    
    weighted_preds_event <- weighted_preds_event |>
      dplyr::left_join(pred_id_lookup, by = "pid") |>
      dplyr::mutate(event_name = event_name)
    weighted_preds_event
  })

# Combine all events
weighted_preds_all <- dplyr::bind_rows(weighted_preds_all_list)

# Wide format, with one column per compound, one row per join_key for each sampling event
weighted_preds_wide <- weighted_preds_all |>
  dplyr::select(-pid) |>
  tidyr::pivot_wider(
    names_from = compound, 
    values_from = pred_val,
    names_glue = "{compound}_pred",
    id_cols = c(event_name, dplyr::all_of(join_key))) |>
  dplyr::distinct(event_name, dplyr::across(dplyr::all_of(join_key)), .keep_all = TRUE)

cat("Weighted predictions wide format:", nrow(weighted_preds_wide), "rows\n")
head(weighted_preds_wide)

# ---- 2.2. Join weighted predictions onto edges per event ---------
edges_with_preds_list <- ssn_obj_list |>
  purrr::imap(function(ssn_data, event_name) {
    
    ssn_obj_obspreds_std <- ssn_data$ssn_obj_obspreds_std
    
    # Get predictions for this event only
    event_preds <- weighted_preds_wide |>
      dplyr::filter(event_name == !!event_name) |>
      dplyr::select(-event_name)  # Remove event_name col for join
    
    if (nrow(event_preds) == 0) {
      cat("Event ", event_name, " — no predictions to join\n")
      return(NULL)
    }
    
    # Join predictions onto edges for this event
    edges_with_preds <- ssn_obj_obspreds_std$edges |>
      dplyr::left_join(event_preds, by = join_key)
    
    pred_cols_present <- grep("_pred$", names(edges_with_preds), value = TRUE)
    cat("Event ", event_name, " — pred columns joined:", length(pred_cols_present), "\n")
    
    edges_with_preds
  })

# Combine edges from all events
edges_with_preds <- dplyr::bind_rows(edges_with_preds_list)

pred_cols_present <- grep("_pred$", names(edges_with_preds), value = TRUE)
cat("\nTotal pred columns in combined edges:", length(pred_cols_present), "\n")

# ---- 2.3. Pre-compute shared scale_max per base_compound (always in ng/L) ----
# Takes the maximum 0.99 quantile across both raw and log1p versions so
# RAW_PFOA and LOG_PFOA maps are on an identical color scale.

shared_scale_max_list <- ssn_obj_list |>
  purrr::imap(function(ssn_data, event_name) {
    
    fit_model_ls      <- fit_model_ls_list[[event_name]]$fit_model_ls
    best_PFAS_models  <- best_models_list[[event_name]]
    obs_covariates_s  <- ssn_data$obs_covariates_s
    
    # Get edges with predictions for this event
    edges_event <- edges_with_preds_list[[event_name]]
    
    if (is.null(edges_event)) {
      cat("Event ", event_name, " — no edges to process\n")
      return(list())
    }
    shared_scale_max <- list()
    
    for (compound in names(fit_model_ls)) {
      is_log1p      <- startsWith(compound, "log1p_")
      base_compound <- ifelse(is_log1p, sub("^log1p_", "", compound), compound)
      pred_col      <- paste0(compound, "_pred")
      
      if (!pred_col %in% names(edges_event)) next
      
      preds <- edges_event[[pred_col]]
      if (is_log1p) preds <- expm1(preds)   # back-transform to ng/L
      
      obs_vals <- obs_covariates_s[[base_compound]]
      vals <- c(preds, obs_vals)
      qval <- quantile(vals, 0.99, na.rm = TRUE)
      
      # Keep the larger of the two (raw vs log1p)
      if (is.null(shared_scale_max[[base_compound]])) {
        shared_scale_max[[base_compound]] <- qval
      } else {
        shared_scale_max[[base_compound]] <- max(shared_scale_max[[base_compound]], qval)
      }
    }
    
    cat("Event ", event_name, " — computed scale_max for ", 
        length(shared_scale_max), " compounds\n", sep = "")
    shared_scale_max
  })

# Example: access for event S3
cat("\nShared scale_max for S3:\n")
print(unlist(shared_scale_max_list[["S3"]]))

# ---- 3. Per-compound maps for each event --------

ssn_obj_list |>
  purrr::imap(function(ssn_data, event_name) {
    
    fit_model_ls      <- fit_model_ls_list[[event_name]]$fit_model_ls
    best_PFAS_models  <- best_models_list[[event_name]]
    obs_covariates_s  <- ssn_data$obs_covariates_s
    shared_scale_max  <- shared_scale_max_list[[event_name]]
    
    edges_event <- edges_with_preds_list[[event_name]]
    
    if (is.null(edges_event)) {
      cat("Event ", event_name, " — no edges to map\n")
      return(NULL)
    }
    
    # Create event-specific output directory
    out_dir_event <- file.path(out_dir, event_name)
    if (!dir.exists(out_dir_event)) dir.create(out_dir_event, showWarnings = FALSE)
    
    pal <- colorRampPalette(c("#012A4A", "#014F86",
                              wes_palette("Zissou1", type = "continuous"),
                              "#8B0000"))(256)
    
    # Get all response vars (both raw and log1p)
    all_response_vars <- c(
      names(fit_model_ls)[!startsWith(names(fit_model_ls), "log1p_")],
      names(fit_model_ls)[startsWith(names(fit_model_ls), "log1p_")]) |> unique()
    
    for (compound in all_response_vars) {
      pred_col <- paste0(compound, "_pred")
      
      # Skip if no prediction for this compound
      if (!pred_col %in% names(edges_event)) {
        next
      }
      
      # Scale detection and display setup -
      is_log1p      <- startsWith(compound, "log1p_")
      base_compound <- ifelse(is_log1p, sub("^log1p_", "", compound), compound)
      obs_col       <- base_compound
      prefix        <- ifelse(is_log1p, "LOG", "RAW")
      
      # Back-transform log1p predictions to ng/L for display
      plot_edges <- edges_event
      if (is_log1p) {
        plot_edges[[pred_col]] <- expm1(plot_edges[[pred_col]])
      }
      
      scale_min <- 0
      scale_max <- shared_scale_max[[base_compound]]
      
      if (is.na(scale_max)) {
        message("Skipping ", compound, " — no scale_max computed")
        next
      }
      
      best_model_info <- best_PFAS_models[[compound]] |> slice(1)
      pseudo_r2  <- best_model_info$pseudo.r.squared
      predictors <- best_model_info$predictor_var
      
      scale_label <- ifelse(is_log1p, "ng/L\n(back-transformed)", "ng/L")
      
      subtitle <- paste0(
        "[", prefix, " model]  ",
        "R\u00b2 = ", round(pseudo_r2, 3),
        "\nPredictors: ", predictors)
      
      p <- ggplot() +
        geom_sf(data = plot_edges,
                aes(color = .data[[pred_col]]), linewidth = 1.5) +
        geom_sf(data = obs_covariates_s,
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
        labs(title    = paste0("[", prefix, "] Predicted ", base_compound),
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
          axis.text.y     = element_text(size = 9))
      
      full_path <- file.path(out_dir_event, paste0(prefix, "_", base_compound, "_HR_map.png"))
      ggsave(filename = full_path, plot = p, width = 8, height = 6, dpi = 300)
      
      if (file.exists(full_path)) {
        message("Saved: ", basename(full_path))
      } else {
        warning("File NOT found after save: ", full_path)
      }
    }
    
    cat("Event ", event_name, " — maps completed\n", sep = "")
    NULL
  })

# ---- 4. Obs vs Predicted scatter plots with R², AIC, RMSPE, Bias ----

ssn_obj_list |>
  purrr::imap(function(ssn_data, event_name) {
    
    best_PFAS_models    <- best_models_list[[event_name]]
    diagnostics_plot_ls <- diagnostics_plot_ls_list[[event_name]]
    final_key_compounds <- final_key_compounds_list[[event_name]]$final_key_compounds
    scale_comparison_df <- final_key_compounds_list[[event_name]]$scale_comparison_df
    
    out_dir_event <- file.path(out_dir, event_name)
    if (!dir.exists(out_dir_event)) dir.create(out_dir_event, showWarnings = FALSE)
    
    for (compound in final_key_compounds) {
      
      is_log1p      <- startsWith(compound, "log1p_")
      base_compound <- ifelse(is_log1p, sub("^log1p_", "", compound), compound)
      prefix        <- ifelse(is_log1p, "LOG", "RAW")
      
      obspreds_raw <- diagnostics_plot_ls[[compound]]
      if (is.null(obspreds_raw)) next
      
      # ---- Pull AIC for BOTH scales of this compound ----------------------
      aic_raw <- tryCatch(
        best_PFAS_models[[base_compound]] |>
          dplyr::slice(1) |> dplyr::pull(AIC),
        error = \(e) NA_real_)
      
      aic_log <- tryCatch(
        best_PFAS_models[[paste0("log1p_", base_compound)]] |>
          dplyr::slice(1) |> dplyr::pull(AIC),
        error = \(e) NA_real_)
      
      # Pull RMSPE and Bias
      scale_row <- scale_comparison_df |>
        dplyr::filter(compound == base_compound)
      
      if (nrow(scale_row) == 0) next
      
      rmspe_display <- if (is_log1p) {
        scale_row |> dplyr::pull(rmspe_log1p_ng_per_L)
      } else {
        scale_row |> dplyr::pull(rmspe_raw)
      }
      
      bias_display <- if (is_log1p) {
        scale_row |> dplyr::pull(bias_log1p)
      } else {
        scale_row |> dplyr::pull(bias_raw)
      }
      
      rmspe_label <- if (is_log1p) "RMSE = " else "LOOCV RMSPE = "
      
      #Get R² and n observations
      pseudo_r2 <- tryCatch(
        best_PFAS_models[[compound]] |> dplyr::slice(1) |> dplyr::pull(pseudo.r.squared),
        error = \(e) NA_real_)
      
      model_avg_df <- obspreds_raw$obspreds_df |>
        dplyr::filter(pred_type == "AIC-weighted predictions")
      
      n_obs <- nrow(model_avg_df)
      
      #annotation with R², n, AIC, RMSPE, Bias
      annot_label <- paste0(
        "R\u00b2 = ",                     format(round(pseudo_r2, 3), nsmall = 3),     "\n",
        "n = ",                           n_obs,                                       "\n",
        "AIC\u2081\u209a = ",             format(round(aic_log, 1), nsmall = 1),       "\n",
        "AIC\u1d63\u2090\u1d67 = ",       format(round(aic_raw, 1), nsmall = 1),       "\n",
        rmspe_label,                      format(round(rmspe_display, 2), nsmall = 2), " ng/L\n",
        "Bias = ",                        format(round(bias_display, 2), nsmall = 2),  " ng/L"
      )
      
      # Back-transform if log1p 
      if (is_log1p) {
        model_avg_df <- model_avg_df |>
          dplyr::mutate(observations = expm1(observations),
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
      
      full_path <- file.path(out_dir_event,
                             paste0(prefix, "_", base_compound, "_obspred_scatter.png"))
      ggsave(filename = full_path, plot = p, width = 6, height = 6, dpi = 300)
      
      if (file.exists(full_path)) {
        message("Saved: ", basename(full_path))
      } else {
        warning("File NOT found: ", full_path)
      }
    }
    
    NULL
  })   

# ---- 5. Assess leverage and influence of extreme sites per event ----

ssn_obj_list |>
  purrr::imap(function(ssn_data, event_name) {
    
    fit_model_ls        <- fit_model_ls_list[[event_name]]$fit_model_ls
    diagnostics_plot_ls <- diagnostics_plot_ls_list[[event_name]]
    
    final_key_compounds <- final_key_compounds_list[[event_name]]$final_key_compounds
    
    cat("\n════ Event: ", event_name, " ════\n", sep = "")
    
    # Process each final key compound
    for (compound in final_key_compounds) {
      
      diagnostics <- diagnostics_plot_ls[[compound]]
      
      if (is.null(diagnostics)) {
        message("  No diagnostics for ", compound)
        next
      }
      
      obspreds_df <- diagnostics$obspreds_df |>
        dplyr::filter(pred_type == "AIC-weighted predictions")
      
      if (nrow(obspreds_df) == 0) {
        message("  No obs/pred data for ", compound)
        next
      }
      
      # Calculate Cook's Distance threshold
      n <- nrow(obspreds_df)
      cooksd_threshold <- 4 / n
      
      # Identify high-leverage points
      high_leverage <- obspreds_df |>
        dplyr::filter(.hat > 0.3) |>
        dplyr::arrange(desc(.hat)) |>
        dplyr::mutate(compound = compound, event = event_name)
      
      if (nrow(high_leverage) > 0) {
        cat("\n  ", compound, " — High leverage points (.hat > 0.3):\n", sep = "")
        print(high_leverage |> 
                dplyr::select(pid, site_id, observations, predictions, .hat, .cooksd) |>
                head(10))
      }
    }
    
    NULL
  })


# ---- 6. Model performance by group ----

# ---- 6.1 Model Performance by Group (R² + CV-RMSPE) ----

ssn_obj_list |>
  purrr::imap(function(ssn_data, event_name) {
    
    best_PFAS_models    <- best_models_list[[event_name]]
    final_key_compounds <- final_key_compounds_list[[event_name]]$final_key_compounds
    scale_comparison_df <- final_key_compounds_list[[event_name]]$scale_comparison_df
    
    out_dir_event <- file.path(out_dir, event_name)
    if (!dir.exists(out_dir_event)) dir.create(out_dir_event, showWarnings = FALSE)
    
    # Build performance_df
    performance_df <- purrr::map_dfr(
      final_key_compounds,
      function(cmp) {
        best_row <- best_PFAS_models[[cmp]] |> dplyr::slice(1)
        is_log1p      <- startsWith(cmp, "log1p_")
        base_compound <- ifelse(is_log1p, sub("^log1p_", "", cmp), cmp)
        
        data.frame(
          compound     = base_compound,
          response_var = cmp,
          better_scale = ifelse(is_log1p, "log1p", "raw"),
          pseudo_r2    = best_row$pseudo.r.squared,
          AIC          = best_row$AIC,
          stringsAsFactors = FALSE)
      }
    ) |>
      dplyr::left_join(pfas_family_lookup, by = "compound") |>
      dplyr::left_join(
        # join on compound only: scale_comparison_df$better_scale is the value BEFORE the
        # "within 5% -> raw" rule, so joining on it silently gave NA CV-RMSPE for those compounds
        scale_comparison_df |>
          dplyr::select(compound, rmspe_raw, rmspe_log1p_ng_per_L, cv_rmspe_raw, cv_rmspe_log1p),
        by = "compound") |>
      dplyr::mutate(
        cv_rmspe_final = ifelse(better_scale == "log1p", cv_rmspe_log1p, cv_rmspe_raw))
    
    # Order by R²
    compound_order <- performance_df |>
      dplyr::arrange(pseudo_r2) |> dplyr::pull(compound)
    
    performance_df <- performance_df |>
      dplyr::mutate(compound = factor(compound, levels = compound_order))
    
    # R² panel
    r2_panel <- ggplot(performance_df, aes(x = compound, y = pseudo_r2,
                                           color = family, shape = better_scale)) +
      geom_point(size = 3.5) +
      coord_flip() +
      scale_color_brewer(palette = "Set1") +
      labs(x = NULL, y = "Pseudo R² (winning model)") +
      theme_bw(base_size = 11) +
      theme(legend.position = "none")
    
    # CV-RMSPE panel
    cv_rmspe_panel <- ggplot(performance_df, aes(x = compound, y = cv_rmspe_final,
                                                 color = family, shape = better_scale)) +
      geom_point(size = 3.5) +
      coord_flip() +
      scale_color_brewer(palette = "Set1") +
      scale_y_continuous(labels = scales::percent) +
      labs(x = NULL, y = "CV-RMSPE\n(error as % of mean)") +
      theme_bw(base_size = 11) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    
    # Combine
    combined_perf <- r2_panel + cv_rmspe_panel +
      patchwork::plot_layout(guides = "collect") +
      patchwork::plot_annotation(
        title = paste0("[", event_name, "] Model fit and normalized error"))
    
    full_path <- file.path(out_dir_event, "ModelPerformance_R2_CVRMSPE.png")
    ggsave(combined_perf, filename = full_path, width = 11, height = 9, dpi = 300)
    message("Saved: ", basename(full_path))
    NULL
  })


# ---- 6.2 Ridge Plot ----

ssn_obj_list |>
  purrr::imap(function(ssn_data, event_name) {
    
    final_key_compounds <- final_key_compounds_list[[event_name]]$final_key_compounds
    diagnostics_plot_ls <- diagnostics_plot_ls_list[[event_name]]
    out_dir_event <- file.path(out_dir, event_name)
    
    # Build obs/pred data
    obspred_dist_df <- purrr::map_dfr(
      final_key_compounds,
      function(cmp) {
        dat <- diagnostics_plot_ls[[cmp]]
        if (is.null(dat)) return(NULL)
        dat <- dat$obspreds_df |> dplyr::filter(pred_type == "AIC-weighted predictions")
        is_log1p <- startsWith(cmp, "log1p_")
        base_compound <- ifelse(is_log1p, sub("^log1p_", "", cmp), cmp)
        if (is_log1p) {
          dat <- dat |> dplyr::mutate(observations = expm1(observations),
                                      predictions  = expm1(predictions))
        }
        dplyr::bind_rows(
          data.frame(compound = base_compound, value = dat$observations, type = "Observed"),
          data.frame(compound = base_compound, value = dat$predictions, type = "Predicted"))
      }
    ) |>
      dplyr::left_join(pfas_family_lookup, by = "compound") |>
      dplyr::mutate(value = pmax(value, 0))
    
    if (nrow(obspred_dist_df) == 0) return(NULL)
    
    ridge_plot <- ggplot(obspred_dist_df, aes(x = value + 0.001, y = compound, fill = type, color = type)) +
      ggridges::geom_density_ridges(alpha = 0.45, scale = 0.85, rel_min_height = 0.005, linewidth = 0.3) +
      scale_x_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                    breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000)) +
      scale_fill_manual(values = c("Observed" = "#2166AC", "Predicted" = "#D6604D"), name = NULL) +
      scale_color_manual(values = c("Observed" = "#2166AC", "Predicted" = "#D6604D"), name = NULL) +
      labs(x = "Concentration (ng/L, log scale)", y = NULL,
           title = paste0("[", event_name, "] Observed vs predicted distributions")) +
      theme_bw(base_size = 11) +
      theme(legend.position = "top", axis.text.y = element_text(size = 8),
            panel.grid.minor = element_blank())
    
    full_path <- file.path(out_dir_event, "RidgePlot_ObsPred.png")
    ggsave(ridge_plot, filename = full_path, width = 10, height = 12, dpi = 300)
    message("Saved: ", basename(full_path))
    NULL
  })


# ---- 6.3 Raincloud Plots (by family) ----

ssn_obj_list |>
  purrr::imap(function(ssn_data, event_name) {
    
    final_key_compounds <- final_key_compounds_list[[event_name]]$final_key_compounds
    diagnostics_plot_ls <- diagnostics_plot_ls_list[[event_name]]
    out_dir_event <- file.path(out_dir, event_name)
    
    # Build data
    obspred_dist_df <- purrr::map_dfr(
      final_key_compounds,
      function(cmp) {
        dat <- diagnostics_plot_ls[[cmp]]
        if (is.null(dat)) return(NULL)
        dat <- dat$obspreds_df |> dplyr::filter(pred_type == "AIC-weighted predictions")
        is_log1p <- startsWith(cmp, "log1p_")
        base_compound <- ifelse(is_log1p, sub("^log1p_", "", cmp), cmp)
        if (is_log1p) {
          dat <- dat |> dplyr::mutate(observations = expm1(observations),
                                      predictions  = expm1(predictions))
        }
        dplyr::bind_rows(
          data.frame(compound = base_compound, value = dat$observations, type = "Observed"),
          data.frame(compound = base_compound, value = dat$predictions, type = "Predicted"))
      }
    ) |>
      dplyr::left_join(pfas_family_lookup, by = "compound") |>
      dplyr::mutate(value = pmax(value, 0.001))
    
    if (nrow(obspred_dist_df) == 0) return(NULL)
    
    families <- unique(obspred_dist_df$family)
    families <- families[!is.na(families)]
    
    for (fam in families) {
      fam_data <- obspred_dist_df |> dplyr::filter(family == fam)
      
      raincloud_plot <- ggplot(fam_data, aes(x = value + 0.001, y = compound, fill = type, color = type)) +
        ggdist::stat_halfeye(aes(slab_alpha = after_stat(pdf)), adjust = 0.7, width = 0.6,
                             point_color = NA, .width = 0, justification = -0.2) +
        geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.4,
                     position = position_nudge(y = -0.15)) +
        ggdist::stat_dots(side = "left", dotsize = 0.5, justification = 1.1,
                          binwidth = 0.05, position = position_nudge(y = -0.15)) +
        scale_x_log10(labels = scales::comma) +
        scale_fill_manual(values = c("Observed" = "#2166AC", "Predicted" = "#D6604D")) +
        scale_color_manual(values = c("Observed" = "#2166AC", "Predicted" = "#D6604D")) +
        coord_flip() +
        labs(x = "Concentration (ng/L, log scale)", y = NULL,
             title = paste0("[", event_name, "] ", fam),
             fill = NULL, color = NULL) +
        theme_bw(base_size = 12) +
        theme(legend.position = "top")
      
      full_path <- file.path(out_dir_event, paste0("Raincloud_", gsub(" ", "_", fam), ".png"))
      ggsave(raincloud_plot, filename = full_path, width = 8, height = 6, dpi = 300)
      message("  Saved: ", basename(full_path))
    }
    NULL
  })


# ---- 6.4 PCA Biplot ----

ssn_obj_list |>
  purrr::imap(function(ssn_data, event_name) {
    
    best_PFAS_models    <- best_models_list[[event_name]]
    final_key_compounds <- final_key_compounds_list[[event_name]]$final_key_compounds
    out_dir_event <- file.path(out_dir, event_name)
    
    # Build importance matrix
    rel_imp_ls <- purrr::map_dfr(final_key_compounds, function(cmp) {
      best_row <- best_PFAS_models[[cmp]] |> dplyr::slice(1)
      pred_str <- best_row$predictor_var
      covariates <- strsplit(pred_str, " + ", fixed = TRUE)[[1]]
      data.frame(response_var = cmp, covariate = covariates,
                 rel_imp = 1 / length(covariates), stringsAsFactors = FALSE)
    })
    
    if (nrow(rel_imp_ls) == 0) return(NULL)
    
    rel_imp_matrix <- rel_imp_ls |>
      tidyr::pivot_wider(names_from = covariate, values_from = rel_imp, values_fill = 0) |>
      tibble::column_to_rownames("response_var") |> as.matrix()
    
    pca_res <- prcomp(rel_imp_matrix, scale. = TRUE)
    scores <- as.data.frame(pca_res$x) |>
      tibble::rownames_to_column("response_var") |>
      dplyr::mutate(compound = sub("^log1p_", "", response_var)) |>
      dplyr::left_join(pfas_family_lookup, by = "compound")
    
    loadings <- as.data.frame(pca_res$rotation) |> tibble::rownames_to_column("covariate")
    
    arrow_scale <- min(
      (max(scores$PC1) - min(scores$PC1)) / (max(loadings$PC1) - min(loadings$PC1)),
      (max(scores$PC2) - min(scores$PC2)) / (max(loadings$PC2) - min(loadings$PC2))) * 0.8
    
    var_explained <- summary(pca_res)$importance[2, 1:2] * 100
    
    biplot_p <- ggplot() +
      geom_hline(yintercept = 0, color = "grey85") +
      geom_vline(xintercept = 0, color = "grey85") +
      geom_segment(data = loadings,
                   aes(x = 0, y = 0, xend = PC1 * arrow_scale, yend = PC2 * arrow_scale),
                   arrow = arrow(length = unit(0.2, "cm")), color = "grey40") +
      ggrepel::geom_text_repel(data = loadings,
                               aes(PC1 * arrow_scale, PC2 * arrow_scale, label = covariate),
                               color = "grey30", size = 3) +
      geom_point(data = scores, aes(PC1, PC2, color = family), size = 3) +
      ggrepel::geom_text_repel(data = scores, aes(PC1, PC2, label = compound, color = family),
                               size = 3, show.legend = FALSE) +
      scale_color_brewer(palette = "Paired") +
      labs(x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
           y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
           color = "PFAS family",
           title = paste0("[", event_name, "] Covariate importance")) +
      theme_bw(base_size = 11) +
      coord_equal()
    
    full_path <- file.path(out_dir_event, "PCA_Biplot.png")
    ggsave(biplot_p, filename = full_path, width = 8, height = 7, dpi = 300)
    message("Saved: ", basename(full_path))
    NULL
  })

