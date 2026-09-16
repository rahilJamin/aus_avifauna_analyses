# =============================================================================
# pipeline_lcbd.R
# Targets for the Gaussian pool-size-standardised LCBD ensemble.
# =============================================================================

# Stage settings ----
lcbd_setting_names <- grep("^lcbd_", names(cfg), value = TRUE)
lcbd_settings <- cfg[setdiff(lcbd_setting_names, "lcbd_results_path")]
lcbd_settings$model_threads <- cfg$model_threads
lcbd_settings$noisy_miner_species <- cfg$noisy_miner_species

lcbd_paths <- cfg[c(
  "candidate_model_df_path",
  "species_matrix_path",
  "beta_matrices_path",
  "lcbd_results_path"
)]

# Target definitions ----

lcbd_targets <- list(
# 04 Gaussian LCBD ensemble ----

  tar_target(
    lcbd_inputs,
    load_lcbd_inputs(
      site_data_path = tracked_file_path(
        model_data_files,
        lcbd_paths$candidate_model_df_path
      ),
      species_matrix_path = tracked_file_path(
        inext_matrix_files,
        lcbd_paths$species_matrix_path
      ),
      hex_path = main_shapefile_path(hex_grid_files),
      focal_predictors = lcbd_settings$lcbd_focal_predictors,
      noisy_miner_species = lcbd_settings$noisy_miner_species,
      reference_pool_crs = lcbd_settings$lcbd_reference_pool_crs
    )
  ),
  tar_target(
    lcbd_distance_matrices,
    calculate_lcbd_beta_matrices(
      lcbd_inputs$lcbd_species_matrix
    )$dist_mats
  ),
  tar_target(
    lcbd_distance_validation,
    validate_lcbd_beta_matrices(
      saved_beta_path = tracked_file_path(
        model_data_files,
        lcbd_paths$beta_matrices_path
      ),
      rebuilt_distance_matrices = lcbd_distance_matrices,
      ordered_site_ids = lcbd_inputs$ordered_site_ids,
      tolerance = lcbd_settings$lcbd_matrix_validation_tolerance
    )
  ),
  tar_target(
    lcbd_gradient_thresholds,
    calculate_lcbd_gradient_thresholds(
      site_data = lcbd_inputs$site_data,
      focal_predictors = lcbd_settings$lcbd_focal_predictors,
      lower_tail_probability = lcbd_settings$lcbd_lower_tail_probability,
      upper_tail_probability = lcbd_settings$lcbd_upper_tail_probability
    )
  ),
  tar_target(
    lcbd_allocation_seeds,
    lcbd_settings$lcbd_allocation_seed_start +
      seq_len(lcbd_settings$lcbd_n_allocations) - 1L
  ),
  tar_target(
    lcbd_allocation_bundle,
    prepare_lcbd_reference_allocations(
      site_data = lcbd_inputs$site_data,
      distance_matrices = lcbd_distance_matrices,
      focal_predictors = lcbd_settings$lcbd_focal_predictors,
      global_gradient_thresholds = lcbd_gradient_thresholds,
      allocation_seeds = lcbd_allocation_seeds,
      reference_pool_k = lcbd_settings$lcbd_reference_pool_k,
      reference_pool_nstart = lcbd_settings$lcbd_reference_pool_nstart,
      minimum_pool_sites = lcbd_settings$lcbd_minimum_pool_sites,
      minimum_tail_sites = lcbd_settings$lcbd_minimum_tail_sites,
      minimum_informative_pools =
        lcbd_settings$lcbd_minimum_informative_pools,
      curve_lower_probability = lcbd_settings$lcbd_curve_lower_probability,
      curve_upper_probability = lcbd_settings$lcbd_curve_upper_probability,
      negative_tolerance = lcbd_settings$lcbd_negative_tolerance,
      lcbd_sum_tolerance = lcbd_settings$lcbd_sum_tolerance
    )
  ),
  tar_target(
    lcbd_formulas,
    define_lcbd_formulas()
  ),
  tar_target(
    lcbd_prediction_grids,
    make_lcbd_prediction_grids(
      common_curve_limits = lcbd_allocation_bundle$common_curve_limits,
      curve_grid_points = lcbd_settings$lcbd_curve_grid_points
    )
  ),
  # HEAVY: fits 12 Gaussian GAMMs for each reference-pool allocation.
  tar_target(
    lcbd_model_bundle,
    fit_lcbd_allocation_models(
      allocation_data = lcbd_allocation_bundle$allocation_data,
      allocation_seeds = lcbd_allocation_seeds,
      focal_predictors = lcbd_settings$lcbd_focal_predictors,
      focal_labels = lcbd_settings$lcbd_focal_labels,
      formulas = lcbd_formulas,
      prediction_grids = lcbd_prediction_grids,
      model_threads = lcbd_settings$model_threads
    )
  ),
  tar_target(
    lcbd_ensemble_bundle,
    summarise_lcbd_ensemble(
      per_allocation_results =
        lcbd_model_bundle$per_allocation_results,
      per_allocation_curves =
        lcbd_model_bundle$per_allocation_curves,
      allocation_summary = lcbd_allocation_bundle$allocation_summary,
      expected_model_count = lcbd_settings$lcbd_n_allocations *
        length(lcbd_settings$lcbd_focal_predictors) *
        length(lcbd_formulas)
    )
  ),
  tar_target(
    lcbd_results_file,
    save_rds_quiet(
      assemble_lcbd_results(
        config = lcbd_settings,
        lcbd_inputs = lcbd_inputs,
        distance_validation = lcbd_distance_validation,
        gradient_thresholds = lcbd_gradient_thresholds,
        allocation_bundle = lcbd_allocation_bundle,
        model_bundle = lcbd_model_bundle,
        ensemble_bundle = lcbd_ensemble_bundle
      ),
      lcbd_paths$lcbd_results_path
    ),
    format = "file"
  )
)
