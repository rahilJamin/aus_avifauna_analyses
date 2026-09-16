# =============================================================================
# pipeline_regional_convergence.R
# Targets for beta-diversity convergence among IBRA7 subregions.
# =============================================================================

# Stage settings ----
regional_setting_names <- grep("^regional_", names(cfg), value = TRUE)
regional_settings <- cfg[
  setdiff(regional_setting_names, "regional_results_path")
]
regional_settings$noisy_miner_species <- cfg$noisy_miner_species

regional_paths <- cfg[c(
  "candidate_model_df_path",
  "species_matrix_path",
  "regional_results_path"
)]

# Target definitions ----

regional_convergence_targets <- list(
# 05 Regional convergence ----

  tar_target(
    regional_inputs,
    load_regional_convergence_inputs(
      site_data_path = tracked_file_path(
        model_data_files,
        regional_paths$candidate_model_df_path
      ),
      landscape_data_path = landscape_metrics_file,
      species_matrix_path = tracked_file_path(
        inext_matrix_files,
        regional_paths$species_matrix_path
      ),
      hex_path = main_shapefile_path(hex_grid_files),
      ibra_path = main_shapefile_path(ibra_files),
      noisy_miner_species = regional_settings$noisy_miner_species,
      regional_crs = regional_settings$regional_crs
    )
  ),
  tar_target(
    regional_gradient_definitions,
    define_regional_gradients(
      site_data = regional_inputs$site_data,
      agriculture_low_maximum =
        regional_settings$regional_agriculture_low_max,
      agriculture_high_minimum =
        regional_settings$regional_agriculture_high_min,
      lower_tail_probability =
        regional_settings$regional_lower_tail_probability,
      upper_tail_probability =
        regional_settings$regional_upper_tail_probability
    )
  ),
  tar_target(
    regional_jaccard_components,
    calculate_regional_jaccard_components(
      species_matrix = regional_inputs$species_matrix
    )
  ),
  tar_target(
    regional_convergence_results,
    analyse_regional_convergence(
      site_data = regional_inputs$site_data,
      species_matrix = regional_inputs$species_matrix,
      gradient_definitions = regional_gradient_definitions,
      component_matrices =
        regional_jaccard_components$component_matrices,
      component_sum_error =
        regional_jaccard_components$component_sum_error,
      minimum_sites_per_endpoint =
        regional_settings$regional_min_sites_per_endpoint,
      minimum_retained_subregions =
        regional_settings$regional_min_retained_subregions,
      bootstrap_draws = regional_settings$regional_bootstrap_draws,
      bootstrap_seed = regional_settings$regional_bootstrap_seed,
      noisy_miner_species = regional_settings$noisy_miner_species,
      species_matrix_source = tracked_file_path(
        inext_matrix_files,
        regional_paths$species_matrix_path
      )
    )
  ),
  tar_target(
    regional_results_file,
    save_rds_quiet(
      regional_convergence_results,
      regional_paths$regional_results_path
    ),
    format = "file"
  )
)
