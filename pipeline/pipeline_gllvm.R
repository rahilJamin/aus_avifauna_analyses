# =============================================================================
# pipeline_gllvm.R
# Targets for GLLVM data preparation, spatial selection, and model fitting.
# =============================================================================

# Stage settings ----
gllvm_setting_names <- setdiff(
  grep("^gllvm_", names(cfg), value = TRUE),
  c("gllvm_metadata_path", "gllvm_env_path", "gllvm_fourth_corner_path")
)
gllvm_settings <- cfg[c(
  gllvm_setting_names,
  "occupancy_threshold",
  "noisy_miner_species",
  "model_coord_crs",
  "model_threads"
)]

gllvm_paths <- cfg[c(
  "dbmem_selection_path",
  "gllvm_env_path",
  "gllvm_fourth_corner_path",
  "gllvm_metadata_path"
)]

# Target definitions ----

gllvm_targets <- list(
# 06 GLLVM/JSDM models ----

  tar_target(
    species_traits,
    readr::read_csv(species_traits_file, show_col_types = FALSE)
  ),
  tar_target(
    gllvm_site_data,
    assign_ibra_subregions(
      model_df = model_data_bundle$model_df,
      ibra_path = main_shapefile_path(ibra_files),
      model_coord_crs = gllvm_settings$model_coord_crs
    )
  ),
  tar_target(
    gllvm_prepared_data,
    prepare_gllvm_model_data(
      model_df = gllvm_site_data,
      species_matrix = inext_matrix_bundle$species_matrix,
      traits_raw = species_traits,
      base_predictors = gllvm_settings$gllvm_base_predictors,
      occupancy_threshold = gllvm_settings$occupancy_threshold,
      noisy_miner_species = gllvm_settings$noisy_miner_species
    )
  ),
  tar_target(
    gllvm_spatial_selection,
    select_gllvm_dbmem(
      species_matrix = gllvm_prepared_data$Y,
      coordinates = gllvm_prepared_data$gllvm_df[c("x", "y")],
      seed = gllvm_settings$gllvm_dbmem_seed,
      permutations = gllvm_settings$gllvm_dbmem_permutations,
      global_alpha = gllvm_settings$gllvm_dbmem_global_alpha,
      selection_alpha = gllvm_settings$gllvm_dbmem_alpha
    )
  ),
  tar_target(
    gllvm_dbmem_file,
    save_rds_quiet(
      gllvm_spatial_selection[c(
        "global_test",
        "adjusted_r2",
        "forward_selection",
        "selected_mem_names"
      )],
      gllvm_paths$dbmem_selection_path
    ),
    format = "file"
  ),
  tar_target(
    gllvm_predictor_data,
    dplyr::bind_cols(
      gllvm_prepared_data$X,
      gllvm_spatial_selection$selected_mems
    )
  ),
  tar_target(
    gllvm_formulas,
    define_gllvm_formulas(gllvm_spatial_selection$selected_mem_names)
  ),
  # HEAVY: environment-only GLLVM checkpoint.
  tar_target(
    gllvm_environment_model_file,
    {
      TMB::openmp(n = gllvm_settings$model_threads)
      fitted_model <- gllvm::gllvm(
        y = as.matrix(gllvm_prepared_data$Y),
        X = gllvm_predictor_data,
        formula = gllvm_formulas$environment,
        family = stats::binomial(link = "probit"),
        num.lv = gllvm_settings$gllvm_env_num_lv,
        row.eff = ~ (1 | ibra_sub),
        studyDesign = data.frame(
          ibra_sub = gllvm_prepared_data$gllvm_df$ibra_sub
        ),
        method = "VA",
        starting.val = "res",
        n.init = gllvm_settings$gllvm_n_init,
        seed = gllvm_settings$gllvm_seed,
        trace = TRUE
      )
      save_rds_quiet(fitted_model, gllvm_paths$gllvm_env_path)
    },
    format = "file"
  ),
  # HEAVY: fourth-corner GLLVM checkpoint.
  tar_target(
    gllvm_fourth_corner_model_file,
    {
      TMB::openmp(n = gllvm_settings$model_threads)
      fitted_model <- gllvm::gllvm(
        y = as.matrix(gllvm_prepared_data$Y),
        X = gllvm_predictor_data,
        TR = gllvm_prepared_data$trait_matrix,
        formula = gllvm_formulas$fourth_corner,
        family = stats::binomial(link = "probit"),
        num.lv = gllvm_settings$gllvm_fourth_corner_num_lv,
        row.eff = ~ (1 | ibra_sub),
        studyDesign = data.frame(
          ibra_sub = gllvm_prepared_data$gllvm_df$ibra_sub
        ),
        method = "EVA",
        n.init = gllvm_settings$gllvm_n_init,
        seed = gllvm_settings$gllvm_seed,
        trace = TRUE
      )
      save_rds_quiet(fitted_model, gllvm_paths$gllvm_fourth_corner_path)
    },
    format = "file"
  ),
  tar_target(
    gllvm_metadata_file,
    {
      # Keep the metadata checkpoint downstream of both fitted-model files.
      invisible(gllvm_environment_model_file)
      invisible(gllvm_fourth_corner_model_file)

      metadata <- list(
        dimensions = dplyr::mutate(
          gllvm_prepared_data$dimensions,
          n_environment_predictors = ncol(gllvm_predictor_data),
          n_selected_dbmem =
            length(gllvm_spatial_selection$selected_mem_names)
        ),
        species_counts = gllvm_prepared_data$species_counts,
        traits = gllvm_prepared_data$traits,
        species_order = colnames(gllvm_prepared_data$Y),
        missing_traits = gllvm_prepared_data$missing_traits,
        minimum_required_presences =
          gllvm_prepared_data$minimum_required_presences,
        occupancy_threshold = gllvm_settings$occupancy_threshold
      )

      save_rds_quiet(metadata, gllvm_paths$gllvm_metadata_path)
    },
    format = "file"
  )
)
