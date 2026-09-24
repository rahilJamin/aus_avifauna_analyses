# =============================================================================
# pipeline_sensitivity.R
# Separate targets for the supplementary robustness analyses.
#
# This pipeline reads completed main checkpoints but never writes to their
# paths. It stores its cache in _targets_sensitivity and exports only below
# outputs/sensitivity.
# =============================================================================

sensitivity_targets <- list(
  # Track main checkpoints and spatial layers explicitly. The sensitivity store
  # must never reuse comparisons after a main result or alignment input changes.
  tar_target(
    sensitivity_main_files,
    c(
      cfg$lcbd_results_path,
      cfg$regional_results_path,
      cfg$chao_table_path,
      cfg$hex_metadata_path,
      cfg$species_matrix_path,
      cfg$candidate_model_df_path,
      cfg$landscape$paths$selected_metrics
    ),
    format = "file"
  ),
  tar_target(
    sensitivity_spatial_files,
    c(
      track_shapefile_files(cfg$hex_2_5_path),
      track_shapefile_files(cfg$ibra_path)
    ),
    format = "file"
  ),
  tar_target(
    sensitivity_occurrence_file,
    cfg$hex_occ_path,
    format = "file"
  ),
  tar_target(
    sensitivity_gllvm_files,
    c(
      cfg$model_df_path,
      cfg$species_matrix_path,
      cfg$trait_path,
      cfg$gllvm_env_path,
      cfg$gllvm_fourth_corner_path
    ),
    format = "file"
  ),

  # Read and validate main checkpoints once. These are dependencies, not copies.
  tar_target(
    sensitivity_baseline,
    read_sensitivity_baseline(cfg, sensitivity_cfg, sensitivity_main_files)
  ),
  tar_target(
    sensitivity_scenarios,
    define_sensitivity_scenarios(cfg, sensitivity_cfg),
    iteration = "list"
  ),
  # The expensive occurrence subset is used only by daily-incidence and
  # structured-protocol branches. The source file itself remains read only.
  tar_target(
    sensitivity_occurrences,
    read_sampling_occurrences(
      sensitivity_occurrence_file,
      sensitivity_cfg
    )
  ),
  tar_target(
    sensitivity_prepared_scenario,
    {
      invisible(sensitivity_spatial_files)
      prepare_sensitivity_scenario(
        sensitivity_scenarios, sensitivity_baseline, sensitivity_occurrences,
        cfg, sensitivity_cfg
      )
    },
    pattern = map(sensitivity_scenarios),
    iteration = "list"
  ),
  # HEAVY: every successful branch, including the main 70% setting, uses the
  # same 30 allocation seeds. Unsupported support is recorded explicitly.
  tar_target(
    sensitivity_lcbd_result,
    fit_sensitivity_lcbd_scenario(
      sensitivity_prepared_scenario, cfg, sensitivity_cfg
    ),
    pattern = map(sensitivity_prepared_scenario),
    iteration = "list"
  ),

  tar_target(
    sensitivity_regional_result,
    fit_sensitivity_regional_scenario(
      sensitivity_prepared_scenario, sensitivity_baseline, cfg
    ),
    pattern = map(sensitivity_prepared_scenario),
    iteration = "list"
  ),

  # Noisy Miner residualisation keeps Y and the same 30 raw-detection pool supports
  # fixed, then refits only the three residualised response curves per allocation.
  tar_target(
    miner_residualisation,
    {
      invisible(sensitivity_spatial_files)
      residualise_noisy_miner(
        sensitivity_baseline$main$analysis_site_data,
        cfg$hex_2_5_path, cfg$ibra_path, cfg, sensitivity_cfg
      )
    }
  ),
  tar_target(
    miner_allocations,
    prepare_miner_allocations(
      miner_residualisation, sensitivity_baseline$main,
      sensitivity_baseline$species_matrix, cfg, sensitivity_cfg
    )
  ),
  # HEAVY: 30 allocations x three LCBD components, with the original formula.
  tar_target(
    miner_lcbd_comparison,
    fit_miner_sensitivity(
      miner_allocations, sensitivity_baseline$main, cfg, sensitivity_cfg
    )
  ),

  # This is a twelve-model allocation-specific check, not a substitute for the
  # main repeated-allocation ensemble.
  tar_target(
    lcbd_spatial_control,
    diagnose_lcbd_models(sensitivity_baseline$main, cfg, sensitivity_cfg)
  ),

  # The one additional large model removes dbMEMs from the completed
  # environmental GLLVM while keeping exactly the same response and base model.
  tar_target(
    sensitivity_gllvm_prepared,
    {
      invisible(sensitivity_gllvm_files)
      invisible(sensitivity_spatial_files)
      prepare_sensitivity_gllvm(cfg)
    }
  ),
  # HEAVY and isolated: child process exits before subsequent diagnostics start.
  tar_target(
    gllvm_without_dbmem_file,
    {
      invisible(sensitivity_gllvm_files)
      sensitivity_worker(
        "fit_gllvm_without_dbmem",
        list(
          prepared = sensitivity_gllvm_prepared,
          main_model_path = cfg$gllvm_env_path,
          config = cfg,
          settings = sensitivity_cfg
        ),
        sensitivity_source_files
      )
    },
    format = "file"
  ),
  tar_target(
    gllvm_diagnostic_definition,
    list(
      list(model_path = cfg$gllvm_env_path, model_id = "environment_with_dbmem"),
      list(model_path = cfg$gllvm_fourth_corner_path, model_id = "fourth_corner"),
      list(model_path = gllvm_without_dbmem_file, model_id = "environment_without_dbmem")
    ),
    iteration = "list"
  ),
  tar_target(
    gllvm_diagnostic_files,
    {
      invisible(sensitivity_gllvm_files)
      sensitivity_worker(
        "diagnose_gllvm_model",
        c(
          gllvm_diagnostic_definition,
          list(prepared = sensitivity_gllvm_prepared, config = cfg,
               settings = sensitivity_cfg)
        ),
        sensitivity_source_files
      )
    },
    pattern = map(gllvm_diagnostic_definition),
    format = "file"
  ),

  # This target writes tables/figures only; all upstream fitted values stay
  # available in the separate targets store and compact summary checkpoint.
  tar_target(
    sensitivity_publication_files,
    produce_sensitivity_outputs(
      sensitivity_baseline, sensitivity_lcbd_result,
      sensitivity_regional_result,
      list(residualisation = miner_residualisation,
           comparison = miner_lcbd_comparison),
      lcbd_spatial_control, gllvm_diagnostic_files, sensitivity_cfg
    ),
    format = "file"
  )
)
