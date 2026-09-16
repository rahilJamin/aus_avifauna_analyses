test_that("configured paths are portable and generated files stay outside raw inputs", {
  source_config()
  source_function("safe_output_writers.R")

  path_values <- unlist(cfg[grepl("(_path|_dir)$", names(cfg))], use.names = TRUE)
  path_values <- path_values[!is.na(path_values)]

  is_absolute <- grepl("^[A-Za-z]:[/\\\\]", path_values) |
    grepl("^/", path_values) |
    grepl("^\\\\\\\\", path_values)
  expect_false(
    any(is_absolute),
    info = paste(names(path_values)[is_absolute], collapse = ", ")
  )

  generated_names <- c(
    "data_dir", "output_dir", "table_dir", "figure_dir", "model_dir",
    "log_dir", "diagnostic_dir", "cleaned_records_path", "hex_occ_path",
    "chao_table_path", "species_matrix_path", "hex_metadata_path",
    "candidate_model_df_path", "model_df_path", "beta_matrices_path",
    "lcbd_results_path", "regional_results_path", "gllvm_metadata_path",
    "gllvm_env_path", "gllvm_fourth_corner_path", "dbmem_selection_path"
  )
  expect_true(all(generated_names %in% names(cfg)))
  expect_false(any(vapply(
    unlist(cfg[generated_names], use.names = FALSE),
    path_is_inside,
    logical(1),
    parent = cfg$raw_data_dir
  )))
})

test_that("configured scientific cutoffs are internally consistent", {
  source_config()

  expect_s3_class(cfg$start_date, "Date")
  expect_s3_class(cfg$end_date, "Date")
  expect_lt(cfg$start_date, cfg$end_date)

  expect_gt(cfg$coord_precision_m, 0)
  expect_gte(cfg$min_sampling_units, 1)
  expect_true(cfg$min_completeness > 0 && cfg$min_completeness <= 1)
  expect_gte(cfg$min_species, 1)
  expect_true(cfg$occupancy_threshold > 0 && cfg$occupancy_threshold <= 1)

  expect_true(cfg$lcbd_lower_tail_probability > 0)
  expect_lt(cfg$lcbd_lower_tail_probability, cfg$lcbd_upper_tail_probability)
  expect_lt(cfg$lcbd_upper_tail_probability, 1)
  expect_true(cfg$lcbd_curve_lower_probability > 0)
  expect_lt(cfg$lcbd_curve_lower_probability, cfg$lcbd_curve_upper_probability)
  expect_lt(cfg$lcbd_curve_upper_probability, 1)
  expect_gte(cfg$lcbd_minimum_pool_sites, 2)
  expect_gte(cfg$lcbd_minimum_tail_sites, 1)
  expect_gte(cfg$lcbd_minimum_informative_pools, 1)
  expect_gte(cfg$lcbd_n_allocations, 2)
  expect_gte(cfg$lcbd_curve_grid_points, 2)
  expect_gt(cfg$lcbd_matrix_validation_tolerance, 0)

  expect_true(cfg$regional_lower_tail_probability > 0)
  expect_lt(
    cfg$regional_lower_tail_probability,
    cfg$regional_upper_tail_probability
  )
  expect_lt(cfg$regional_upper_tail_probability, 1)
  expect_lt(
    cfg$regional_agriculture_low_max,
    cfg$regional_agriculture_high_min
  )
  expect_gte(cfg$regional_bootstrap_draws, 2)

  expect_true(cfg$gllvm_dbmem_global_alpha > 0)
  expect_lte(cfg$gllvm_dbmem_global_alpha, 1)
  expect_true(cfg$gllvm_dbmem_alpha > 0 && cfg$gllvm_dbmem_alpha <= 1)
  expect_gte(cfg$model_threads, 1)
  expect_identical(
    cfg$gllvm_base_predictors,
    c(
      "pland_agriculture_2500m",
      "area_mn_woodland_2500m",
      "cohesion_woodland_2500m",
      "log_effort"
    )
  )

  expect_true(all(vapply(
    cfg[c("input_lonlat_crs", "model_coord_crs", "regional_crs")],
    is.numeric,
    logical(1)
  )))

  expect_named(cfg$lcbd_focal_predictors)
  expect_named(cfg$lcbd_focal_labels)
  expect_identical(names(cfg$lcbd_focal_predictors), names(cfg$lcbd_focal_labels))
  expect_setequal(
    names(cfg$lcbd_focal_predictors),
    c("agriculture", "patch_area", "cohesion", "noisy_miner")
  )
  expect_identical(
    cfg$lcbd_focal_predictors[["noisy_miner"]],
    "miner_detection_rate"
  )
})
