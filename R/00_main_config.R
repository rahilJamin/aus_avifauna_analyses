# =============================================================================
# 00_main_config.R
# Central configuration for the 2020-onward iNEXT main analysis.
#
# Sourcing this file only defines `cfg`; it does not create folders, read data,
# fit models, or write outputs. The run scripts read these values so that paths,
# filters, thresholds, seeds, and labels are kept in one place.
# =============================================================================

cfg <- local({
  # Project roots ----
  # Keep these relative to the repository root so the workflow can be run on
  # another computer without editing machine-specific absolute paths.
  workflow_dir <- "."
  script_dir <- "R"

  # Input and output folders ----
  # `raw_data_dir` stores required external inputs. Scripts should read from it
  # but never write to it. `data_dir` stores reproducible derived data created by
  # the workflow before model fitting.
  raw_data_dir <- file.path("inputs", "data_raw")
  data_dir <- file.path("inputs", "data_processed")
  output_dir <- "outputs"

  # Generated outputs are separated by purpose so tables, figures, fitted model
  # objects, diagnostics, and logs can be cleaned or shared independently.
  table_dir <- file.path(output_dir, "tables")
  figure_dir <- file.path(output_dir, "figures")
  model_dir <- file.path(output_dir, "models")
  log_dir <- file.path("docs", "logs")
  diagnostic_dir <- file.path(output_dir, "diagnostics")

  list(
    # Analysis identity ----
    # Added to selected tables/logs so outputs can be traced to this workflow.
    analysis_id = "main_2020_onward_inext",

    # Directory paths ----
    workflow_dir = workflow_dir,
    script_dir = script_dir,
    raw_data_dir = raw_data_dir,
    data_dir = data_dir,
    output_dir = output_dir,
    table_dir = table_dir,
    figure_dir = figure_dir,
    model_dir = model_dir,
    log_dir = log_dir,
    diagnostic_dir = diagnostic_dir,

    # Required external inputs ----
    # These files are treated as fixed inputs. They should exist before running
    # the main workflow and should not be overwritten by analysis scripts.
    raw_ala_path = file.path(raw_data_dir, "species", "ala_df_raw.csv"),
    hex_2_5_path = file.path(raw_data_dir, "spatial_grids", "hex_2_5.shp"),
    selected_metrics_2_5_path = file.path(
      raw_data_dir,
      "lsm_processed",
      "selected_metrics_2_5km.rds"
    ),
    trait_path = file.path(raw_data_dir, "species", "trait_set.csv"),
    ibra_path = file.path(
      raw_data_dir,
      "spatial_grids",
      "IBRA7_subregions",
      "IBRA7_subregions.shp"
    ),

    # Derived data from preprocessing stages ----
    # 01_clean_records.R writes `cleaned_records_path`.
    cleaned_records_path = file.path(
      data_dir,
      "ala_cleaned_2020_inext.csv"
    ),
    # 02_assign_hexes.R writes `hex_occ_path`.
    hex_occ_path = file.path(
      data_dir,
      "hex_occ_2_5_2020_inext.csv"
    ),
    # 03_build_inext_matrix.R writes Chao2 estimates, retained hex metadata,
    # and the retained binary species matrix.
    chao_table_path = file.path(
      data_dir,
      "chao2_inext_2020.csv"
    ),
    species_matrix_path = file.path(
      data_dir,
      "species_matrix_2_5km_2020_inext.rds"
    ),
    hex_metadata_path = file.path(
      data_dir,
      "hex_metadata_2_5km_2020_inext.csv"
    ),

    # Derived data for LCBD and downstream model stages ----
    # 03b_prepare_model_data.R writes these objects. `candidate_model_df_path`
    # is the site table before reference-pool allocation; `model_df_path` is the
    # aligned modelling table for later analyses; `beta_matrices_path` stores the
    # Podani-Jaccard distance matrices used by the LCBD workflow.
    candidate_model_df_path = file.path(
      data_dir,
      "model_df_2_5km_2020_inext_candidate_pools.rds"
    ),
    model_df_path = file.path(
      data_dir,
      "model_df_2_5km_2020_inext.rds"
    ),
    beta_matrices_path = file.path(
      data_dir,
      "beta_matrices_2020_inext.rds"
    ),

    # Main model result checkpoints ----
    # These are generated model objects, not manual inputs. The user runs the
    # corresponding analysis scripts to create them.
    lcbd_results_path = file.path(
      model_dir,
      "04_lcbd_results.rds"
    ),
    regional_results_path = file.path(
      model_dir,
      "05_regional_convergence_results.rds"
    ),
    gllvm_metadata_path = file.path(
      model_dir,
      "06_gllvm_metadata.rds"
    ),
    gllvm_env_path = file.path(
      model_dir,
      "gllvm_env_spatial_2020_inext.rds"
    ),
    gllvm_fourth_corner_path = file.path(
      model_dir,
      "gllvm_fourth_corner_2020_inext.rds"
    ),
    dbmem_selection_path = file.path(
      model_dir,
      "dbmem_forward_selection_2020_inext.rds"
    ),

    # Coordinate reference systems ----
    # Used in stages 02, 03b, 04a, and spatial diagnostics. Keep CRS choices here
    # so all spatial joins, coordinates, and distance-based pool allocations are
    # consistent across scripts.
    input_lonlat_crs = 4326,      # WGS84 lon/lat for raw occurrence coordinates.
    model_coord_crs = 32655,      # Projected CRS used for model x/y coordinates.

    # Record-cleaning settings: 01_clean_records.R ----
    # Retain only records inside this date window.
    start_date = as.Date("2020-01-01"),
    end_date = as.Date("2026-02-28"),

    # Protocol filter passed to the occurrence-cleaning function. The main
    # analysis keeps all retained protocols after the protocol screen.
    protocol_filter = "all_retained",

    # Maximum accepted coordinate uncertainty in metres. Records with missing
    # uncertainty can still be kept when their protocol gives enough spatial
    # structure, as defined in clean_occurrence_records().
    coord_precision_m = 2500,

    # Sampling-completeness settings: 03_build_inext_matrix.R ----
    # Minimum number of incidence units needed before a hexagon can be evaluated.
    min_sampling_units = 6,

    # Minimum estimated Chao2 completeness required for a hexagon to be retained.
    min_completeness = 0.70,

    # Minimum observed richness required for a hexagon to enter the species matrix.
    min_species = 5,

    # Checklist/incidence definition used by iNEXT and Noisy Miner detection.
    incidence_unit = "month_year",

    # LCBD species-matrix settings: 03b_prepare_model_data.R ----
    # Noisy Miner is removed from the LCBD response matrix so its detection rate
    # can be used as an explanatory gradient rather than part of the response.
    noisy_miner_species = "Manorina melanocephala",

    # Remove strongly urban hexagons before LCBD and GLLVM model preparation.
    # This threshold applies to the processed landscape metric percentage.
    urban_cover_max = 10,

    # LCBD reference-pool settings: 04_fit_lcbd.R ----
    # Australian Albers is used intentionally for continent-wide spatial
    # clustering of reference pools.
    lcbd_reference_pool_crs = 3577,

    # Number of spatial reference pools created by k-means for each allocation.
    lcbd_reference_pool_k = 25L,

    # Number of repeated reference-pool allocations in the ensemble.
    lcbd_n_allocations = 100L,

    # First random seed for allocation-specific k-means runs. Allocation seeds are
    # this value plus 0, 1, 2, ... so every allocation can be reproduced.
    lcbd_allocation_seed_start = 42L,

    # Number of random starts used inside k-means for each allocation.
    lcbd_reference_pool_nstart = 100L,

    # Minimum number of sites required for a reference pool to be included in the
    # LCBD calculation.
    lcbd_minimum_pool_sites = 30L,

    # Minimum number of sites required in both tails of a focal gradient for a
    # pool to be informative for that gradient.
    lcbd_minimum_tail_sites = 5L,

    # Minimum number of informative pools required before fitting a GAMM for a
    # focal gradient in an allocation.
    lcbd_minimum_informative_pools = 5L,

    # Global gradient quantiles used to define low and high tails when checking
    # whether each reference pool covers enough of a focal gradient.
    lcbd_lower_tail_probability = 0.25,
    lcbd_upper_tail_probability = 0.75,

    # Central quantile range used for prediction curves. The final shared curve
    # range is the overlap of these ranges across all allocations.
    lcbd_curve_lower_probability = 0.05,
    lcbd_curve_upper_probability = 0.95,

    # Number of x-values used in each fitted LCBD response curve.
    lcbd_curve_grid_points = 100L,

    # Numerical tolerances for LCBD diagnostics. `lcbd_negative_tolerance`
    # defines which negative component LCBD values are counted and reported;
    # these do not stop the workflow. `lcbd_sum_tolerance` defines the accepted
    # deviation from a within-pool LCBD sum of one and is enforced.
    lcbd_negative_tolerance = 1e-10,
    lcbd_sum_tolerance = 1e-6,

    # Maximum accepted element-wise difference when stage 04 rebuilds the three
    # beta-diversity matrices saved during model-data preparation.
    lcbd_matrix_validation_tolerance = 1e-12,

    # Focal LCBD gradients. Names are stable identifiers used in tables; values
    # are columns in the prepared model dataframe.
    lcbd_focal_predictors = c(
      agriculture = "pland_agriculture_2500m",
      patch_area = "area_mn_woodland_2500m",
      cohesion = "cohesion_woodland_2500m",
      noisy_miner = "miner_detection_rate"
    ),

    # Human-readable labels for LCBD tables and plots.
    lcbd_focal_labels = c(
      agriculture = "Agricultural cover",
      patch_area = "Mean woodland patch area",
      cohesion = "Woodland cohesion",
      noisy_miner = "Noisy Miner detection rate"
    ),


    # Regional convergence settings: 05_regional_convergence.R ----
    # Agriculture uses fixed percentage endpoints for intact and disturbed
    # habitat comparisons.
    regional_agriculture_low_max = 20,
    regional_agriculture_high_min = 80,

    # Other gradients use study-wide lower and upper tails for endpoint sites.
    regional_lower_tail_probability = 0.20,
    regional_upper_tail_probability = 0.80,

    # Retain an IBRA subregion for a gradient only when it has enough sites at
    # both low and high endpoints; retain a gradient only when enough subregions
    # pass that support rule.
    regional_min_sites_per_endpoint = 3L,
    regional_min_retained_subregions = 4L,

    # Bootstrap sites within each retained subregion and endpoint. This is 99 by
    # default for practical manual reruns; increase for final publication checks
    # if needed.
    regional_bootstrap_draws = 99L,
    regional_bootstrap_seed = 4202L,

    # Projected CRS used for regional spatial joins against IBRA subregions.
    regional_crs = 3577,
    # GLLVM/JSDM settings: main GLLVM analysis ----
    # Species must occur in at least this proportion of retained sites before
    # entering the GLLVM/JSDM stage.
    occupancy_threshold = 0.02,

    # Environmental columns supplied to both GLLVM models. Spatial dbMEM
    # variables selected later are appended to this fixed set.
    gllvm_base_predictors = c(
      "pland_agriculture_2500m",
      "area_mn_woodland_2500m",
      "cohesion_woodland_2500m",
      "log_effort"
    ),

    # Seed for reproducible GLLVM/JSDM model fitting.
    gllvm_seed = 42,

    # Seed used before dbMEM permutation tests and forward selection.
    gllvm_dbmem_seed = 4206L,

    # Number of permutations used by the dbMEM global test and forward selection.
    gllvm_dbmem_permutations = 199L,

    # Significance level for deciding whether the global dbMEM test supports
    # running forward selection.
    gllvm_dbmem_global_alpha = 0.05,

    # Forward-selection alpha for retaining dbMEM spatial eigenvectors.
    gllvm_dbmem_alpha = 0.01,

    # Number of latent variables in the environment-only GLLVM.
    gllvm_env_num_lv = 1L,

    # Number of latent variables in the fourth-corner GLLVM.
    gllvm_fourth_corner_num_lv = 0L,

    # Number of starting values used by each GLLVM fit.
    gllvm_n_init = 3L,

    # Number of CPU threads used by model-fitting scripts. Keep two cores free
    # so the computer remains responsive during long manual runs.
    model_threads = max(1L, parallel::detectCores() - 2L),

    # Console progress messages/progress bars for long data-preparation steps.
    show_progress = TRUE
  )
})


