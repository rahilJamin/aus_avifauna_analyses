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
    study_area_path = file.path(raw_data_dir, "spatial_grids", "study_area.shp"),
    trait_path = file.path(raw_data_dir, "species", "trait_set.csv"),
    ibra_path = file.path(
      raw_data_dir,
      "spatial_grids",
      "IBRA7_subregions",
      "ibra7_subregions.shp"
    ),

    # Derived data from preprocessing stages ----
    # Created from study_area_path; every spatial stage uses this generated grid.
    hex_2_5_path = file.path(data_dir, "spatial_grids", "hex_2_5.shp"),
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

    # Hexagon grid ----
    # WORKFLOW CHANGE: construct the grid in Australian Albers, using metres.
    # The former supplied grid used an approximate kilometres-to-degrees size.
    # The generated cells therefore have new geometries and site identifiers.
    hex_grid = list(
      crs = 3577,
      # Opposite-edge spacing, not radius or area. A full cell has area
      # sqrt(3) / 2 * cellsize_m^2 (about 5.413 square kilometres at 2500 m).
      cellsize_m = 2500,
      flat_topped = FALSE
    ),

    # Landscape metrics ----
    # WORKFLOW CHANGE: sample the hexagons retained by this same 2020-onward
    # iNEXT/Chao2 workflow. No externally supplied valid-site list is needed.
    # Raster classification remains independent of occurrence processing;
    # extraction depends on the generated grid and retained Chao2 identifiers.
    landscape = list(
      paths = list(
        landuse = file.path(raw_data_dir, "raw_rasters", "clum_50m_2023_v2.tif"),
        vegetation = file.path(raw_data_dir, "raw_rasters", "nvis_mvg.tif"),
        vegetation_classes = file.path(data_dir, "rasters", "veg_reclass.tif"),
        vegetation_binary = file.path(data_dir, "rasters", "veg_reclass_binary.tif"),
        landuse_classes = file.path(data_dir, "rasters", "landuse_reclass.tif"),
        selected_metrics = file.path(
          data_dir, "lsm_processed", "selected_metrics_2_5km.rds"
        )
      ),
      classification = list(
        # Common projected CRS for all derived categorical rasters. Each source
        # raster keeps the output grid chosen by terra::project(); no additional
        # common resolution or alignment template is imposed.
        crs = "+proj=utm +zone=55 +datum=WGS84 +units=m +no_defs",
        resampling = "near",
        # CLUM SIMPN groups -> intact, native production, plantation,
        # agriculture, urban. Water (19) and other unmapped categories are NA.
        landuse_groups = c(
          "1" = 1L, "2" = 1L, "3" = 1L, "4" = 2L, "5" = 2L,
          "6" = 3L, "7" = 4L, "8" = 4L, "9" = 4L, "10" = 4L,
          "11" = 4L, "12" = 4L, "13" = 4L, "14" = 5L, "15" = 4L,
          "16" = 5L, "17" = 5L, "18" = 5L, "19" = NA_integer_
        ),
        # NVIS Value codes -> wet/dense forest, sclerophyll woodland,
        # mallee/arid woodland, open/cleared. Water and sea (24, 28) are NA.
        vegetation_groups = c(
          "1" = 1L, "2" = 1L, "15" = 1L,
          "3" = 2L, "4" = 2L, "5" = 2L, "7" = 2L, "8" = 2L,
          "9" = 2L, "10" = 2L, "11" = 2L, "12" = 2L, "29" = 2L,
          "6" = 3L, "13" = 3L, "14" = 3L, "31" = 3L, "32" = 3L,
          "16" = 4L, "17" = 4L, "18" = 4L, "19" = 4L, "20" = 4L,
          "21" = 4L, "22" = 4L, "23" = 4L, "25" = 4L, "26" = 4L,
          "27" = 4L, "99" = 4L, "24" = NA_integer_, "28" = NA_integer_
        ),
        # Fallbacks apply to numeric values absent from the category table. A
        # listed category with no group mapping remains NA. The reconstruction
        # explicitly restores the source vegetation NA mask after classification
        # because current terra versions otherwise assign those cells class 5.
        landuse_other = NA_integer_,
        vegetation_other = 5L,
        # Woodland classes 1--3 become 1; nonhabitat classes 4--5 become 2.
        vegetation_binary = c("1" = 1L, "2" = 1L, "3" = 1L, "4" = 2L, "5" = 2L)
      ),
      extraction = list(
        # Circular sampling windows centred on projected hexagon centroids.
        # `size` in sample_lsm() is the circle radius, in projected metres.
        buffer_m = 2500,
        vegetation_metrics = c("lsm_c_pland", "lsm_c_area_mn"),
        woodland_metrics = c("lsm_c_ed", "lsm_c_cohesion", "lsm_c_area_mn"),
        landuse_metrics = "lsm_c_pland",
        vegetation_labels = c(
          "1" = "wet_dense", "2" = "sclerophyll", "3" = "mallee_arid",
          "4" = "nonhabitat", "5" = "nonhabitat_other"
        ),
        landuse_labels = c(
          "1" = "intact", "2" = "native_prod", "3" = "plantation",
          "4" = "agriculture", "5" = "urban"
        ),
        vegetation_drop_classes = c(4L, 5L),
        woodland_class = 1L,
        # Apply one ten-column selection and complete-case filter, including
        # predictors not subsequently used in the final model formulae.
        selected_metrics = c(
          "pland_wet_dense", "pland_sclerophyll", "pland_mallee_arid",
          "PLAND_habitat", "area_mn_woodland", "cohesion_woodland",
          "pland_native_prod", "pland_plantation", "pland_agriculture", "pland_urban"
        )
      )
    ),

    # Coordinate reference systems ----
    # The grid stays in hex_grid$crs. Occurrence coordinates are transformed
    # from input_lonlat_crs to that grid before joining. Landscape sampling and
    # model preparation explicitly transform the same grid to their own metric
    # CRSs; a CRS is never assigned to coordinates to simulate a transformation.
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


