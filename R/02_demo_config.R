# Demo-only paths and reduced computational settings.
# Source R/00_main_config.R first. This replaces cfg only for _targets_demo.R.

demo_data_dir <- file.path("demo", "data")
demo_input_dir <- file.path(demo_data_dir, "inputs")
demo_output_dir <- file.path("demo", "outputs")
demo_derived_dir <- file.path(demo_output_dir, "derived_data")

demo_settings <- list(
  store = "_targets_demo",
  data_dir = demo_data_dir,
  input_dir = demo_input_dir,
  output_dir = demo_output_dir,
  extent_lonlat = c(xmin = 145, xmax = 149, ymin = -38, ymax = -34.5),
  n_sites = 150L,
  n_species = 50L,
  raster_buffer_m = 5000,
  occurrence_chunk_size = 250000L,
  occupancy_threshold = 0.02
)

demo_lcbd_settings <- list(
  n_allocations = 3L,
  reference_pool_k = 3L,
  minimum_pool_sites = 10L,
  minimum_tail_sites = 1L,
  minimum_informative_pools = 2L
)

demo_regional_settings <- list(
  minimum_sites_per_endpoint = 1L,
  minimum_retained_subregions = 2L
)

cfg$analysis_id <- "demo_2020_onward"
cfg$raw_data_dir <- demo_input_dir
cfg$data_dir <- demo_derived_dir
cfg$output_dir <- demo_output_dir
cfg$table_dir <- file.path(demo_output_dir, "tables")
cfg$figure_dir <- file.path(demo_output_dir, "figures")
cfg$model_dir <- file.path(demo_output_dir, "models")
cfg$diagnostic_dir <- file.path(cfg$figure_dir, "diagnostics")
cfg$log_dir <- file.path(demo_output_dir, "logs")

cfg$raw_ala_path <- file.path(
  demo_input_dir, "species", "ala_df_raw.csv.gz"
)
cfg$study_area_path <- file.path(demo_input_dir, "spatial_grids", "study_area.shp")
cfg$trait_path <- file.path(demo_input_dir, "species", "trait_set.csv")
cfg$ibra_path <- file.path(
  demo_input_dir, "spatial_grids", "ibra7_subregions.shp"
)

cfg$hex_2_5_path <- file.path(demo_derived_dir, "spatial_grids", "hex_2_5.shp")
cfg$cleaned_records_path <- file.path(demo_derived_dir, "ala_cleaned_demo.csv")
cfg$hex_occ_path <- file.path(demo_derived_dir, "hex_occ_demo.csv")
cfg$chao_table_path <- file.path(demo_derived_dir, "chao2_demo.csv")
cfg$species_matrix_path <- file.path(demo_derived_dir, "species_matrix_demo.rds")
cfg$hex_metadata_path <- file.path(demo_derived_dir, "hex_metadata_demo.csv")
cfg$candidate_model_df_path <- file.path(demo_derived_dir, "candidate_model_df_demo.rds")
cfg$model_df_path <- file.path(demo_derived_dir, "model_df_demo.rds")
cfg$beta_matrices_path <- file.path(demo_derived_dir, "beta_matrices_demo.rds")

cfg$lcbd_results_path <- file.path(cfg$model_dir, "lcbd_demo_results.rds")
cfg$regional_results_path <- file.path(cfg$model_dir, "regional_demo_results.rds")
cfg$gllvm_metadata_path <- file.path(cfg$model_dir, "gllvm_demo_metadata.rds")
cfg$gllvm_env_path <- file.path(cfg$model_dir, "gllvm_demo_environment.rds")
cfg$gllvm_fourth_corner_path <- file.path(cfg$model_dir, "gllvm_demo_fourth_corner.rds")
cfg$dbmem_selection_path <- file.path(cfg$model_dir, "gllvm_demo_dbmem_selection.rds")

cfg$landscape$paths$landuse <- file.path(
  demo_input_dir, "raw_rasters", "clum_50m_2023_v2.tif"
)
cfg$landscape$paths$vegetation <- file.path(
  demo_input_dir, "raw_rasters", "nvis_mvg.tif"
)
cfg$landscape$paths$vegetation_classes <- file.path(
  demo_derived_dir, "rasters", "veg_reclass.tif"
)
cfg$landscape$paths$vegetation_binary <- file.path(
  demo_derived_dir, "rasters", "veg_reclass_binary.tif"
)
cfg$landscape$paths$landuse_classes <- file.path(
  demo_derived_dir, "rasters", "landuse_reclass.tif"
)
cfg$landscape$paths$selected_metrics <- file.path(
  demo_derived_dir, "lsm_processed", "selected_metrics_2_5km.rds"
)

# Preserve the main sampling and model filters/formulas. The demo-only LCBD
# ensemble uses fewer spatial pools and allocations to keep the tutorial short.
cfg$occupancy_threshold <- demo_settings$occupancy_threshold
cfg$lcbd_n_allocations <- demo_lcbd_settings$n_allocations
cfg$lcbd_reference_pool_k <- demo_lcbd_settings$reference_pool_k
cfg$lcbd_minimum_pool_sites <- demo_lcbd_settings$minimum_pool_sites
cfg$lcbd_minimum_tail_sites <- demo_lcbd_settings$minimum_tail_sites
cfg$lcbd_minimum_informative_pools <-
  demo_lcbd_settings$minimum_informative_pools
cfg$regional_min_sites_per_endpoint <-
  demo_regional_settings$minimum_sites_per_endpoint
cfg$regional_min_retained_subregions <-
  demo_regional_settings$minimum_retained_subregions
