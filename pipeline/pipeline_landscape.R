# =============================================================================
# Script: pipeline/pipeline_landscape.R
# Purpose: Build landscape metrics for the current Chao2-retained hexagons.
# Workflow stage: Landscape metrics.
# Inputs: CLUM and NVIS rasters with category sidecars; shared study extent;
#         generated hexagon grid; identifiers retained by the iNEXT stage.
# Outputs: Three reclassified rasters and the selected landscape-metric table.
# Execution note: The user runs landscape_metrics_file or a downstream target.
# =============================================================================

# Stage settings ----

landscape_paths <- cfg$landscape$paths
landscape_classification <- cfg$landscape$classification
landscape_extraction <- cfg$landscape$extraction
landscape_progress <- cfg$show_progress

# Target definitions ----

landscape_targets <- list(
  tar_target(
    landscape_site_ids,
    as.character(inext_matrix_bundle$valid$hex_id)
  ),
  tar_target(
    landuse_raster_files,
    track_raster_files(landscape_paths$landuse),
    format = "file"
  ),
  tar_target(
    vegetation_raster_files,
    track_raster_files(landscape_paths$vegetation),
    format = "file"
  ),
  tar_target(
    landscape_reclassified_files,
    reclassify_landscape_rasters(
      landuse_path = tracked_file_path(landuse_raster_files, landscape_paths$landuse),
      vegetation_path = tracked_file_path(vegetation_raster_files, landscape_paths$vegetation),
      study_area_path = main_shapefile_path(study_area_files),
      classification = landscape_classification,
      output_paths = landscape_paths[c("vegetation_classes", "vegetation_binary", "landuse_classes")]
    ),
    format = "file",
    packages = c("sf", "terra")
  ),
  # Use the same site eligibility as the species-matrix stage. If occurrence
  # filters or completeness settings change the IDs, extraction updates too.
  tar_target(
    landscape_sample_sites,
    prepare_landscape_sites(
      hex_path = main_shapefile_path(hex_grid_files),
      site_ids = landscape_site_ids,
      raster_crs = terra::crs(terra::rast(tracked_file_path(
        landscape_reclassified_files, landscape_paths$landuse_classes
      )))
    ),
    packages = c("sf", "terra")
  ),
  tar_target(
    landscape_vegetation_metrics,
    extract_landscape_metrics(
      tracked_file_path(landscape_reclassified_files, landscape_paths$vegetation_classes),
      landscape_sample_sites, landscape_extraction$buffer_m,
      landscape_extraction$vegetation_metrics, landscape_progress
    ),
    packages = c("sf", "terra", "landscapemetrics")
  ),
  tar_target(
    landscape_woodland_metrics,
    extract_landscape_metrics(
      tracked_file_path(landscape_reclassified_files, landscape_paths$vegetation_binary),
      landscape_sample_sites, landscape_extraction$buffer_m,
      landscape_extraction$woodland_metrics, landscape_progress
    ),
    packages = c("sf", "terra", "landscapemetrics")
  ),
  tar_target(
    landscape_landuse_metrics,
    extract_landscape_metrics(
      tracked_file_path(landscape_reclassified_files, landscape_paths$landuse_classes),
      landscape_sample_sites, landscape_extraction$buffer_m,
      landscape_extraction$landuse_metrics, landscape_progress
    ),
    packages = c("sf", "terra", "landscapemetrics")
  ),
  tar_target(
    landscape_metric_bundle,
    assemble_landscape_metrics(
      landscape_vegetation_metrics, landscape_woodland_metrics,
      landscape_landuse_metrics, landscape_extraction
    )
  ),
  tar_target(
    landscape_metrics_file,
    save_rds_quiet(landscape_metric_bundle$selected_metrics,
                   landscape_paths$selected_metrics),
    format = "file"
  )
)
