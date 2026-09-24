#' Read completed main checkpoints for sensitivity comparisons
#' @param config Main configuration (read only).
#' @param settings Sensitivity configuration.
#' @param tracked_main_files Main checkpoints tracked as file targets so changed
#'   publication results invalidate the separate sensitivity cache.
#' @return Validated main objects and inputs; no GLLVM model is loaded here.
read_sensitivity_baseline <- function(config, settings, tracked_main_files) {
  # `tracked_main_files` is supplied by the separate targets pipeline. Keeping
  # it as an explicit argument makes a changed main checkpoint invalidate this
  # sensitivity baseline instead of silently reusing a stale cached object.
  invisible(tracked_main_files)
  main <- readRDS(config$lcbd_results_path)
  seeds <- main$settings$allocation_seeds
  sensitivity_seed_count <- length(settings$allocation_seeds)
  if (length(seeds) < sensitivity_seed_count ||
      !identical(as.numeric(seeds[seq_len(sensitivity_seed_count)]),
                 as.numeric(settings$allocation_seeds)) ||
      main$settings$reference_pool_k != config$lcbd_reference_pool_k ||
      !identical(main$settings$model_family,
                 "Gaussian identity for all three LCBD components")) {
    stop("The main LCBD checkpoint cannot supply the requested sensitivity seeds and settings.",
         call. = FALSE)
  }
  if (length(unique(main$model_results_by_allocation$allocation)) != length(seeds) ||
      nrow(main$model_results_by_allocation) != length(seeds) * 12L) {
    stop("The main checkpoint does not contain the complete four-gradient LCBD ensemble.", call. = FALSE)
  }
  chao <- readr::read_csv(config$chao_table_path, show_col_types = FALSE)
  metadata <- readr::read_csv(config$hex_metadata_path, show_col_types = FALSE)
  matrix <- readRDS(config$species_matrix_path)
  retained <- retain_sampled_hexagons(chao, config$min_sampling_units,
                                    config$min_species, config$min_completeness)
  if (!setequal(as.character(retained$hex_id), rownames(matrix)) ||
      !setequal(as.character(metadata$hex_id), rownames(matrix))) {
    stop("The main Chao2, metadata and species checkpoints describe different retained sites.", call. = FALSE)
  }
  candidate <- readRDS(config$candidate_model_df_path)
  shared_columns <- names(candidate)[names(candidate) %in% names(main$analysis_site_data)]
  reference <- main$analysis_site_data
  aligned <- candidate[match(reference$hex_id, candidate$hex_id), shared_columns, drop = FALSE]
  if (!setequal(candidate$hex_id, reference$hex_id) ||
      !isTRUE(all.equal(as.data.frame(aligned), as.data.frame(reference[shared_columns]),
                        check.attributes = FALSE, tolerance = 1e-12))) {
    stop("The LCBD fit and prepared candidate data are inconsistent; finish the main run first.", call. = FALSE)
  }
  regional <- readRDS(config$regional_results_path)
  if (!is.list(regional) || !all(c("settings", "endpoint_contrasts") %in% names(regional)) ||
      !nrow(regional$endpoint_contrasts)) {
    stop("The completed regional-convergence checkpoint is missing usable endpoint contrasts.",
         call. = FALSE)
  }
  list(main = main, chao = chao, metadata = metadata, species_matrix = matrix,
       landscape = readRDS(config$landscape$paths$selected_metrics),
       regional = regional)
}

#' Prepare one completeness alternative using the main scientific functions
#' @param threshold Completeness cutoff, between the main cutoff and one.
#' @param baseline Validated main checkpoint bundle.
#' @param config Main settings.
#' @return Candidate model data, aligned LCBD response and beta matrices.
prepare_completeness_data <- function(threshold, baseline, config) {
  ids <- sensitivity_site_ids(baseline$chao, baseline$species_matrix, threshold, config)
  # Landscape coverage follows Chao2 retention in the main pipeline. Subsetting
  # BEFORE preparation reproduces its scaling population for each alternative.
  landscape <- baseline$landscape[baseline$landscape$hex_id %in% ids, , drop = FALSE]
  metadata <- baseline$metadata[match(ids, baseline$metadata$hex_id), , drop = FALSE]
  prepared <- prepare_lcbd_model_data(
    species_matrix = baseline$species_matrix[ids, , drop = FALSE],
    hex_metadata = metadata, landscape_data = landscape,
    hex_path = config$hex_2_5_path, noisy_miner_species = config$noisy_miner_species,
    urban_cover_max = config$urban_cover_max, model_coord_crs = config$model_coord_crs
  )
  dat <- add_lcbd_pool_coordinates(
    prepared$candidate_model_df, config$hex_2_5_path, config$lcbd_reference_pool_crs
  )
  list(threshold = threshold, site_data = dat,
       species_matrix = prepared$lcbd_species_matrix,
       distance_matrices = prepared$beta_matrices$dist_mats,
       attrition = prepared$attrition, n_chao_sites = length(ids))
}

#' Add reference-pool coordinates to any prepared LCBD site table
#' @param site_data Prepared candidate table with `hex_id`.
#' @param hex_path Generated analysis grid.
#' @param crs Australian Albers CRS used for the main LCBD pool allocation.
#' @return The site table with finite `x_km` and `y_km` columns.
add_lcbd_pool_coordinates <- function(site_data, hex_path, crs) {
  assert_unique_ids(site_data, "hex_id", "sensitivity LCBD site table")
  grid <- sf::st_read(hex_path, quiet = TRUE)
  grid <- grid[as.character(grid$hex_id) %in% as.character(site_data$hex_id), "hex_id"]
  grid <- sf::st_transform(grid, crs)
  points <- suppressWarnings(sf::st_point_on_surface(grid))
  xy <- sf::st_coordinates(points) / 1000
  lookup <- data.frame(hex_id = as.character(points$hex_id),
                       x_km = xy[, 1], y_km = xy[, 2])
  assert_unique_ids(lookup, "hex_id", "sensitivity LCBD coordinate lookup")
  site_data$hex_id <- as.character(site_data$hex_id)
  site_data <- dplyr::left_join(site_data, lookup, by = "hex_id")
  if (anyNA(site_data[c("x_km", "y_km")])) {
    stop("A sensitivity LCBD site is missing a reference-pool coordinate.", call. = FALSE)
  }
  site_data
}
