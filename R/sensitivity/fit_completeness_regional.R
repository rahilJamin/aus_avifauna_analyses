#' Prepare a scenario-specific regional-convergence input bundle
#'
#' Uses the scenario's retained sites and species matrix, then attaches the same
#' unscaled landscape values and IBRA7 subregions used by the main analysis. This
#' permits daily and protocol alternatives to gain or lose sites transparently.
#' @param prepared Scenario-specific LCBD preparation result.
#' @param landscape Unscaled selected landscape metrics.
#' @param hex_path,ibra_path Analysis grid and IBRA7 subregion layer.
#' @param config Main configuration.
#' @return Regional site data and a species matrix in exactly the same row order.
prepare_sensitivity_regional_inputs <- function(prepared, landscape, hex_path,
                                                ibra_path, config) {
  required <- c("hex_id", "pland_agriculture_2500m", "area_mn_woodland_2500m",
                "cohesion_woodland_2500m")
  if (!all(required %in% names(landscape))) {
    stop("Landscape metrics cannot supply regional sensitivity predictors.", call. = FALSE)
  }
  ids <- as.character(prepared$site_data$hex_id)
  site_data <- data.frame(
    hex_id = ids,
    noisy_miner_detection = prepared$site_data$miner_detection_rate,
    stringsAsFactors = FALSE
  )
  landscape$hex_id <- as.character(landscape$hex_id)
  landscape <- landscape[match(ids, as.character(landscape$hex_id)), required, drop = FALSE]
  names(landscape) <- c("hex_id", "agricultural_cover", "mean_woodland_patch_area",
                        "woodland_cohesion")
  site_data <- dplyr::left_join(site_data, landscape, by = "hex_id")
  site_data <- site_data[stats::complete.cases(site_data), , drop = FALSE]
  y <- prepared$species_matrix[site_data$hex_id, , drop = FALSE]

  grid <- sf::st_read(hex_path, quiet = TRUE)
  grid <- grid[as.character(grid$hex_id) %in% site_data$hex_id, "hex_id"]
  grid <- sf::st_transform(grid, config$regional_crs)
  ibra <- sf::st_transform(sf::st_read(ibra_path, quiet = TRUE), config$regional_crs)
  if (!"SUB_NAME_7" %in% names(ibra)) stop("The IBRA layer lacks SUB_NAME_7.", call. = FALSE)
  points <- suppressWarnings(sf::st_point_on_surface(grid))
  lookup <- sf::st_drop_geometry(sf::st_join(
    sf::st_sf(hex_id = as.character(grid$hex_id), geometry = sf::st_geometry(points)),
    ibra["SUB_NAME_7"], left = TRUE
  ))
  lookup <- unique(lookup)
  assert_unique_ids(lookup, "hex_id", "sensitivity IBRA subregion lookup")
  site_data <- dplyr::left_join(site_data, lookup, by = "hex_id")
  names(site_data)[names(site_data) == "SUB_NAME_7"] <- "ibra_subregion"
  site_data <- site_data[!is.na(site_data$ibra_subregion), , drop = FALSE]
  y <- y[site_data$hex_id, , drop = FALSE]
  y <- y[, colSums(y) > 0, drop = FALSE]
  if (nrow(y) < 2L || !ncol(y) || any(rowSums(y) == 0L)) {
    stop("Insufficient regional sensitivity sites or species after IBRA alignment.", call. = FALSE)
  }
  site_data$site_row <- seq_len(nrow(site_data))
  if (!identical(rownames(y), site_data$hex_id)) {
    stop("Regional sensitivity sites and species rows are misaligned.", call. = FALSE)
  }
  list(site_data = site_data, species_matrix = y)
}

#' Propagate a sampling scenario to regional convergence
#' @param prepared Alternative model-site data.
#' @param landscape Unscaled landscape data from the main analysis.
#' @param config Main settings, including 99 bootstrap draws.
#' @return Per-gradient status and regional endpoint contrasts.
fit_completeness_regional <- function(prepared, landscape, config) {
  regional <- prepare_sensitivity_regional_inputs(
    prepared, landscape, config$hex_2_5_path, config$ibra_path, config
  )
  dat <- regional$site_data
  y <- regional$species_matrix
  definitions <- define_regional_gradients(
    dat, config$regional_agriculture_low_max, config$regional_agriculture_high_min,
    config$regional_lower_tail_probability, config$regional_upper_tail_probability
  )
  components <- calculate_regional_jaccard_components(y)
  results <- statuses <- vector("list", nrow(definitions))
  for (i in seq_len(nrow(definitions))) {
    # Biological endpoint support is an informative sensitivity outcome. Other
    # errors signal malformed data and continue to stop the pipeline.
    result <- tryCatch(analyse_regional_gradient(
      definitions[i, ], dat, y, components$component_matrices, i,
      config$regional_min_sites_per_endpoint, config$regional_min_retained_subregions,
      config$regional_bootstrap_draws, config$regional_bootstrap_seed
    ), error = function(e) {
      if (!grepl("^Fewer than .* IBRA subregions contain at least", conditionMessage(e))) stop(e)
      list(unsupported = conditionMessage(e))
    })
    statuses[[i]] <- data.frame(
      gradient = definitions$gradient[i],
      status = if (is.null(result$unsupported)) "complete" else "unsupported",
      reason = if (is.null(result$unsupported)) "" else result$unsupported
    )
    results[[i]] <- result
  }
  supported <- Filter(function(x) is.null(x$unsupported), results)
  list(status = dplyr::bind_rows(statuses), definitions = definitions,
       contrasts = dplyr::bind_rows(lapply(supported, `[[`, "endpoint_contrasts")),
       dimensions = dplyr::bind_rows(lapply(supported, `[[`, "dimensions")))
}

#' Label regional-convergence sensitivity outputs
#' @param prepared Scenario-specific matrix and model sites.
#' @param baseline Main sensitivity input bundle.
#' @param config Main configuration.
#' @return Status and regional results with a stable scenario identifier.
fit_sensitivity_regional_scenario <- function(prepared, baseline, config) {
  # The completed main regional result is inserted directly into the output
  # tables. The extra `main_70` branch exists only to refit the LCBD ensemble
  # with the matched sensitivity allocation count.
  if (identical(prepared$scenario$type, "main")) {
    return(list(status = data.frame(), definitions = data.frame(),
                contrasts = data.frame(), dimensions = data.frame()))
  }
  if (identical(prepared$scenario$type, "pools")) {
    return(list(
      status = data.frame(gradient = NA_character_, status = "not_applicable",
        reason = "Reference-pool number affects LCBD only.", scenario_id = prepared$scenario$id),
      definitions = data.frame(), contrasts = data.frame(), dimensions = data.frame()
    ))
  }
  out <- fit_completeness_regional(prepared, baseline$landscape, config)
  for (name in c("status", "definitions", "contrasts", "dimensions")) {
    if (nrow(out[[name]])) out[[name]]$scenario_id <- prepared$scenario$id
  }
  out
}
