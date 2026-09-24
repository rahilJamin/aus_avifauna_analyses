#' Reuse the main Noisy Miner pools and calculate their unchanged LCBD responses
#' @param residualised Count-model result with one residual per main candidate site.
#' @param main Main LCBD ensemble checkpoint.
#' @param species_matrix Main retained presence-absence matrix.
#' @param config,settings Main and sensitivity settings.
#' @return Matched allocation data and a common residual prediction grid.
prepare_miner_allocations <- function(residualised, main, species_matrix, config, settings) {
  dat <- residualised$analysis_site_data
  if (!identical(as.character(dat$hex_id), as.character(main$analysis_site_data$hex_id))) {
    stop("Residualisation changed the main LCBD site order.", call. = FALSE)
  }
  if (!all(dat$hex_id %in% rownames(species_matrix))) stop("LCBD sites are missing from the response matrix.", call. = FALSE)
  y <- species_matrix[dat$hex_id, colnames(species_matrix) != config$noisy_miner_species, drop = FALSE]
  y <- y[, colSums(y) > 0, drop = FALSE]
  distances <- calculate_lcbd_beta_matrices(y)$dist_mats
  allocations <- support <- vector("list", length(settings$allocation_seeds))
  for (i in seq_along(allocations)) {
    membership <- main$allocation_membership
    membership <- membership[membership$allocation == i & membership$complete_pool, ]
    assert_unique_ids(membership, "hex_id", "saved pool membership")
    current <- dat[match(membership$hex_id, dat$hex_id), , drop = FALSE]
    current$reference_pool <- factor(membership$reference_pool)
    current$reference_pool_size <- as.integer(table(current$reference_pool)[current$reference_pool])
    current <- current[order(current$reference_pool, current$hex_id), , drop = FALSE]

    # LCBD must be calculated within ALL sites in each complete pool first.
    # Selecting only a residual tail before this calculation would change Y.
    lcbd <- calculate_reference_pool_lcbd(
      current, distances, i, settings$allocation_seeds[i],
      config$lcbd_negative_tolerance, config$lcbd_sum_tolerance
    )$allocation_lcbd_data
    pools <- main$gradient_support_by_pool
    pools <- pools[pools$allocation == i & pools$focal_predictor == "miner_detection_rate" &
                     pools$informative_pool, ]
    lcbd <- lcbd[as.character(lcbd$reference_pool) %in% as.character(pools$reference_pool), ]
    lcbd$reference_pool <- droplevels(factor(lcbd$reference_pool))
    expected <- main$model_results_by_allocation
    expected <- expected[expected$allocation == i & expected$focal_gradient == "noisy_miner", ]
    if (nrow(expected) != 3L || any(expected$n_sites != nrow(lcbd)) ||
        any(expected$n_informative_pools != nlevels(lcbd$reference_pool))) {
      stop("Reconstructed Noisy Miner support differs from the saved main fits.", call. = FALSE)
    }
    allocations[[i]] <- list(noisy_miner = lcbd)
    limits <- stats::quantile(lcbd$miner_residuals,
      c(config$lcbd_curve_lower_probability, config$lcbd_curve_upper_probability), names = FALSE)
    support[[i]] <- data.frame(allocation = i, seed = settings$allocation_seeds[i],
      n_sites = nrow(lcbd), n_informative_pools = nlevels(lcbd$reference_pool),
      central_lower = limits[1], central_upper = limits[2])
  }
  support <- dplyr::bind_rows(support)
  lower <- max(support$central_lower)
  upper <- min(support$central_upper)
  if (!is.finite(lower) || !is.finite(upper) || lower >= upper) {
    stop("No shared central Noisy Miner residual range across the matched allocations.", call. = FALSE)
  }
  list(allocation_data = allocations, support = support,
       grid = seq(lower, upper, length.out = config$lcbd_curve_grid_points))
}
