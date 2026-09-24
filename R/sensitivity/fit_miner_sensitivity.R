#' Fit the original residualised Noisy Miner LCBD comparison
#' @param allocations Exact main raw-detection pools with added Miner residuals.
#' @param main Main LCBD checkpoint, supplying the matched raw-detection fits.
#' @param config,settings Main and sensitivity settings.
#' @return Raw and residualised allocation curves, summaries and paired support.
fit_miner_sensitivity <- function(allocations, main, config, settings) {
  # Only the Noisy Miner smooth changes. Landscape, effort, GP, pool random
  # intercept, Gaussian family and all fitting options match the main models.
  s <- mgcv::s
  total <- relative_lcbd_total ~
    s(pland_agriculture_2500m, k = 5) + s(area_mn_woodland_2500m, k = 5) +
    s(cohesion_woodland_2500m, k = 5) + s(miner_residuals, k = 5) +
    log_effort + s(x_km, y_km, bs = "gp", k = 50) + s(reference_pool, bs = "re")
  formulas <- list(total = total,
    replacement = stats::update(total, relative_lcbd_replacement ~ .),
    richness_difference = stats::update(total, relative_lcbd_richness_difference ~ .))
  fits <- fit_lcbd_allocation_models(
    allocation_data = allocations$allocation_data,
    allocation_seeds = settings$allocation_seeds,
    focal_predictors = c(noisy_miner = "miner_residuals"),
    focal_labels = c(noisy_miner = "Landscape-residualised Noisy Miner detection"),
    formulas = formulas, prediction_grids = list(noisy_miner = allocations$grid),
    model_threads = config$model_threads
  )
  # Use the same first 30 canonical allocations for both raw and residualised
  # detection. The saved main checkpoint contains 100 allocations, but mixing
  # that full ensemble with a reduced sensitivity ensemble would be unmatched.
  allocation_ids <- seq_along(settings$allocation_seeds)
  partitions <- unique(as.data.frame(
    main$model_results_by_allocation[
      main$model_results_by_allocation$allocation %in% allocation_ids,
      c("allocation", "partition_id")
    ]
  ))
  if (nrow(partitions) != length(allocation_ids) ||
      !identical(as.integer(sort(partitions$allocation)), allocation_ids)) {
    stop("The main checkpoint does not contain every requested Miner allocation.",
         call. = FALSE)
  }
  ensemble <- summarise_lcbd_ensemble(
    fits$per_allocation_results, fits$per_allocation_curves, partitions,
    length(settings$allocation_seeds) * 3L
  )

  raw_allocation_results <- main$model_results_by_allocation[
    main$model_results_by_allocation$allocation %in% allocation_ids &
      main$model_results_by_allocation$focal_gradient == "noisy_miner",
    setdiff(names(main$model_results_by_allocation), "partition_id"),
    drop = FALSE
  ]
  raw_allocation_curves <- main$curves_by_allocation[
    main$curves_by_allocation$allocation %in% allocation_ids &
      main$curves_by_allocation$focal_gradient == "noisy_miner",
    setdiff(names(main$curves_by_allocation), "partition_id"),
    drop = FALSE
  ]
  raw_ensemble <- summarise_lcbd_ensemble(
    raw_allocation_results, raw_allocation_curves, partitions,
    length(settings$allocation_seeds) * 3L
  )

  raw_results <- raw_ensemble$main_results
  raw_results$scenario <- "Raw detection rate"
  residual_results <- ensemble$main_results
  residual_results$scenario <- "Landscape-residualised detection"
  raw_curves <- raw_ensemble$per_allocation_curves
  raw_curves$scenario <- "Raw detection rate"
  residual_curves <- ensemble$per_allocation_curves
  residual_curves$scenario <- "Landscape-residualised detection"
  curves <- dplyr::bind_rows(raw_curves, residual_curves)
  # Rates and deviance residuals have different units. Percent along each
  # central predictor range permits a SHAPE comparison, not a unit-effect ratio.
  curves <- dplyr::mutate(dplyr::group_by(curves, scenario, allocation, lcbd_component),
    gradient_position_percent = 100 * (gradient_value - min(gradient_value)) /
      (max(gradient_value) - min(gradient_value)))
  list(summary = dplyr::bind_rows(raw_results, residual_results),
       curves = dplyr::ungroup(curves), support = allocations$support)
}
