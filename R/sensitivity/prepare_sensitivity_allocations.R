#' Prepare matched-size LCBD ensembles and retain unsupported scenarios
#'
#' Uses the same pool assignment, LCBD calculation and gradient-support functions
#' as the main pipeline. A scenario that cannot support ALL planned allocations
#' is reported as unsupported, rather than being rescued by dropping allocations
#' or relaxing the minimum-pool rule.
#' @param prepared Output of prepare_completeness_data().
#' @param config Main settings.
#' @param settings Sensitivity settings with the matched reduced allocation seeds.
#' @return Allocation data, support tables, prediction limits and a status.
prepare_sensitivity_allocations <- function(prepared, config, settings) {
  thresholds <- calculate_lcbd_gradient_thresholds(
    prepared$site_data, config$lcbd_focal_predictors,
    config$lcbd_lower_tail_probability, config$lcbd_upper_tail_probability
  )
  allocation_data <- support_rows <- summary_rows <- vector("list", length(settings$allocation_seeds))
  pool_k <- prepared$scenario$pool_k
  for (i in seq_along(settings$allocation_seeds)) {
    seed <- settings$allocation_seeds[i]
    message(prepared$scenario$label, ": allocation ", i,
            " of ", length(settings$allocation_seeds))
    pools <- tryCatch(
      assign_lcbd_reference_pools(
        prepared$site_data, i, seed, pool_k,
        config$lcbd_reference_pool_nstart, config$lcbd_minimum_pool_sites
      ),
      error = function(e) {
        # These two errors describe inadequate scenario support. Other errors
        # remain hard failures because they indicate a data or code problem.
        if (grepl("^(No reference pool meets|The number of reference pools exceeds)",
                  conditionMessage(e))) {
          return(structure(list(message = conditionMessage(e)), class = "unsupported_pools"))
        }
        stop(e)
      }
    )
    if (inherits(pools, "unsupported_pools") || !nrow(pools$complete_pool_data)) {
      reason <- if (inherits(pools, "unsupported_pools")) pools$message else
        paste("No complete pools in allocation", i)
      return(list(status = "unsupported", reason = reason,
                  support = dplyr::bind_rows(support_rows)))
    }
    lcbd <- calculate_reference_pool_lcbd(
      pools$complete_pool_data, prepared$distance_matrices, i, seed,
      config$lcbd_negative_tolerance, config$lcbd_sum_tolerance
    )
    support <- calculate_lcbd_gradient_support(
      lcbd$allocation_lcbd_data, config$lcbd_focal_predictors, thresholds,
      i, seed, config$lcbd_minimum_pool_sites, config$lcbd_minimum_tail_sites
    )
    allocation_data[[i]] <- support$gradient_data
    summary <- summarise_lcbd_allocation_support(
      i, seed, pools$allocation_site_data, pools$geographical_kmeans,
      pools$pool_sizes, lcbd$pool_checks, support$gradient_data,
      support$gradient_pools, config$lcbd_focal_predictors,
      config$lcbd_curve_lower_probability, config$lcbd_curve_upper_probability
    )
    support_rows[[i]] <- summary$model_support
    summary_rows[[i]] <- summary$allocation_summary
  }
  model_support <- dplyr::bind_rows(support_rows)
  allocation_summary <- dplyr::bind_rows(summary_rows)
  allocation_summary$partition_id <- match(
    allocation_summary$partition_signature, unique(allocation_summary$partition_signature)
  )
  if (any(model_support$n_informative_pools < config$lcbd_minimum_informative_pools)) {
    return(list(status = "unsupported", reason = "At least one gradient/allocation has too few informative pools.",
                support = model_support))
  }
  limits <- dplyr::summarise(
    dplyr::group_by(model_support, focal_gradient),
    curve_lower = max(central_lower), curve_upper = min(central_upper), .groups = "drop"
  )
  if (any(!is.finite(as.matrix(limits[c("curve_lower", "curve_upper")]))) ||
      any(limits$curve_lower >= limits$curve_upper)) {
    return(list(status = "unsupported", reason = "No shared central prediction range across allocations.",
                support = model_support))
  }
  list(status = "ready", reason = "", allocation_data = allocation_data,
       allocation_summary = allocation_summary, support = model_support,
       common_curve_limits = limits)
}

#' Fit one supported completeness ensemble with unchanged main GAMM formulae
#' @param allocations Prepared allocations, including unsupported status.
#' @param config,settings Main and sensitivity settings.
#' @return Completed ensemble summaries or an explicit unsupported result.
fit_completeness_lcbd <- function(allocations, config, settings) {
  if (allocations$status != "ready") return(allocations)
  fits <- fit_lcbd_allocation_models(
    allocation_data = allocations$allocation_data,
    allocation_seeds = settings$allocation_seeds,
    focal_predictors = config$lcbd_focal_predictors,
    focal_labels = config$lcbd_focal_labels,
    formulas = define_lcbd_formulas(),
    prediction_grids = make_lcbd_prediction_grids(
      allocations$common_curve_limits, config$lcbd_curve_grid_points
    ), model_threads = config$model_threads
  )
  ensemble <- summarise_lcbd_ensemble(
    fits$per_allocation_results, fits$per_allocation_curves,
    allocations$allocation_summary, length(settings$allocation_seeds) * 12L
  )
  c(list(status = "complete", reason = "", support = allocations$support), ensemble)
}

#' Run and label one LCBD sensitivity scenario
#' @param prepared Scenario-specific LCBD inputs.
#' @param config,settings Main and sensitivity settings.
#' @return A compact result that always contains a scenario status row.
fit_sensitivity_lcbd_scenario <- function(prepared, config, settings) {
  fitted <- fit_completeness_lcbd(
    prepare_sensitivity_allocations(prepared, config, settings), config, settings
  )
  status <- data.frame(
    scenario_id = prepared$scenario$id,
    scenario_label = prepared$scenario$label,
    scenario_type = prepared$scenario$type,
    completeness = prepared$scenario$completeness,
    minimum_sampling_units = prepared$scenario$minimum_units,
    reference_pool_k = prepared$scenario$pool_k,
    n_chao_sites = prepared$n_chao_sites,
    n_model_sites = nrow(prepared$site_data),
    n_response_species = ncol(prepared$species_matrix),
    status = fitted$status,
    reason = fitted$reason,
    n_allocations_requested = length(settings$allocation_seeds),
    stringsAsFactors = FALSE
  )
  if (!identical(fitted$status, "complete")) {
    support <- fitted$support
    if (nrow(support)) support$scenario_id <- prepared$scenario$id
    return(list(status = status, support = support,
                results = data.frame(), curves = data.frame()))
  }
  results <- fitted$main_results
  curves <- fitted$main_curves
  support <- fitted$support
  results$scenario_id <- prepared$scenario$id
  curves$scenario_id <- prepared$scenario$id
  support$scenario_id <- prepared$scenario$id
  list(status = status, support = support, results = results, curves = curves)
}
