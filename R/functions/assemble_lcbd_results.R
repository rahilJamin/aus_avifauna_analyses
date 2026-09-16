#' Assemble the main LCBD checkpoint
#'
#' Collects the validated inputs, allocation support, fitted-model summaries,
#' prediction curves, and representative models into the single documented
#' checkpoint consumed by the LCBD output stage.
#'
#' @param config Central workflow configuration list.
#' @param lcbd_inputs Aligned LCBD inputs returned by `load_lcbd_inputs()`.
#' @param distance_validation Comparison of saved and rebuilt beta matrices.
#' @param gradient_thresholds Fixed global focal-gradient thresholds.
#' @param allocation_bundle Reference-pool allocation results.
#' @param model_bundle Allocation-specific GAMM results and curves.
#' @param ensemble_bundle Summaries across reference-pool allocations.
#'
#' @return The organised main LCBD result list.
assemble_lcbd_results <- function(config,
                                  lcbd_inputs,
                                  distance_validation,
                                  gradient_thresholds,
                                  allocation_bundle,
                                  model_bundle,
                                  ensemble_bundle) {
  allocation_seeds <- config$lcbd_allocation_seed_start +
    seq_len(config$lcbd_n_allocations) - 1L

  list(
    settings = list(
      reference_pool_k = config$lcbd_reference_pool_k,
      n_reference_pool_allocations = config$lcbd_n_allocations,
      allocation_seeds = allocation_seeds,
      reference_pool_nstart = config$lcbd_reference_pool_nstart,
      minimum_pool_sites = config$lcbd_minimum_pool_sites,
      minimum_tail_sites = config$lcbd_minimum_tail_sites,
      minimum_informative_pools = config$lcbd_minimum_informative_pools,
      lower_tail_probability = config$lcbd_lower_tail_probability,
      upper_tail_probability = config$lcbd_upper_tail_probability,
      negative_tolerance = config$lcbd_negative_tolerance,
      lcbd_sum_tolerance = config$lcbd_sum_tolerance,
      curve_lower_probability = config$lcbd_curve_lower_probability,
      curve_upper_probability = config$lcbd_curve_upper_probability,
      model_family = "Gaussian identity for all three LCBD components",
      response_standardisation =
        "relative_lcbd = raw_lcbd * reference_pool_size",
      crs = "GDA94 / Australian Albers (EPSG:3577)"
    ),
    dimensions = tibble::tibble(
      n_lcbd_sites = lcbd_inputs$lcbd_site_count,
      n_lcbd_species = lcbd_inputs$lcbd_species_count,
      n_allocations = config$lcbd_n_allocations,
      n_unique_partitions = dplyr::n_distinct(
        allocation_bundle$allocation_summary$partition_id
      ),
      n_fitted_gamms = nrow(ensemble_bundle$per_allocation_results)
    ),
    distance_validation = distance_validation,
    analysis_site_data = lcbd_inputs$site_data,
    global_gradient_thresholds = gradient_thresholds,
    common_curve_limits = allocation_bundle$common_curve_limits,
    allocation_validation = allocation_bundle$allocation_validation,
    allocation_membership = allocation_bundle$allocation_membership,
    pool_lcbd_checks = allocation_bundle$pool_lcbd_checks,
    gradient_support_by_pool = allocation_bundle$gradient_support_by_pool,
    model_support_by_allocation = allocation_bundle$model_support,
    model_results_by_allocation = ensemble_bundle$per_allocation_results,
    curves_by_allocation = ensemble_bundle$per_allocation_curves,
    main_results = ensemble_bundle$main_results,
    main_curves = ensemble_bundle$main_curves,
    representative_models = model_bundle$representative_models
  )
}
