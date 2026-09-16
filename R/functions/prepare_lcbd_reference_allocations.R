import::from(dplyr, bind_rows, distinct, filter, group_by, left_join, mutate,
             n_distinct, relocate, row_number, select, summarise)
import::from(magrittr, `%>%`)

#' Prepare all LCBD reference-pool allocations
#'
#' Coordinates the LCBD allocation workflow across seeds: assign spatial
#' reference pools, calculate within-pool LCBD, identify gradient-specific
#' informative pools, summarise allocation support, validate minimum support, and
#' calculate shared prediction ranges. Repeating allocations checks whether LCBD
#' patterns are robust to the exact spatial pooling of nearby sites.
#'
#' @param site_data Aligned LCBD site table with focal predictors and projected
#'   kilometre coordinates.
#' @param distance_matrices List with total, replacement, and richness-difference
#'   dissimilarity matrices.
#' @param focal_predictors Named character vector of focal predictor columns.
#' @param global_gradient_thresholds Global low/high gradient thresholds.
#' @param allocation_seeds Integer seeds for repeated reference-pool allocations.
#' @param reference_pool_k Number of spatial reference pools per allocation.
#' @param reference_pool_nstart Number of random starts passed to `kmeans()`.
#' @param minimum_pool_sites Minimum sites required for a complete reference pool.
#' @param minimum_tail_sites Minimum sites required in both gradient tails.
#' @param minimum_informative_pools Minimum informative pools required before
#'   model fitting.
#' @param curve_lower_probability Lower quantile for central prediction ranges.
#' @param curve_upper_probability Upper quantile for central prediction ranges.
#' @param negative_tolerance Threshold for meaningful negative LCBD values.
#' @param lcbd_sum_tolerance Maximum allowed raw LCBD sum error.
#'
#' @return A list containing allocation-specific model data, support summaries,
#'   validation tables, membership tables, and shared curve limits.

prepare_lcbd_reference_allocations <- function(site_data,
                                               distance_matrices,
                                               focal_predictors,
                                               global_gradient_thresholds,
                                               allocation_seeds,
                                               reference_pool_k,
                                               reference_pool_nstart,
                                               minimum_pool_sites,
                                               minimum_tail_sites,
                                               minimum_informative_pools,
                                               curve_lower_probability,
                                               curve_upper_probability,
                                               negative_tolerance,
                                               lcbd_sum_tolerance) {
  n_reference_pool_allocations <- length(allocation_seeds)

  allocation_data <- vector("list", n_reference_pool_allocations)
  allocation_summary_rows <- vector("list", n_reference_pool_allocations)
  model_support_rows <- vector("list", n_reference_pool_allocations)
  pool_lcbd_check_rows <- vector("list", n_reference_pool_allocations)
  gradient_support_rows <- vector("list", n_reference_pool_allocations)
  allocation_membership_rows <- vector("list", n_reference_pool_allocations)

  for (allocation in seq_len(n_reference_pool_allocations)) {
    allocation_seed <- allocation_seeds[allocation]
    message(
      "Preparing reference-pool allocation ",
      allocation,
      " of ",
      n_reference_pool_allocations,
      " (seed ",
      allocation_seed,
      ")"
    )

    pool_assignment <- assign_lcbd_reference_pools(
      site_data = site_data,
      allocation = allocation,
      allocation_seed = allocation_seed,
      reference_pool_k = reference_pool_k,
      reference_pool_nstart = reference_pool_nstart,
      minimum_pool_sites = minimum_pool_sites
    )

    lcbd_by_pool <- calculate_reference_pool_lcbd(
      complete_pool_data = pool_assignment$complete_pool_data,
      distance_matrices = distance_matrices,
      allocation = allocation,
      allocation_seed = allocation_seed,
      negative_tolerance = negative_tolerance,
      lcbd_sum_tolerance = lcbd_sum_tolerance
    )

    gradient_support <- calculate_lcbd_gradient_support(
      allocation_lcbd_data = lcbd_by_pool$allocation_lcbd_data,
      focal_predictors = focal_predictors,
      global_gradient_thresholds = global_gradient_thresholds,
      allocation = allocation,
      allocation_seed = allocation_seed,
      minimum_pool_sites = minimum_pool_sites,
      minimum_tail_sites = minimum_tail_sites
    )

    allocation_data[[allocation]] <- gradient_support$gradient_data

    allocation_support <- summarise_lcbd_allocation_support(
      allocation = allocation,
      allocation_seed = allocation_seed,
      allocation_site_data = pool_assignment$allocation_site_data,
      geographical_kmeans = pool_assignment$geographical_kmeans,
      pool_sizes = pool_assignment$pool_sizes,
      pool_checks = lcbd_by_pool$pool_checks,
      gradient_data = gradient_support$gradient_data,
      gradient_pools = gradient_support$gradient_pools,
      focal_predictors = focal_predictors,
      curve_lower_probability = curve_lower_probability,
      curve_upper_probability = curve_upper_probability
    )

    allocation_summary_rows[[allocation]] <- allocation_support$allocation_summary
    model_support_rows[[allocation]] <- allocation_support$model_support
    pool_lcbd_check_rows[[allocation]] <- lcbd_by_pool$pool_checks
    gradient_support_rows[[allocation]] <- gradient_support$gradient_support
    allocation_membership_rows[[allocation]] <- pool_assignment$allocation_membership
  }

  allocation_summary <- bind_rows(allocation_summary_rows)
  model_support <- bind_rows(model_support_rows)
  pool_lcbd_checks <- bind_rows(pool_lcbd_check_rows)
  gradient_support_by_pool <- bind_rows(gradient_support_rows)
  allocation_membership <- bind_rows(allocation_membership_rows)

  partition_lookup <- allocation_summary %>%
    distinct(partition_signature) %>%
    mutate(partition_id = row_number())

  allocation_summary <- allocation_summary %>%
    left_join(partition_lookup, by = "partition_signature") %>%
    relocate(partition_id, .after = seed)

  insufficient_support <- model_support %>%
    filter(n_informative_pools < minimum_informative_pools)

  allocation_validation <- allocation_summary %>%
    select(-partition_signature) %>%
    mutate(
      n_unique_partitions_across_all_allocations = n_distinct(partition_id),
      all_lcbd_sums_valid =
        maximum_absolute_lcbd_sum_error <= lcbd_sum_tolerance,
      all_four_gradients_supported =
        agriculture_informative_pools >= minimum_informative_pools &
        patch_area_informative_pools >= minimum_informative_pools &
        cohesion_informative_pools >= minimum_informative_pools &
        noisy_miner_informative_pools >= minimum_informative_pools
    )

  if (nrow(insufficient_support) > 0L) {
    print(insufficient_support)
    stop(
      paste(
        "At least one allocation has insufficient focal-gradient support.",
        "The failing support rows are printed above; no GAMMs were fitted."
      ),
      call. = FALSE
    )
  }

  common_curve_limits <- model_support %>%
    group_by(focal_gradient) %>%
    summarise(
      curve_lower = max(central_lower),
      curve_upper = min(central_upper),
      .groups = "drop"
    )

  if (any(common_curve_limits$curve_lower >= common_curve_limits$curve_upper)) {
    print(common_curve_limits)
    stop("A focal gradient has no shared central prediction range.", call. = FALSE)
  }

  list(
    allocation_data = allocation_data,
    allocation_summary = allocation_summary,
    model_support = model_support,
    pool_lcbd_checks = pool_lcbd_checks,
    gradient_support_by_pool = gradient_support_by_pool,
    allocation_membership = allocation_membership,
    allocation_validation = allocation_validation,
    common_curve_limits = common_curve_limits
  )
}
