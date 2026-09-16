import::from(dplyr, bind_rows, filter, n_distinct)
import::from(magrittr, `%>%`)
import::from(stats, median, quantile)
import::from(tibble, tibble)

#' Summarise LCBD allocation support
#'
#' Builds allocation-level diagnostics and gradient-specific model-support rows,
#' including complete-pool sizes, LCBD validation counts, informative-pool counts,
#' and central curve ranges used to define shared prediction grids.
#'
#' @param allocation Allocation index.
#' @param allocation_seed Seed used for the allocation.
#' @param allocation_site_data Site data with reference-pool labels for all sites.
#' @param geographical_kmeans K-means object returned by
#'   `assign_lcbd_reference_pools()`.
#' @param pool_sizes Pool-size summary from `assign_lcbd_reference_pools()`.
#' @param pool_checks LCBD sum/negative-value checks from
#'   `calculate_reference_pool_lcbd()`.
#' @param gradient_data Named list of gradient-specific modelling datasets.
#' @param gradient_pools Named list of informative reference pools by gradient.
#' @param focal_predictors Named character vector of focal predictor columns.
#' @param curve_lower_probability Lower quantile for the central prediction
#'   range.
#' @param curve_upper_probability Upper quantile for the central prediction
#'   range.
#'
#' @return A list with `model_support` and `allocation_summary`.

summarise_lcbd_allocation_support <- function(allocation,
                                              allocation_seed,
                                              allocation_site_data,
                                              geographical_kmeans,
                                              pool_sizes,
                                              pool_checks,
                                              gradient_data,
                                              gradient_pools,
                                              focal_predictors,
                                              curve_lower_probability,
                                              curve_upper_probability) {
  model_support <- bind_rows(lapply(
    names(focal_predictors),
    function(focal_gradient) {
      focal_data <- gradient_data[[focal_gradient]]
      focal_variable <- focal_predictors[[focal_gradient]]

      tibble(
        allocation,
        seed = allocation_seed,
        focal_gradient = focal_gradient,
        n_sites = nrow(focal_data),
        n_informative_pools = n_distinct(focal_data$reference_pool),
        central_lower = if (nrow(focal_data) > 0L) {
          quantile(
            focal_data[[focal_variable]],
            curve_lower_probability,
            na.rm = TRUE,
            names = FALSE
          )
        } else {
          NA_real_
        },
        central_upper = if (nrow(focal_data) > 0L) {
          quantile(
            focal_data[[focal_variable]],
            curve_upper_probability,
            na.rm = TRUE,
            names = FALSE
          )
        } else {
          NA_real_
        }
      )
    }
  ))

  ordered_pool_labels <- allocation_site_data$reference_pool[
    order(allocation_site_data$hex_id)
  ]
  partition_signature <- paste(ordered_pool_labels, collapse = ",")

  complete_pool_sizes <- pool_sizes$reference_pool_size[pool_sizes$complete_pool]

  allocation_summary <- tibble(
    allocation,
    seed = allocation_seed,
    partition_signature,
    total_within_cluster_sum_of_squares = geographical_kmeans$tot.withinss,
    n_complete_pools = sum(pool_sizes$complete_pool),
    minimum_complete_pool_size = min(complete_pool_sizes),
    median_complete_pool_size = median(complete_pool_sizes),
    maximum_complete_pool_size = max(complete_pool_sizes),
    maximum_absolute_lcbd_sum_error = max(abs(pool_checks$raw_lcbd_sum - 1)),
    n_meaningful_negative_total = sum(
      pool_checks$n_values_below_negative_tolerance[
        pool_checks$lcbd_component == "total"
      ]
    ),
    n_meaningful_negative_replacement = sum(
      pool_checks$n_values_below_negative_tolerance[
        pool_checks$lcbd_component == "replacement"
      ]
    ),
    n_meaningful_negative_richness_difference = sum(
      pool_checks$n_values_below_negative_tolerance[
        pool_checks$lcbd_component == "richness_difference"
      ]
    ),
    agriculture_informative_pools = length(gradient_pools$agriculture),
    agriculture_sites = nrow(gradient_data$agriculture),
    patch_area_informative_pools = length(gradient_pools$patch_area),
    patch_area_sites = nrow(gradient_data$patch_area),
    cohesion_informative_pools = length(gradient_pools$cohesion),
    cohesion_sites = nrow(gradient_data$cohesion),
    noisy_miner_informative_pools = length(gradient_pools$noisy_miner),
    noisy_miner_sites = nrow(gradient_data$noisy_miner)
  )

  list(
    model_support = model_support,
    allocation_summary = allocation_summary
  )
}
