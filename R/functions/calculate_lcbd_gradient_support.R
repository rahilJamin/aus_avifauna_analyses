import::from(dplyr, filter, first, group_by, left_join, mutate, n, summarise, select)
import::from(magrittr, `%>%`)
import::from(tidyr, pivot_longer)
import::from(tidyselect, all_of)

#' Identify informative LCBD reference pools by gradient
#'
#' For each focal habitat/disturbance gradient, marks a reference pool as
#' informative when it has enough total sites and enough sites in both global
#' gradient tails. The returned gradient-specific datasets retain all sites from
#' informative pools, including intermediate-gradient sites.
#'
#' @param allocation_lcbd_data Site-level LCBD data for one allocation.
#' @param focal_predictors Named character vector of focal predictor columns.
#' @param global_gradient_thresholds Global low/high thresholds by predictor.
#' @param allocation Allocation index.
#' @param allocation_seed Seed used for the allocation.
#' @param minimum_pool_sites Minimum sites required in a reference pool.
#' @param minimum_tail_sites Minimum low-tail and high-tail sites required for a
#'   pool to support a gradient.
#'
#' @return A list with `gradient_support`, `gradient_pools`, and `gradient_data`.

calculate_lcbd_gradient_support <- function(allocation_lcbd_data,
                                            focal_predictors,
                                            global_gradient_thresholds,
                                            allocation,
                                            allocation_seed,
                                            minimum_pool_sites,
                                            minimum_tail_sites) {
  if (!is.character(focal_predictors) || length(focal_predictors) < 1L ||
      is.null(names(focal_predictors)) || anyDuplicated(focal_predictors) ||
      anyDuplicated(names(focal_predictors)) || anyNA(focal_predictors) ||
      any(!nzchar(focal_predictors)) || any(!nzchar(names(focal_predictors)))) {
    stop("Focal predictors must be uniquely named columns.", call. = FALSE)
  }
  required_columns <- c("hex_id", "reference_pool", unname(focal_predictors))
  missing_columns <- setdiff(required_columns, names(allocation_lcbd_data))
  if (length(missing_columns) > 0L) {
    stop(
      "The allocation LCBD table lacks: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  assert_unique_ids(allocation_lcbd_data, "hex_id", "allocation LCBD table")
  predictor_data <- allocation_lcbd_data[unname(focal_predictors)]
  if (!all(vapply(predictor_data, is.numeric, logical(1))) ||
      anyNA(predictor_data) || any(!is.finite(as.matrix(predictor_data))) ||
      anyNA(allocation_lcbd_data$reference_pool)) {
    stop("Allocation gradients and reference pools must be complete and finite.", call. = FALSE)
  }

  valid_count <- function(value) {
    is.numeric(value) && length(value) == 1L && !is.na(value) &&
      is.finite(value) && value >= 1 && value == floor(value)
  }
  if (!valid_count(minimum_pool_sites) || !valid_count(minimum_tail_sites)) {
    stop("LCBD pool and tail minimums must be positive whole numbers.", call. = FALSE)
  }

  required_threshold_columns <- c(
    "focal_predictor", "global_low_threshold", "global_high_threshold"
  )
  threshold_value_columns <- c("global_low_threshold", "global_high_threshold")
  if (!all(required_threshold_columns %in% names(global_gradient_thresholds)) ||
      anyDuplicated(global_gradient_thresholds$focal_predictor) ||
      !setequal(
        as.character(global_gradient_thresholds$focal_predictor),
        unname(focal_predictors)
      ) ||
      anyNA(global_gradient_thresholds[required_threshold_columns]) ||
      !all(vapply(
        global_gradient_thresholds[threshold_value_columns],
        is.numeric,
        logical(1)
      )) ||
      any(!is.finite(as.matrix(
        global_gradient_thresholds[threshold_value_columns]
      ))) ||
      any(global_gradient_thresholds$global_low_threshold >=
          global_gradient_thresholds$global_high_threshold)) {
    stop("Global LCBD gradient thresholds are incomplete or invalid.", call. = FALSE)
  }

  # A pool is informative for a focal gradient only when it has enough sites in
  # both global tails of that gradient. This prevents a fitted curve from being
  # driven by pools that cover only a narrow slice of the gradient.
  gradient_support <- allocation_lcbd_data %>%
    select(hex_id, reference_pool, all_of(unname(focal_predictors))) %>%
    pivot_longer(
      cols = all_of(unname(focal_predictors)),
      names_to = "focal_predictor",
      values_to = "predictor_value"
    ) %>%
    left_join(global_gradient_thresholds, by = "focal_predictor") %>%
    group_by(focal_predictor, reference_pool) %>%
    summarise(
      n_sites = n(),
      global_low_threshold = first(global_low_threshold),
      global_high_threshold = first(global_high_threshold),
      within_pool_minimum = min(predictor_value),
      within_pool_maximum = max(predictor_value),
      n_low = sum(predictor_value <= first(global_low_threshold)),
      n_high = sum(predictor_value >= first(global_high_threshold)),
      .groups = "drop"
    ) %>%
    mutate(
      informative_pool =
        n_sites >= minimum_pool_sites &
        n_low >= minimum_tail_sites &
        n_high >= minimum_tail_sites,
      allocation,
      seed = allocation_seed,
      .before = 1
    )

  gradient_pools <- lapply(focal_predictors, function(focal_column) {
    gradient_support$reference_pool[
      gradient_support$focal_predictor == focal_column &
        gradient_support$informative_pool
    ] %>%
      as.character()
  })

  names(gradient_pools) <- names(focal_predictors)

  # All sites from informative pools are retained, including intermediate sites.
  # This preserves the original model-support definition.
  gradient_data <- lapply(
    gradient_pools,
    function(pools) {
      allocation_lcbd_data %>%
        filter(as.character(reference_pool) %in% pools) %>%
        mutate(reference_pool = droplevels(factor(reference_pool)))
    }
  )

  list(
    gradient_support = gradient_support,
    gradient_pools = gradient_pools,
    gradient_data = gradient_data
  )
}
