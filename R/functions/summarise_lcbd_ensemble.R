import::from(dplyr, case_when, first, group_by, left_join, mutate, n,
             n_distinct, select, summarise)
import::from(magrittr, `%>%`)
import::from(stats, median, sd)

#' Summarise the LCBD allocation ensemble
#'
#' Combines model summaries and prediction curves across reference-pool
#' allocations. Reported uncertainty combines mean within-allocation uncertainty
#' with between-allocation variation.
#'
#' @param per_allocation_results Extracted model summaries for every allocation,
#'   gradient, and LCBD component.
#' @param per_allocation_curves Extracted partial-effect curves for every
#'   allocation, gradient, and LCBD component.
#' @param allocation_summary Allocation summary table containing `partition_id`.
#' @param expected_model_count Expected number of fitted model summaries.
#'
#' @return A list with allocation-level results/curves and ensemble-level
#'   `main_results` and `main_curves`.

summarise_lcbd_ensemble <- function(per_allocation_results,
                                    per_allocation_curves,
                                    allocation_summary,
                                    expected_model_count) {
  required_result_columns <- c(
    "allocation", "focal_gradient", "focal_label", "focal_variable",
    "lcbd_component", "n_sites", "n_informative_pools", "endpoint_low",
    "endpoint_high", "endpoint_difference", "endpoint_se", "edf",
    "deviance_explained", "model_converged"
  )
  required_curve_columns <- c(
    "allocation", "focal_gradient", "focal_label", "focal_variable",
    "lcbd_component", "gradient_value", "partial_effect", "partial_effect_se"
  )
  required_allocation_columns <- c("allocation", "partition_id")

  missing_results <- setdiff(required_result_columns, names(per_allocation_results))
  missing_curves <- setdiff(required_curve_columns, names(per_allocation_curves))
  missing_allocations <- setdiff(required_allocation_columns, names(allocation_summary))
  if (length(c(missing_results, missing_curves, missing_allocations)) > 0L) {
    stop("The LCBD ensemble inputs are missing required columns.", call. = FALSE)
  }
  if (!is.numeric(expected_model_count) || length(expected_model_count) != 1L ||
      is.na(expected_model_count) || !is.finite(expected_model_count) ||
      expected_model_count < 1L || expected_model_count != floor(expected_model_count)) {
    stop("`expected_model_count` must be a positive whole number.", call. = FALSE)
  }
  if (nrow(allocation_summary) < 2L || anyNA(allocation_summary[required_allocation_columns]) ||
      anyDuplicated(allocation_summary$allocation)) {
    stop("The LCBD allocation summary must contain unique repeated allocations.", call. = FALSE)
  }

  result_numeric_columns <- c(
    "n_sites", "n_informative_pools", "endpoint_low", "endpoint_high",
    "endpoint_difference", "endpoint_se", "edf", "deviance_explained"
  )
  curve_numeric_columns <- c(
    "gradient_value", "partial_effect", "partial_effect_se"
  )
  if (anyNA(per_allocation_results[required_result_columns]) ||
      anyNA(per_allocation_curves[required_curve_columns]) ||
      !all(vapply(
        per_allocation_results[result_numeric_columns], is.numeric, logical(1)
      )) ||
      !all(vapply(
        per_allocation_curves[curve_numeric_columns], is.numeric, logical(1)
      )) ||
      any(!is.finite(as.matrix(per_allocation_results[result_numeric_columns]))) ||
      any(!is.finite(as.matrix(per_allocation_curves[curve_numeric_columns]))) ||
      any(per_allocation_results$endpoint_se < 0) ||
      any(per_allocation_curves$partial_effect_se < 0)) {
    stop("LCBD ensemble results or curves contain invalid values.", call. = FALSE)
  }

  allocation_ids <- allocation_summary$allocation
  if (!setequal(unique(per_allocation_results$allocation), allocation_ids) ||
      !setequal(unique(per_allocation_curves$allocation), allocation_ids)) {
    stop("LCBD results and curves do not cover the same allocations.", call. = FALSE)
  }

  result_key <- with(
    per_allocation_results,
    paste(allocation, focal_gradient, focal_variable, lcbd_component, sep = "\r")
  )
  curve_key <- with(
    per_allocation_curves,
    paste(
      allocation, focal_gradient, focal_variable, lcbd_component,
      format(gradient_value, digits = 17),
      sep = "\r"
    )
  )
  if (anyDuplicated(result_key) || anyDuplicated(curve_key)) {
    stop("LCBD ensemble inputs contain duplicated model or curve rows.", call. = FALSE)
  }

  result_group <- with(
    per_allocation_results,
    paste(focal_gradient, focal_variable, lcbd_component, sep = "\r")
  )
  curve_group <- with(
    per_allocation_curves,
    paste(focal_gradient, focal_variable, lcbd_component, sep = "\r")
  )
  n_allocations <- length(allocation_ids)
  if (any(table(result_group) != n_allocations) ||
      !setequal(unique(result_group), unique(curve_group))) {
    stop("Each LCBD model group must contain one result per allocation.", call. = FALSE)
  }
  curve_point_group <- paste(
    curve_group,
    format(per_allocation_curves$gradient_value, digits = 17),
    sep = "\r"
  )
  if (any(table(curve_point_group) != n_allocations)) {
    stop("LCBD prediction grids differ among allocations.", call. = FALSE)
  }

  endpoint_group <- split(per_allocation_results, result_group)
  endpoints_are_shared <- vapply(
    endpoint_group,
    function(rows) {
      length(unique(rows$endpoint_low)) == 1L &&
        length(unique(rows$endpoint_high)) == 1L &&
        rows$endpoint_low[1] < rows$endpoint_high[1]
    },
    logical(1)
  )
  if (!all(endpoints_are_shared)) {
    stop("LCBD model endpoints differ among allocations or are not ordered.", call. = FALSE)
  }

  # Attach canonical partition IDs before combining within-model uncertainty
  # with allocation-to-allocation variation.
  per_allocation_results <- per_allocation_results %>%
    left_join(
      select(allocation_summary, allocation, partition_id),
      by = "allocation"
    )

  per_allocation_curves <- per_allocation_curves %>%
    left_join(
      select(allocation_summary, allocation, partition_id),
      by = "allocation"
    )

  if (anyNA(per_allocation_results$partition_id) ||
      anyNA(per_allocation_curves$partition_id)) {
    stop("At least one LCBD result could not be linked to a partition.", call. = FALSE)
  }

  if (nrow(per_allocation_results) != expected_model_count) {
    stop(
      "Expected ",
      expected_model_count,
      " fitted model summaries but found ",
      nrow(per_allocation_results),
      ".",
      call. = FALSE
    )
  }

  ensemble_results <- per_allocation_results %>%
    group_by(
      focal_gradient,
      focal_label,
      focal_variable,
      lcbd_component
    ) %>%
    summarise(
      n_allocations = n(),
      n_unique_partitions = n_distinct(partition_id),
      minimum_n_sites = min(n_sites),
      mean_n_sites = mean(n_sites),
      maximum_n_sites = max(n_sites),
      minimum_n_informative_pools = min(n_informative_pools),
      mean_n_informative_pools = mean(n_informative_pools),
      maximum_n_informative_pools = max(n_informative_pools),
      endpoint_low = first(endpoint_low),
      endpoint_high = first(endpoint_high),
      mean_endpoint_difference = mean(endpoint_difference),
      mean_within_allocation_variance = mean(endpoint_se^2),
      between_allocation_sd = sd(endpoint_difference),
      total_se = sqrt(
        mean_within_allocation_variance +
          between_allocation_sd^2
      ),
      lower_95 = mean_endpoint_difference - 1.96 * total_se,
      upper_95 = mean_endpoint_difference + 1.96 * total_se,
      proportion_positive = mean(endpoint_difference > 0),
      proportion_negative = mean(endpoint_difference < 0),
      median_edf = median(edf),
      median_deviance_explained = median(deviance_explained),
      all_models_converged = all(model_converged),
      .groups = "drop"
    ) %>%
    mutate(
      direction_agreement = pmax(proportion_positive, proportion_negative),
      interpretation = case_when(
        lower_95 > 0 ~ "higher relative LCBD at the high endpoint",
        upper_95 < 0 ~ "lower relative LCBD at the high endpoint",
        TRUE ~ "endpoint difference uncertain"
      )
    )

  ensemble_curves <- per_allocation_curves %>%
    group_by(
      focal_gradient,
      focal_label,
      focal_variable,
      lcbd_component,
      gradient_value
    ) %>%
    summarise(
      n_allocations = n(),
      mean_partial_effect = mean(partial_effect),
      mean_within_allocation_variance = mean(partial_effect_se^2),
      between_allocation_sd = sd(partial_effect),
      total_se = sqrt(
        mean_within_allocation_variance +
          between_allocation_sd^2
      ),
      lower_95 = mean_partial_effect - 1.96 * total_se,
      upper_95 = mean_partial_effect + 1.96 * total_se,
      .groups = "drop"
    )

  list(
    per_allocation_results = per_allocation_results,
    per_allocation_curves = per_allocation_curves,
    main_results = ensemble_results,
    main_curves = ensemble_curves
  )
}
