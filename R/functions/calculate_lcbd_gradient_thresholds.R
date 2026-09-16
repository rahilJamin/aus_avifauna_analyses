import::from(dplyr, group_by, summarise, select)
import::from(magrittr, `%>%`)
import::from(stats, quantile)
import::from(tidyr, pivot_longer)
import::from(tidyselect, all_of, everything)

#' Calculate global LCBD gradient-tail thresholds
#'
#' Calculates low and high quantiles for each focal predictor across the aligned
#' LCBD site table. These thresholds remain fixed across reference-pool
#' allocations so only pool assignment changes across seeds.
#'
#' @param site_data LCBD site table containing focal predictor columns.
#' @param focal_predictors Named character vector of focal predictor columns.
#' @param lower_tail_probability Lower-tail quantile probability.
#' @param upper_tail_probability Upper-tail quantile probability.
#'
#' @return A tibble with focal predictor names and global low/high thresholds.


calculate_lcbd_gradient_thresholds <- function(site_data,
                                               focal_predictors,
                                               lower_tail_probability,
                                               upper_tail_probability) {
  if (!is.character(focal_predictors) || length(focal_predictors) < 1L ||
      is.null(names(focal_predictors)) || anyNA(focal_predictors) ||
      any(!nzchar(focal_predictors)) || any(!nzchar(names(focal_predictors))) ||
      anyDuplicated(focal_predictors) || anyDuplicated(names(focal_predictors))) {
    stop(
      "Focal predictors must be a uniquely named character vector of unique columns.",
      call. = FALSE
    )
  }
  if (!is.numeric(lower_tail_probability) ||
      !is.numeric(upper_tail_probability) ||
      length(lower_tail_probability) != 1L ||
      length(upper_tail_probability) != 1L ||
      is.na(lower_tail_probability) || is.na(upper_tail_probability) ||
      !is.finite(lower_tail_probability) ||
      !is.finite(upper_tail_probability) ||
      lower_tail_probability <= 0 || upper_tail_probability >= 1 ||
      lower_tail_probability >= upper_tail_probability) {
    stop("LCBD tail probabilities must be ordered values between zero and one.", call. = FALSE)
  }

  missing_columns <- setdiff(unname(focal_predictors), names(site_data))
  if (length(missing_columns) > 0L) {
    stop(
      "The LCBD site table lacks focal predictors: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  predictor_columns <- site_data[unname(focal_predictors)]
  numeric_columns <- vapply(predictor_columns, is.numeric, logical(1))
  if (!all(numeric_columns) || anyNA(predictor_columns) ||
      any(!is.finite(as.matrix(predictor_columns)))) {
    stop("LCBD focal predictors must be numeric, present, and finite.", call. = FALSE)
  }

  # These thresholds are fixed across all allocations so that a different seed
  # changes only the geographic pool assignment, not the gradient definition.
  thresholds <- site_data %>%
    select(all_of(unname(focal_predictors))) %>%
    pivot_longer(
      cols = everything(),
      names_to = "focal_predictor",
      values_to = "predictor_value"
    ) %>%
    group_by(focal_predictor) %>%
    summarise(
      global_low_threshold = as.numeric(quantile(
        predictor_value,
        lower_tail_probability,
        names = FALSE
      )),
      global_high_threshold = as.numeric(quantile(
        predictor_value,
        upper_tail_probability,
        names = FALSE
      )),
      .groups = "drop"
    )

  if (any(thresholds$global_low_threshold >= thresholds$global_high_threshold)) {
    stop("At least one LCBD focal predictor has no distinct low and high tail.", call. = FALSE)
  }

  thresholds
}
