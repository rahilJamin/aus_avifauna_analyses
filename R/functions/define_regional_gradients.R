import::from(stats, quantile)
import::from(tibble, tibble)

#' Define regional convergence endpoint gradients
#'
#' Builds low/high endpoint definitions for the four regional convergence
#' gradients. An endpoint is the low or high end of a habitat/disturbance
#' gradient used for regional comparison. Agriculture uses fixed percentage
#' thresholds; patch area, cohesion, and Noisy Miner detection use study-wide
#' lower and upper quantiles.
#'
#' @param site_data Regional site table with unscaled gradient columns.
#' @param agriculture_low_maximum Maximum agricultural cover for the low
#'   endpoint.
#' @param agriculture_high_minimum Minimum agricultural cover for the high
#'   endpoint.
#' @param lower_tail_probability Lower-tail quantile probability for non-fixed
#'   gradients.
#' @param upper_tail_probability Upper-tail quantile probability for non-fixed
#'   gradients.
#'
#' @return A tibble describing gradient names, predictor columns, thresholds,
#'   endpoint definitions, and plot labels.

define_regional_gradients <- function(site_data,
                                      agriculture_low_maximum,
                                      agriculture_high_minimum,
                                      lower_tail_probability,
                                      upper_tail_probability) {
  required_columns <- c(
    "agricultural_cover",
    "mean_woodland_patch_area",
    "woodland_cohesion",
    "noisy_miner_detection"
  )
  missing_columns <- setdiff(required_columns, names(site_data))
  if (length(missing_columns) > 0L) {
    stop(
      "The regional site table lacks: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  numeric_columns <- vapply(site_data[required_columns], is.numeric, logical(1))
  if (!all(numeric_columns)) {
    stop("Regional gradient columns must be numeric.", call. = FALSE)
  }
  if (anyNA(site_data[required_columns]) ||
      any(!is.finite(as.matrix(site_data[required_columns])))) {
    stop("Regional gradient values must be present and finite.", call. = FALSE)
  }
  if (!is.numeric(lower_tail_probability) ||
      !is.numeric(upper_tail_probability) ||
      length(lower_tail_probability) != 1L ||
      length(upper_tail_probability) != 1L ||
      is.na(lower_tail_probability) || is.na(upper_tail_probability) ||
      !is.finite(lower_tail_probability) || !is.finite(upper_tail_probability) ||
      lower_tail_probability <= 0 || upper_tail_probability >= 1 ||
      lower_tail_probability >= upper_tail_probability) {
    stop("Regional tail probabilities must be ordered values between zero and one.", call. = FALSE)
  }
  if (!is.numeric(agriculture_low_maximum) ||
      !is.numeric(agriculture_high_minimum) ||
      length(agriculture_low_maximum) != 1L ||
      length(agriculture_high_minimum) != 1L ||
      is.na(agriculture_low_maximum) || is.na(agriculture_high_minimum) ||
      !is.finite(agriculture_low_maximum) ||
      !is.finite(agriculture_high_minimum) ||
      agriculture_low_maximum >= agriculture_high_minimum) {
    stop("Agriculture endpoint thresholds must be finite and ordered.", call. = FALSE)
  }

  lower_percentile <- 100 * lower_tail_probability
  upper_percentile <- 100 * upper_tail_probability
  low_tail_share <- 100 * lower_tail_probability
  high_tail_share <- 100 * (1 - upper_tail_probability)

  patch_area_low_threshold <- quantile(
    site_data$mean_woodland_patch_area,
    lower_tail_probability,
    na.rm = TRUE,
    names = FALSE
  )

  patch_area_high_threshold <- quantile(
    site_data$mean_woodland_patch_area,
    upper_tail_probability,
    na.rm = TRUE,
    names = FALSE
  )

  cohesion_low_threshold <- quantile(
    site_data$woodland_cohesion,
    lower_tail_probability,
    na.rm = TRUE,
    names = FALSE
  )

  cohesion_high_threshold <- quantile(
    site_data$woodland_cohesion,
    upper_tail_probability,
    na.rm = TRUE,
    names = FALSE
  )

  noisy_miner_low_threshold <- quantile(
    site_data$noisy_miner_detection,
    lower_tail_probability,
    na.rm = TRUE,
    names = FALSE
  )

  noisy_miner_high_threshold <- quantile(
    site_data$noisy_miner_detection,
    upper_tail_probability,
    na.rm = TRUE,
    names = FALSE
  )

  gradient_definitions <- tibble(
    gradient = c("agriculture", "patch_area", "cohesion", "noisy_miner"),
    gradient_label = c(
      "Agricultural cover",
      "Mean woodland patch area",
      "Woodland cohesion",
      "Noisy Miner detection rate"
    ),
    predictor = c(
      "agricultural_cover",
      "mean_woodland_patch_area",
      "woodland_cohesion",
      "noisy_miner_detection"
    ),
    low_threshold = c(
      agriculture_low_maximum,
      patch_area_low_threshold,
      cohesion_low_threshold,
      noisy_miner_low_threshold
    ),
    high_threshold = c(
      agriculture_high_minimum,
      patch_area_high_threshold,
      cohesion_high_threshold,
      noisy_miner_high_threshold
    ),
    threshold_rule = c(
      "Fixed percentage of land cover",
      rep(
        sprintf(
          "Study-wide %.0fth and %.0fth percentiles",
          lower_percentile,
          upper_percentile
        ),
        3
      )
    ),
    low_endpoint_definition = c(
      sprintf("Agricultural cover < %g%%", agriculture_low_maximum),
      sprintf("Lowest %g%% of mean woodland patch area", low_tail_share),
      sprintf("Lowest %g%% of woodland cohesion", low_tail_share),
      sprintf("Lowest %g%% of Noisy Miner detection rate", low_tail_share)
    ),
    high_endpoint_definition = c(
      sprintf("Agricultural cover > %g%%", agriculture_high_minimum),
      sprintf("Highest %g%% of mean woodland patch area", high_tail_share),
      sprintf("Highest %g%% of woodland cohesion", high_tail_share),
      sprintf("Highest %g%% of Noisy Miner detection rate", high_tail_share)
    ),
    low_axis_label = c(
      sprintf("Low end\n(<%g%% agriculture)", agriculture_low_maximum),
      sprintf("Low end\n(lowest %g%%)", low_tail_share),
      sprintf("Low end\n(lowest %g%%)", low_tail_share),
      sprintf("Low end\n(lowest %g%%)", low_tail_share)
    ),
    high_axis_label = c(
      sprintf("High end\n(>%g%% agriculture)", agriculture_high_minimum),
      sprintf("High end\n(highest %g%%)", high_tail_share),
      sprintf("High end\n(highest %g%%)", high_tail_share),
      sprintf("High end\n(highest %g%%)", high_tail_share)
    )
  )

  if (any(gradient_definitions$low_threshold >=
          gradient_definitions$high_threshold)) {
    print(gradient_definitions)
    stop(
      "At least one gradient does not have distinct low and high endpoints.",
      call. = FALSE
    )
  }

  gradient_definitions
}
