import::from(dplyr, first, last, n_distinct)
import::from(magrittr, `%>%`)
import::from(stats, coef, median, predict, vcov)
import::from(tibble, tibble)
import::from(mgcv, bam)

#' Extract a focal LCBD smooth from a fitted GAMM
#'
#' Builds a prediction grid for one focal gradient, extracts the partial smooth
#' effect and uncertainty, and calculates the endpoint contrast using the model's
#' linear-predictor matrix so endpoint covariance is retained.
#'
#' @param fitted_model A fitted `mgcv::bam()` model.
#' @param fitted_data Data used to fit the model.
#' @param focal_variable Column name of the focal predictor.
#' @param focal_gradient Stable gradient identifier.
#' @param focal_label Readable gradient label.
#' @param lcbd_component LCBD component identifier.
#' @param gradient_values Numeric prediction grid for the focal predictor.
#' @param allocation Reference-pool allocation index.
#' @param allocation_seed Seed used for the allocation.
#'
#' @return A list with `curve`, the predicted partial-effect curve, and
#'   `result`, the endpoint contrast and model summary row.

extract_lcbd_effect <- function(
    fitted_model,
    fitted_data,
    focal_variable,
    focal_gradient,
    focal_label,
    lcbd_component,
    gradient_values,
    allocation,
    allocation_seed) {

  prediction_data <- fitted_data[
    rep(1L, length(gradient_values)),
    ,
    drop = FALSE
  ]

  prediction_data$pland_agriculture_2500m <- median(
    fitted_data$pland_agriculture_2500m,
    na.rm = TRUE
  )
  prediction_data$area_mn_woodland_2500m <- median(
    fitted_data$area_mn_woodland_2500m,
    na.rm = TRUE
  )
  prediction_data$cohesion_woodland_2500m <- median(
    fitted_data$cohesion_woodland_2500m,
    na.rm = TRUE
  )
  prediction_data$miner_detection_rate <- median(
    fitted_data$miner_detection_rate,
    na.rm = TRUE
  )
  prediction_data$log_effort <- median(
    fitted_data$log_effort,
    na.rm = TRUE
  )
  prediction_data$x_km <- median(fitted_data$x_km, na.rm = TRUE)
  prediction_data$y_km <- median(fitted_data$y_km, na.rm = TRUE)
  prediction_data$reference_pool <- factor(
    levels(fitted_data$reference_pool)[1],
    levels = levels(fitted_data$reference_pool)
  )
  prediction_data[[focal_variable]] <- gradient_values

  smooth_term <- paste0("s(", focal_variable, ")")

  partial_prediction <- predict(
    fitted_model,
    newdata = prediction_data,
    type = "terms",
    terms = smooth_term,
    se.fit = TRUE,
    unconditional = TRUE
  )

  partial_fit <- as.numeric(partial_prediction$fit)
  partial_se <- as.numeric(partial_prediction$se.fit)

  # The linear-predictor matrix retains the covariance between the two
  # endpoints, giving a valid uncertainty estimate for their difference.
  linear_predictor_matrix <- predict(
    fitted_model,
    newdata = prediction_data,
    type = "lpmatrix"
  )
  focal_coefficient_columns <- startsWith(
    colnames(linear_predictor_matrix),
    paste0(smooth_term, ".")
  )

  if (!any(focal_coefficient_columns)) {
    stop("No coefficients found for focal smooth: ", smooth_term, call. = FALSE)
  }

  endpoint_contrast <-
    linear_predictor_matrix[nrow(linear_predictor_matrix), ] -
    linear_predictor_matrix[1, ]
  endpoint_contrast[!focal_coefficient_columns] <- 0

  endpoint_difference <- sum(
    endpoint_contrast * coef(fitted_model)
  )
  endpoint_variance <- as.numeric(
    t(endpoint_contrast) %*%
      vcov(fitted_model, unconditional = TRUE) %*%
      endpoint_contrast
  )
  endpoint_se <- sqrt(max(endpoint_variance, 0))

  fitted_summary <- summary(fitted_model)
  smooth_table <- fitted_summary$s.table

  if (!smooth_term %in% rownames(smooth_table)) {
    stop("Focal smooth not found: ", smooth_term, call. = FALSE)
  }

  curve <- tibble(
    allocation,
    seed = allocation_seed,
    focal_gradient,
    focal_label,
    focal_variable,
    lcbd_component,
    gradient_value = gradient_values,
    partial_effect = partial_fit,
    partial_effect_se = partial_se
  )

  result <- tibble(
    allocation,
    seed = allocation_seed,
    focal_gradient,
    focal_label,
    focal_variable,
    lcbd_component,
    n_sites = nrow(fitted_data),
    n_informative_pools = n_distinct(fitted_data$reference_pool),
    endpoint_low = first(gradient_values),
    endpoint_high = last(gradient_values),
    endpoint_difference,
    endpoint_se,
    edf = smooth_table[smooth_term, "edf"],
    test_statistic = smooth_table[
      smooth_term,
      ncol(smooth_table) - 1L
    ],
    p_value = smooth_table[smooth_term, ncol(smooth_table)],
    deviance_explained = fitted_summary$dev.expl,
    model_converged = isTRUE(fitted_model$converged)
  )

  list(curve = curve, result = result)
}
