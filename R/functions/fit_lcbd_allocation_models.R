import::from(dplyr, bind_rows)
import::from(magrittr, `%>%`)
import::from(mgcv, bam)
import::from(stats, gaussian)

#' Fit LCBD GAMMs across reference-pool allocations
#'
#' Fits Gaussian GAMMs for every allocation, focal gradient, and LCBD component.
#' Extracted summaries and curves are retained for all allocations; full fitted
#' model objects are retained only for the first allocation to keep checkpoints
#' manageable.
#'
#' @param allocation_data Nested list of allocation and gradient-specific model
#'   data returned by `prepare_lcbd_reference_allocations()`.
#' @param allocation_seeds Integer seeds used for allocations.
#' @param focal_predictors Named character vector of focal predictor columns.
#' @param focal_labels Named character vector of focal predictor labels.
#' @param formulas Named list of GAMM formulas by LCBD component.
#' @param prediction_grids Named list of focal-gradient prediction grids.
#' @param model_threads Number of threads passed to `mgcv::bam()`.
#'
#' @return A list with per-allocation result rows, per-allocation curve rows,
#'   and representative fitted models from the first allocation.

fit_lcbd_allocation_models <- function(allocation_data,
                                       allocation_seeds,
                                       focal_predictors,
                                       focal_labels,
                                       formulas,
                                       prediction_grids,
                                       model_threads) {
  n_reference_pool_allocations <- length(allocation_seeds)
  component_order <- c("total", "replacement", "richness_difference")

  per_allocation_result_rows <- vector("list", n_reference_pool_allocations)
  per_allocation_curve_rows <- vector("list", n_reference_pool_allocations)
  representative_models <- NULL

  for (allocation in seq_len(n_reference_pool_allocations)) {
    allocation_seed <- allocation_seeds[allocation]
    message(
      "Fitting LCBD allocation ",
      allocation,
      " of ",
      n_reference_pool_allocations,
      " (seed ",
      allocation_seed,
      ")"
    )

    allocation_results <- list()
    allocation_curves <- list()
    fitted_models <- list()

    # Preserve the original fitting order: each gradient is fitted for total,
    # replacement and richness-difference LCBD before moving to the next gradient.
    for (focal_gradient in names(focal_predictors)) {
      focal_data <- allocation_data[[allocation]][[focal_gradient]]
      focal_variable <- focal_predictors[[focal_gradient]]
      focal_label <- focal_labels[[focal_gradient]]
      gradient_values <- prediction_grids[[focal_gradient]]
      fitted_models[[focal_gradient]] <- list()

      for (lcbd_component in component_order) {
        fitted_model <- bam(
          formulas[[lcbd_component]],
          data = focal_data,
          family = gaussian("identity"),
          method = "fREML",
          discrete = TRUE,
          select = TRUE,
          nthreads = model_threads
        )

        fitted_models[[focal_gradient]][[lcbd_component]] <- fitted_model

        extracted <- extract_lcbd_effect(
          fitted_model,
          focal_data,
          focal_variable,
          focal_gradient,
          focal_label,
          lcbd_component,
          gradient_values,
          allocation,
          allocation_seed
        )

        allocation_results[[length(allocation_results) + 1L]] <-
          extracted$result
        allocation_curves[[length(allocation_curves) + 1L]] <-
          extracted$curve
      }
    }

    per_allocation_result_rows[[allocation]] <- bind_rows(allocation_results)
    per_allocation_curve_rows[[allocation]] <- bind_rows(allocation_curves)

    # Keep full fitted objects for the first allocation only. Other allocations
    # are represented by extracted results and curves to avoid huge checkpoints.
    if (allocation == 1L) {
      representative_models <- c(list(seed = allocation_seed), fitted_models)
    }

    rm(fitted_models)
    invisible(gc())
  }

  list(
    per_allocation_results = bind_rows(per_allocation_result_rows),
    per_allocation_curves = bind_rows(per_allocation_curve_rows),
    representative_models = representative_models
  )
}
