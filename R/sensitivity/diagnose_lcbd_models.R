#' Compare LCBD spatial control on the saved first allocation
#'
#' This diagnostic fits twelve matched models without the GP smooth. It is an
#' allocation-specific model check, not a reduced ensemble used for inference.
#' Pool random intercepts and the Gaussian main model specification are retained.
#' @param main Main LCBD checkpoint with representative_models.
#' @param config,settings Main and sensitivity settings.
#' @return Model summaries, focal effects and signed residual correlograms.
diagnose_lcbd_models <- function(main, config, settings) {
  s <- mgcv::s
  total <- relative_lcbd_total ~ s(pland_agriculture_2500m, k = 5) +
    s(area_mn_woodland_2500m, k = 5) + s(cohesion_woodland_2500m, k = 5) +
    s(miner_detection_rate, k = 5) + log_effort + s(reference_pool, bs = "re")
  formulas <- list(total = total,
    replacement = stats::update(total, relative_lcbd_replacement ~ .),
    richness_difference = stats::update(total, relative_lcbd_richness_difference ~ .))
  checks <- effects <- correlograms <- list()
  for (gradient in names(config$lcbd_focal_predictors)) {
    for (component in names(formulas)) {
      reference <- main$representative_models[[gradient]][[component]]
      dat <- reference$model
      if (!all(c("x_km", "y_km") %in% names(dat))) stop("Representative LCBD model lacks coordinates.", call. = FALSE)
      reduced <- mgcv::bam(formulas[[component]], data = dat,
        family = stats::gaussian("identity"), method = "fREML", discrete = TRUE,
        select = TRUE, nthreads = config$model_threads)
      # Fit from the saved model frame: same rows, responses and factor levels.
      if (!identical(rownames(reference$model), rownames(reduced$model)) ||
          !identical(stats::model.response(reference$model),
                     stats::model.response(reduced$model))) {
        stop("The LCBD spatial comparison changed response values or row order.", call. = FALSE)
      }
      grids <- main$curves_by_allocation
      grid <- sort(unique(grids$gradient_value[grids$focal_gradient == gradient]))
      for (model_id in c("with_gp", "without_gp")) {
        fit <- if (model_id == "with_gp") reference else reduced
        # `bam()` stores model-specific prediction metadata in its own model
        # frame. In particular, the no-GP fit has a different smooth structure
        # from the saved main fit. Passing the main fit's frame to the reduced
        # fit can therefore fail inside `predict.bam()` even though the rows and
        # predictor values are identical.
        fit_data <- fit$model
        if (!identical(rownames(fit_data), rownames(dat))) {
          stop("LCBD diagnostic model rows are not aligned with their coordinates.",
               call. = FALSE)
        }
        residuals <- matrix(as.numeric(stats::residuals(fit, type = "pearson")), ncol = 1,
          dimnames = list(rownames(fit_data), component))
        xy <- as.matrix(dat[c("x_km", "y_km")]) * 1000
        rownames(xy) <- rownames(residuals)
        moran <- calculate_residual_moran(residuals, xy, settings$distance_breaks_m,
          settings$lcbd_moran_style, settings$lcbd_moran_alternative)
        moran$p_bh <- stats::p.adjust(moran$p_value, "BH")
        moran$model <- model_id
        moran$focal_gradient <- gradient
        correlograms[[length(correlograms) + 1L]] <- moran
        extracted <- extract_lcbd_effect(fit, fit_data,
          config$lcbd_focal_predictors[[gradient]],
          gradient, config$lcbd_focal_labels[[gradient]], component, grid, 1L,
          settings$allocation_seeds[1])$result
        extracted$model <- model_id
        effects[[length(effects) + 1L]] <- extracted
        checks[[length(checks) + 1L]] <- data.frame(focal_gradient = gradient,
          component = component, model = model_id, n_sites = nrow(fit_data),
          converged = isTRUE(fit$converged), deviance_explained = summary(fit)$dev.expl,
          residual_sd = stats::sd(as.numeric(residuals)))
      }
    }
  }
  list(summary = dplyr::bind_rows(checks), effects = dplyr::bind_rows(effects),
       moran = dplyr::bind_rows(correlograms))
}
