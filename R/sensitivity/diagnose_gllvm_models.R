#' Diagnose one completed GLLVM without fitting it again
#'
#' Both main environmental and fourth-corner models use this function. Each is
#' read in a separate process, avoiding simultaneous multi-gigabyte model loads.
#' @param model_path Completed GLLVM checkpoint.
#' @param model_id Stable identifier used in output names and tables.
#' @param prepared Main GLLVM preparation bundle, for exact response checks.
#' @param config,settings Main and sensitivity settings.
#' @return Tracked compact diagnostic RDS and native diagnostic figure files.
diagnose_gllvm_model <- function(model_path, model_id, prepared, config, settings) {
  fit <- readRDS(model_path)
  assert_same_response(as.matrix(prepared$Y), as.matrix(fit$y))
  set.seed(settings$residual_seed)
  residuals <- align_gllvm_residuals(stats::residuals(fit), as.matrix(fit$y))
  xy <- as.matrix(prepared$gllvm_df[c("x", "y")])
  rownames(xy) <- prepared$gllvm_df$hex_id
  moran <- calculate_residual_moran(residuals, xy, settings$distance_breaks_m,
    settings$gllvm_moran_style, settings$gllvm_moran_alternative)
  moran <- dplyr::mutate(dplyr::group_by(moran, lower_km, upper_km),
                         p_bh = stats::p.adjust(p_value, "BH"))
  moran <- dplyr::ungroup(moran)
  moran$model <- model_id

  # Export all environmental coefficients with their model-based uncertainty;
  # species are never selected because one comparison crosses p = 0.05.
  coef_table <- summary(fit)$Coef.tableX
  coefficients <- data.frame()
  if (!is.null(coef_table)) {
    coefficients <- data.frame(term = rownames(coef_table),
      estimate = coef_table[, "Estimate"], se = coef_table[, "Std. Error"],
      p_value = coef_table[, "Pr(>|z|)"], model = model_id, row.names = NULL)
  }
  ll <- stats::logLik(fit)
  method <- if (model_id == "fourth_corner") "EVA" else "VA"
  n_mem <- if (is.null(fit$X)) NA_integer_ else
    sum(!colnames(fit$X) %in% config$gllvm_base_predictors)
  diagnostic <- data.frame(model = model_id, n_sites = nrow(residuals),
    n_species = ncol(residuals), method = method,
    n_latent_variables = if (model_id == "fourth_corner") config$gllvm_fourth_corner_num_lv else config$gllvm_env_num_lv,
    n_dbmem = n_mem, convergence_code = paste(fit$convergence, collapse = ";"),
    # Some gllvm versions store one convergence code per optimisation stage.
    # Treat the fit as converged only when every reported code is zero.
    converged = length(fit$convergence) > 0L &&
      all(as.character(fit$convergence) == "0"),
    log_likelihood = as.numeric(ll), aic = as.numeric(stats::AIC(fit)),
    n_parameters = if (is.null(attr(ll, "df"))) NA_real_ else as.numeric(attr(ll, "df")),
    residual_mean = mean(residuals), residual_sd = stats::sd(as.numeric(residuals)),
    residual_seed = settings$residual_seed)

  # Preserve the original package-native five-panel diagnostic and Q-Q envelope.
  # It is generated for BOTH main models, and for the matched no-dbMEM model.
  figure <- sensitivity_path(settings$diagnostic_dir, paste0(model_id, "_diagnostics.svg"))
  dir.create(dirname(figure), recursive = TRUE, showWarnings = FALSE)
  grDevices::svg(figure, width = 10, height = 8.5)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(3, 2), mar = c(4, 4, 2.2, 1))
  set.seed(settings$residual_seed)
  graphics::plot(fit, which = 1:5,
    caption = c("Residuals vs linear predictor", "Normal Q-Q", "Residuals vs site index",
                "Residuals vs species index", "Scale-location"),
    var.colors = "grey35", add.smooth = TRUE, envelopes = TRUE,
    reps = settings$qq_envelope_simulations, envelope.col = c("#2F6F73", "#C9DFDE"))
  path <- sensitivity_path(settings$diagnostic_dir, paste0(model_id, "_diagnostics.rds"))
  saveRDS(list(summary = diagnostic, moran = moran, coefficients = coefficients,
               figure = figure), path)
  c(path, figure)
}
