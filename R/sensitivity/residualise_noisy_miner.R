#' Residualise Noisy Miner detection using the original count model
#'
#' Preserves 05a_noisy_miner_residualisation.R: quasibinomial logit GAM,
#' landscape smooths, log effort, spatial GP and broad IBRA random intercept.
#' @param site_data Completed main LCBD candidate sites, before pool filtering.
#' @param hex_path,ibra_path Grid and IBRA shapefiles.
#' @param config,settings Main and sensitivity settings.
#' @return Fitted count model, aligned deviance residuals and association checks.
residualise_noisy_miner <- function(site_data, hex_path, ibra_path, config, settings) {
  validate_miner_counts(site_data)
  assert_unique_ids(site_data, "hex_id", "Noisy Miner residualisation sites")

  # Join on representative points in Australian Albers, preserving every site.
  # An unmapped coastal site receives the original explicit outside-region level.
  grid <- sf::st_read(hex_path, quiet = TRUE)
  grid <- grid[as.character(grid$hex_id) %in% site_data$hex_id, "hex_id"]
  grid <- sf::st_transform(grid, settings$miner_region_crs)
  ibra <- sf::st_transform(sf::st_read(ibra_path, quiet = TRUE), settings$miner_region_crs)
  if (!settings$miner_region_field %in% names(ibra)) stop("IBRA broad-region field is missing.", call. = FALSE)
  points <- sf::st_point_on_surface(sf::st_geometry(grid))
  lookup <- sf::st_drop_geometry(sf::st_join(
    sf::st_sf(hex_id = as.character(grid$hex_id), geometry = points),
    ibra[settings$miner_region_field], left = TRUE
  ))
  # Several subregions can share a broad region; remove only identical pairs.
  # Conflicting region assignments require attention instead of arbitrary choice.
  lookup <- unique(lookup)
  assert_unique_ids(lookup, "hex_id", "Noisy Miner broad-region lookup")
  idx <- match(site_data$hex_id, lookup$hex_id)
  if (anyNA(idx)) stop("The grid does not cover all main LCBD sites.", call. = FALSE)
  region <- as.character(lookup[[settings$miner_region_field]][idx])
  n_unmapped <- sum(is.na(region))
  region[is.na(region)] <- "Outside mapped IBRA polygon"
  site_data$ibra_region <- factor(region)

  # WORKFLOW CHANGE: the current main pipeline defines a checklist as a sampled
  # month-year within a grid cell. Both numerator and denominator use that same
  # definition; the former date/resource proxy is not reintroduced here.
  # Counts admit exact zero/one rates and quasibinomial dispersion accommodates
  # extra-binomial variation. No squeeze or Gaussian rate model is used.
  s <- mgcv::s
  model <- mgcv::bam(
    cbind(n_miner_checklists, n_checklists - n_miner_checklists) ~
      s(pland_agriculture_2500m, k = 5) +
      s(area_mn_woodland_2500m, k = 5) +
      s(cohesion_woodland_2500m, k = 5) + log_effort +
      s(x_km, y_km, bs = "gp", k = 50) + s(ibra_region, bs = "re"),
    family = stats::quasibinomial(link = "logit"), data = site_data,
    method = "fREML", discrete = TRUE, select = TRUE, nthreads = config$model_threads
  )
  site_data$miner_residuals <- as.numeric(stats::residuals(model, type = "deviance"))
  if (any(!is.finite(site_data$miner_residuals))) stop("Non-finite Noisy Miner residuals.", call. = FALSE)
  predictors <- unname(config$lcbd_focal_predictors[names(config$lcbd_focal_predictors) != "noisy_miner"])
  association <- dplyr::bind_rows(lapply(predictors, function(variable) data.frame(
    predictor = variable,
    raw_spearman_rho = stats::cor(site_data$miner_detection_rate, site_data[[variable]], method = "spearman"),
    residual_spearman_rho = stats::cor(site_data$miner_residuals, site_data[[variable]], method = "spearman")
  )))
  list(analysis_site_data = site_data, model = model, association = association,
       diagnostics = data.frame(n_sites = nrow(site_data),
         n_zero_detections = sum(site_data$n_miner_checklists == 0),
         n_all_detections = sum(site_data$n_miner_checklists == site_data$n_checklists),
         n_unmapped_ibra = n_unmapped, converged = isTRUE(model$converged),
         dispersion = summary(model)$dispersion, deviance_explained = summary(model)$dev.expl))
}
