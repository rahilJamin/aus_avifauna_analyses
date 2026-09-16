import::from(mgcv, s)

#' Define the three main LCBD GAMM formulae
#'
#' Returns the Gaussian GAMM specifications used for total, replacement, and
#' richness-difference LCBD. The models share the four focal gradients, effort
#' adjustment, two-dimensional spatial smooth, and reference-pool random effect;
#' only the response changes among components.
#'
#' @return A named list of formulae for total, replacement, and
#'   richness-difference LCBD.
define_lcbd_formulas <- function() {
  shared_terms <- ~
    s(pland_agriculture_2500m, k = 5) +
    s(area_mn_woodland_2500m, k = 5) +
    s(cohesion_woodland_2500m, k = 5) +
    s(miner_detection_rate, k = 5) +
    log_effort +
    s(x_km, y_km, bs = "gp", k = 50) +
    s(reference_pool, bs = "re")

  add_response <- function(response) {
    stats::update.formula(shared_terms, paste(response, "~ ."))
  }

  list(
    total = add_response("relative_lcbd_total"),
    replacement = add_response("relative_lcbd_replacement"),
    richness_difference = add_response("relative_lcbd_richness_difference")
  )
}
