#' Calculate residual Moran's I in non-overlapping distance bands
#'
#' Empty bands and constant residual series are reported explicitly. They are
#' never interpreted as absence of spatial autocorrelation. Tests are analytic
#' randomisation tests, matching the original scripts, not permutation p-values.
#' @param residual_matrix Named site-by-series matrix (species or LCBD responses).
#' @param coordinates Two-column projected coordinates in metres, in site order.
#' @param distance_breaks_m Strictly increasing non-negative band boundaries.
#' @param style Spatial weights style, B for LCBD or W for GLLVM.
#' @param alternative Test alternative used in the original diagnostic.
#' @return One row per series and band; signed I, expectation, p and support.
calculate_residual_moran <- function(residual_matrix, coordinates, distance_breaks_m,
                                     style, alternative) {
  residual_matrix <- as.matrix(residual_matrix)
  coordinates <- as.matrix(coordinates)
  if (nrow(residual_matrix) != nrow(coordinates) || ncol(coordinates) != 2L ||
      is.null(rownames(residual_matrix)) || is.null(colnames(residual_matrix)) ||
      !identical(rownames(residual_matrix), rownames(coordinates)) ||
      any(!is.finite(residual_matrix)) || any(!is.finite(coordinates)) ||
      anyDuplicated(as.data.frame(coordinates))) {
    stop("Spatial residuals and unique projected coordinates must be finite and aligned.", call. = FALSE)
  }
  if (length(distance_breaks_m) < 2L || any(!is.finite(distance_breaks_m)) ||
      min(distance_breaks_m) < 0 || any(diff(distance_breaks_m) <= 0)) {
    stop("Distance boundaries must be finite, non-negative and strictly increasing.", call. = FALSE)
  }
  rows <- list()
  for (b in seq_len(length(distance_breaks_m) - 1L)) {
    nb <- spdep::dnearneigh(coordinates, distance_breaks_m[b], distance_breaks_m[b + 1L],
                           bounds = c("GT", "LE"), longlat = FALSE)
    counts <- spdep::card(nb)
    weights <- if (sum(counts)) spdep::nb2listw(nb, style = style, zero.policy = TRUE) else NULL
    for (j in seq_len(ncol(residual_matrix))) {
      x <- residual_matrix[, j]
      status <- if (is.null(weights)) "no_pairs" else if (sum(counts > 0) < 4L) {
        "too_few_connected_sites"
      } else if (stats::sd(x) == 0) "constant_residuals" else "ok"
      result <- NULL
      if (status == "ok") {
        result <- spdep::moran.test(x, weights, randomisation = TRUE,
          alternative = alternative, zero.policy = TRUE, adjust.n = TRUE)
        if (!is.finite(result$p.value) || any(!is.finite(result$estimate))) status <- "undefined_test"
      }
      rows[[length(rows) + 1L]] <- data.frame(
        series = colnames(residual_matrix)[j], lower_km = distance_breaks_m[b] / 1000,
        upper_km = distance_breaks_m[b + 1L] / 1000,
        midpoint_km = mean(distance_breaks_m[c(b, b + 1L)]) / 1000,
        n_pairs = sum(counts) / 2, n_connected_sites = sum(counts > 0),
        status = status,
        moran_i = if (is.null(result)) NA_real_ else unname(result$estimate[1]),
        expected_i = if (is.null(result)) NA_real_ else unname(result$estimate[2]),
        p_value = if (status != "ok") NA_real_ else result$p.value
      )
    }
  }
  dplyr::bind_rows(rows)
}
