import::from(dplyr, filter)
import::from(magrittr, `%>%`)

#' Make LCBD prediction grids
#'
#' Creates one evenly spaced prediction grid for each focal gradient using the
#' shared central range available across all reference-pool allocations.
#'
#' @param common_curve_limits Data frame with `focal_gradient`, `curve_lower`,
#'   and `curve_upper` columns.
#' @param curve_grid_points Number of grid points per focal gradient.
#'
#' @return A named list of numeric prediction grids.

make_lcbd_prediction_grids <- function(common_curve_limits,
                                       curve_grid_points) {
  required_columns <- c("focal_gradient", "curve_lower", "curve_upper")
  missing_columns <- setdiff(required_columns, names(common_curve_limits))
  if (length(missing_columns) > 0L) {
    stop(
      "The shared LCBD curve-limit table lacks: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(common_curve_limits) < 1L ||
      anyNA(common_curve_limits[required_columns]) ||
      any(!nzchar(as.character(common_curve_limits$focal_gradient))) ||
      anyDuplicated(common_curve_limits$focal_gradient) ||
      !is.numeric(common_curve_limits$curve_lower) ||
      !is.numeric(common_curve_limits$curve_upper) ||
      any(!is.finite(common_curve_limits$curve_lower)) ||
      any(!is.finite(common_curve_limits$curve_upper)) ||
      any(common_curve_limits$curve_lower >= common_curve_limits$curve_upper)) {
    stop("Shared LCBD curve limits must be unique, finite, and ordered.", call. = FALSE)
  }
  if (!is.numeric(curve_grid_points) || length(curve_grid_points) != 1L ||
      is.na(curve_grid_points) || !is.finite(curve_grid_points) ||
      curve_grid_points < 2L || curve_grid_points != floor(curve_grid_points)) {
    stop("`curve_grid_points` must be a whole number of at least two.", call. = FALSE)
  }

  # Build one shared prediction grid per focal gradient using the central range
  # that is available in every allocation.
  grid_list <- lapply(common_curve_limits$focal_gradient, function(gradient) {
    limits <- common_curve_limits %>%
      filter(focal_gradient == gradient)

    seq(
      limits$curve_lower,
      limits$curve_upper,
      length.out = curve_grid_points
    )
  })

  names(grid_list) <- common_curve_limits$focal_gradient
  grid_list
}
