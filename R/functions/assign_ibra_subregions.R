#' Assign model sites to IBRA7 subregions
#'
#' Converts the projected model coordinates to points and attaches the IBRA7
#' subregion containing each point. The join must return exactly one row for
#' every input site so duplicated polygon matches cannot silently duplicate
#' observations before GLLVM fitting.
#'
#' @param model_df Site-level model data containing unique `hex_id`, `x`, and
#'   `y` columns. Coordinates must use `model_coord_crs`.
#' @param ibra_path Path to the main IBRA7 shapefile.
#' @param model_coord_crs Coordinate reference system of `x` and `y`.
#'
#' @return `model_df` with an `ibra_sub` column.
assign_ibra_subregions <- function(model_df, ibra_path, model_coord_crs) {
  required_columns <- c("hex_id", "x", "y")
  missing_columns <- setdiff(required_columns, names(model_df))

  if (length(missing_columns) > 0L) {
    stop(
      "The model site table lacks: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  assert_unique_ids(model_df, "hex_id", "GLLVM site table")

  if (!all(vapply(model_df[c("x", "y")], is.numeric, logical(1))) ||
      anyNA(model_df[c("x", "y")]) ||
      any(!is.finite(as.matrix(model_df[c("x", "y")])))) {
    stop("Model coordinates must be numeric, present, and finite.", call. = FALSE)
  }

  if (is.na(sf::st_crs(model_coord_crs))) {
    stop("`model_coord_crs` is not a valid coordinate reference system.", call. = FALSE)
  }

  ibra <- sf::st_read(ibra_path, quiet = TRUE) |>
    sf::st_transform(crs = model_coord_crs)

  if (!"SUB_NAME_7" %in% names(ibra)) {
    stop("The IBRA spatial layer lacks the SUB_NAME_7 column.", call. = FALSE)
  }

  model_points <- sf::st_as_sf(
    model_df,
    coords = c("x", "y"),
    crs = model_coord_crs
  )

  joined_sites <- sf::st_join(
    model_points,
    ibra["SUB_NAME_7"],
    left = TRUE
  )

  if (nrow(joined_sites) != nrow(model_df) ||
      !identical(
        as.character(joined_sites$hex_id),
        as.character(model_df$hex_id)
      )) {
    stop(
      "The IBRA join did not return exactly one subregion result per model site.",
      call. = FALSE
    )
  }

  model_df$ibra_sub <- as.character(joined_sites$SUB_NAME_7)
  model_df
}
