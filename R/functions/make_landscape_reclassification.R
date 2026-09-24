#' Translate raster category tables into analysis classes
#'
#' The raster Value field is the pixel code. CLUM groups those codes using
#' SIMPN, whereas NVIS rules refer directly to Value. Keeping every row of the
#' source category table leaves unmapped rows as NA; terra's `others` fallback
#' applies only to codes absent from that table.
#'
#' @param categories Raster category data frame, with Value and, for CLUM, SIMPN.
#' @param group_column Column whose codes are named in `groups`.
#' @param groups Named numeric vector mapping source groups to output classes;
#'   NA values explicitly remove a group such as water.
#' @return Two-column numeric matrix (`from`, `to`) for terra::classify().
make_landscape_reclassification <- function(categories, group_column, groups) {
  # GDAL can round-trip the same pixel-code field as `Value` or `value`.
  # Match field names case-insensitively without changing the category codes.
  value_column <- names(categories)[tolower(names(categories)) == "value"]
  group_field <- names(categories)[tolower(names(categories)) == tolower(group_column)]
  if (!is.data.frame(categories) || nrow(categories) == 0L ||
      length(value_column) != 1L || length(group_field) != 1L) {
    stop("Raster category table must contain Value and ", group_column,
         " and at least one row.", call. = FALSE)
  }
  values <- categories[[value_column]]
  if (!is.numeric(values) || anyNA(values) ||
      any(!is.finite(values)) || anyDuplicated(values)) {
    stop("Raster category Value codes must be unique, finite numbers.", call. = FALSE)
  }
  if (!is.numeric(groups) || is.null(names(groups)) ||
      anyNA(names(groups)) || any(!nzchar(names(groups))) ||
      anyDuplicated(names(groups)) || any(!is.finite(groups[!is.na(groups)]))) {
    stop("Classification groups must be numeric with unique, non-empty names.",
         call. = FALSE)
  }

  # Unmatched category-table groups intentionally map to NA, not to `others`.
  output_class <- unname(groups[as.character(categories[[group_field]])])
  cbind(from = values, to = output_class)
}
