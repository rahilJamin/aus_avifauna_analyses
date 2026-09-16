#' Validate an identifier column
#'
#' Checks that a workflow table contains the requested identifier column and
#' that every identifier is present and occurs exactly once. Silent duplicate
#' removal can misalign site tables and species matrices, so duplicated or
#' missing identifiers are treated as data errors.
#'
#' @param data A data frame, tibble, or `sf` object.
#' @param id_column Name of the identifier column.
#' @param label Human-readable name used in error messages.
#'
#' @return `data`, invisibly, when the identifier is valid.
assert_unique_ids <- function(data, id_column = "hex_id", label = "data") {
  if (!id_column %in% names(data)) {
    stop(
      "The ", label, " lacks the required `", id_column, "` column.",
      call. = FALSE
    )
  }

  ids <- as.character(data[[id_column]])

  if (anyNA(ids) || any(!nzchar(ids))) {
    stop(
      "The ", label, " contains missing or empty `", id_column, "` values.",
      call. = FALSE
    )
  }

  duplicated_ids <- unique(ids[duplicated(ids)])
  if (length(duplicated_ids) > 0L) {
    preview <- paste(utils::head(duplicated_ids, 5L), collapse = ", ")
    stop(
      "The ", label, " contains duplicated `", id_column, "` values: ",
      preview,
      if (length(duplicated_ids) > 5L) " ..." else "",
      call. = FALSE
    )
  }

  invisible(data)
}
