#' Validate paired analyses and isolated sensitivity exports
#'
#' These checks prevent changed response data or misordered site identifiers from
#' being mistaken for a change caused by a modelling choice.
#' @param reference,candidate Named site-by-species response matrices.
#' @return Invisibly returns TRUE when response values AND names match exactly.
assert_same_response <- function(reference, candidate) {
  valid <- function(x) is.matrix(x) && is.numeric(x) &&
    !is.null(rownames(x)) && !is.null(colnames(x)) &&
    !anyDuplicated(rownames(x)) && !anyDuplicated(colnames(x)) &&
    !anyNA(x) && all(x %in% c(0, 1))
  if (!valid(reference) || !valid(candidate) ||
      !identical(dimnames(reference), dimnames(candidate)) ||
      !identical(dim(reference), dim(candidate)) ||
      any(reference != candidate)) {
    stop("Paired models must use identical response values, sites and species in the same order.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Check the Noisy Miner numerator, denominator and reported detection rate
#' @param site_data One row per grid cell with monthly checklist counts and rate.
#' @return The input data, invisibly; zero and one detection rates are valid.
validate_miner_counts <- function(site_data) {
  columns <- c("n_checklists", "n_miner_checklists", "miner_detection_rate")
  if (!all(columns %in% names(site_data)) || !nrow(site_data) ||
      !all(vapply(site_data[columns], is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(site_data[columns])))) {
    stop("Monthly Noisy Miner counts and rates must be present and finite.", call. = FALSE)
  }
  n <- site_data$n_checklists
  y <- site_data$n_miner_checklists
  if (any(n <= 0 | y < 0 | y > n | n != round(n) | y != round(y)) ||
      any(abs(site_data$miner_detection_rate - y / n) > 1e-10)) {
    stop("Noisy Miner counts must be integers with 0 <= detections <= checklists and rate = detections/checklists.",
         call. = FALSE)
  }
  invisible(site_data)
}

#' Resolve an export below the dedicated sensitivity directory
#' @param directory Sensitivity subdirectory.
#' @param filename A single filename, without directory components.
#' @return A path; does not create a directory or write a file.
sensitivity_path <- function(directory, filename) {
  path <- normalizePath(file.path(directory, filename), winslash = "/", mustWork = FALSE)
  root <- normalizePath(file.path("outputs", "sensitivity"), winslash = "/", mustWork = FALSE)
  # Normalise dot segments even when the output does not yet exist.
  if (grepl("(^|[/\\\\])\\.\\.([/\\\\]|$)", directory) ||
      grepl("[/\\\\]", filename) || !nzchar(filename) ||
      !startsWith(tolower(path), paste0(tolower(root), "/"))) {
    stop("Sensitivity exports must stay inside outputs/sensitivity.", call. = FALSE)
  }
  file.path(directory, filename)
}

#' Select higher-completeness sites from the completed main matrix
#' @param chao Chao2 table, including unsuccessful estimates.
#' @param species_matrix Main retained incidence matrix.
#' @param threshold Alternative completeness threshold.
#' @param config Main settings, used without modification.
#' @return Character identifiers in the original matrix row order.
sensitivity_site_ids <- function(chao, species_matrix, threshold, config) {
  if (length(threshold) != 1L || !is.finite(threshold) ||
      threshold <= config$min_completeness || threshold > 1) {
    stop("Alternative completeness must exceed the main threshold and be at most one.", call. = FALSE)
  }
  retained <- retain_sampled_hexagons(
    chao, config$min_sampling_units, config$min_species, threshold
  )
  ids <- as.character(retained$hex_id)
  if (!all(ids %in% rownames(species_matrix))) {
    stop("Higher-completeness sites are absent from the main matrix; regenerate consistent main inputs first.",
         call. = FALSE)
  }
  rownames(species_matrix)[rownames(species_matrix) %in% ids]
}

#' Align a GLLVM residual matrix without silently relabelling mismatched rows
#' @param residuals Matrix, or the list returned by residuals.gllvm().
#' @param response Named response matrix stored in the fitted model.
#' @return A finite site-by-species residual matrix in response order.
align_gllvm_residuals <- function(residuals, response) {
  if (is.list(residuals)) residuals <- residuals$residuals
  if (is.null(residuals)) stop("The model returned no residual matrix.", call. = FALSE)
  out <- as.matrix(residuals)
  row_labels <- rownames(out)
  column_labels <- colnames(out)
  response_row_labels <- rownames(response)
  response_column_labels <- colnames(response)
  # Test names before dimensions: a square response matrix is compatible with
  # both orientations numerically, but its species/site labels still identify
  # which orientation gllvm returned.
  transposed_names <- !is.null(row_labels) && !is.null(column_labels) &&
    identical(row_labels, response_column_labels) &&
    identical(column_labels, response_row_labels)
  if (identical(dim(out), rev(dim(response))) && transposed_names) {
    out <- t(out)
  }
  if (!identical(dim(out), dim(response)) || !is.numeric(out) || any(!is.finite(out))) {
    stop("GLLVM residuals have invalid dimensions or non-finite values.", call. = FALSE)
  }
  for (axis in 1:2) {
    labels <- dimnames(out)[[axis]]
    if (!is.null(labels) && !identical(labels, dimnames(response)[[axis]])) {
      stop("GLLVM residual names do not match the response ordering.", call. = FALSE)
    }
  }
  dimnames(out) <- dimnames(response)
  out
}
