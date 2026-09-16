import::from(adespatial, beta.div.comp)

#' Calculate regional Jaccard dissimilarity components
#'
#' Calculates pairwise total Jaccard, replacement, and richness-difference
#' matrices used to compare bird assemblage dissimilarity among IBRA7 subregions
#' at low and high gradient endpoints.
#'
#' @param species_matrix Binary site-by-species matrix with aligned site rows.
#' @param sum_tolerance Maximum allowed absolute difference between total Jaccard
#'   and replacement plus richness difference.
#'
#' @return A list with `component_matrices` and `component_sum_error`.


calculate_regional_jaccard_components <- function(species_matrix,
                                                  sum_tolerance = 1e-10) {
  species_matrix <- as.matrix(species_matrix)
  if (is.null(rownames(species_matrix)) || is.null(colnames(species_matrix))) {
    stop("The regional species matrix must have row and column names.", call. = FALSE)
  }
  if (nrow(species_matrix) < 2L || ncol(species_matrix) < 1L) {
    stop("Regional beta diversity requires at least two sites and one species.", call. = FALSE)
  }
  if (anyDuplicated(rownames(species_matrix)) ||
      anyDuplicated(colnames(species_matrix)) ||
      anyNA(species_matrix) || !is.numeric(species_matrix) ||
      !all(species_matrix %in% c(0, 1))) {
    stop("The regional species matrix must be uniquely named binary data.", call. = FALSE)
  }
  if (any(rowSums(species_matrix) == 0)) {
    stop("Every regional site must contain at least one species.", call. = FALSE)
  }
  if (!is.numeric(sum_tolerance) || length(sum_tolerance) != 1L ||
      is.na(sum_tolerance) || !is.finite(sum_tolerance) || sum_tolerance < 0) {
    stop("`sum_tolerance` must be one non-negative number.", call. = FALSE)
  }

  jaccard_components <- beta.div.comp(
    species_matrix,
    coef = "J",
    quant = FALSE
  )

  component_matrices <- list(
    "Total Jaccard" = as.matrix(jaccard_components$D),
    "Replacement" = as.matrix(jaccard_components$repl),
    "Richness difference" = as.matrix(jaccard_components$rich)
  )

  component_sum_error <- max(abs(
    component_matrices[["Total Jaccard"]] -
      component_matrices[["Replacement"]] -
      component_matrices[["Richness difference"]]
  ))

  if (component_sum_error > sum_tolerance) {
    stop(
      "Replacement and richness difference do not sum to total Jaccard.",
      call. = FALSE
    )
  }

  list(
    component_matrices = component_matrices,
    component_sum_error = component_sum_error
  )
}
