#' Calculate Podani-Jaccard beta-diversity matrices
#'
#' Calculates total Jaccard dissimilarity, replacement, and richness-difference
#' components from a binary site-by-species matrix. These components describe
#' how bird assemblages differ between pairs of sites. Site-level mean
#' replacement and richness-difference values exclude the diagonal
#' self-comparison.
#'
#' @param species_matrix Binary site-by-species matrix with site IDs as row names
#'   and species names as column names.
#'
#' @return A list with `dist_mats`, `site_mean_replacement`, and
#'   `site_mean_richness_difference`.

calculate_lcbd_beta_matrices <- function(species_matrix) {
  species_matrix <- as.matrix(species_matrix)

  if (is.null(rownames(species_matrix))) {
    stop("The species matrix must have row names containing hex IDs.", call. = FALSE)
  }

  if (is.null(colnames(species_matrix))) {
    stop("The species matrix must have species column names.", call. = FALSE)
  }

  if (anyDuplicated(rownames(species_matrix))) {
    stop("The species matrix contains duplicated hex IDs.", call. = FALSE)
  }

  if (anyDuplicated(colnames(species_matrix))) {
    stop("The species matrix contains duplicated species names.", call. = FALSE)
  }

  if (anyNA(species_matrix) || !is.numeric(species_matrix) ||
      !all(is.finite(species_matrix)) ||
      !all(species_matrix %in% c(0, 1))) {
    stop("The species matrix must contain only presence-absence values.", call. = FALSE)
  }

  species_matrix <- species_matrix[, colSums(species_matrix) > 0, drop = FALSE]

  if (nrow(species_matrix) < 2L) {
    stop("At least two sites are required to calculate beta diversity.", call. = FALSE)
  }

  if (ncol(species_matrix) < 1L) {
    stop("At least one retained species is required to calculate beta diversity.", call. = FALSE)
  }

  if (any(rowSums(species_matrix) == 0)) {
    stop("Every retained site must contain at least one species.", call. = FALSE)
  }

  shared_species <- species_matrix %*% t(species_matrix)
  site_richness <- diag(shared_species)
  not_shared <- abs(sweep(shared_species, 2, site_richness))
  sum_not_shared <- not_shared + t(not_shared)
  minimum_not_shared <- pmin(not_shared, t(not_shared))
  maximum_not_shared <- pmax(not_shared, t(not_shared))
  jaccard_denominator <- shared_species + sum_not_shared

  if (any(jaccard_denominator == 0)) {
    stop("At least one site pair has an empty species union.", call. = FALSE)
  }

  total <- sum_not_shared / jaccard_denominator
  replacement <- (2 * minimum_not_shared) / jaccard_denominator
  richness_difference <- abs(maximum_not_shared - minimum_not_shared) /
    jaccard_denominator

  diag(total) <- 0
  diag(replacement) <- 0
  diag(richness_difference) <- 0

  dimnames(total) <- dimnames(shared_species)
  dimnames(replacement) <- dimnames(shared_species)
  dimnames(richness_difference) <- dimnames(shared_species)

  replacement_for_means <- replacement
  richness_difference_for_means <- richness_difference
  diag(replacement_for_means) <- NA_real_
  diag(richness_difference_for_means) <- NA_real_

  list(
    dist_mats = list(
      total = total,
      replacement = replacement,
      richness_difference = richness_difference
    ),
    site_mean_replacement = rowMeans(replacement_for_means, na.rm = TRUE),
    site_mean_richness_difference = rowMeans(
      richness_difference_for_means,
      na.rm = TRUE
    )
  )
}
