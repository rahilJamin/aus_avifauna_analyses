import::from(tibble, tibble)

#' Validate saved LCBD beta matrices against rebuilt matrices
#'
#' Compares the beta-diversity matrices saved by the preparation stage with
#' matrices rebuilt from the aligned LCBD response matrix. This is a refactor
#' guard against stale or misaligned checkpoints; model fitting still uses the
#' rebuilt matrices.
#'
#' @param saved_beta_path Path to the saved beta-matrix bundle.
#' @param rebuilt_distance_matrices List of rebuilt distance matrices with
#'   `total`, `replacement`, and `richness_difference` components.
#' @param ordered_site_ids Site IDs in the order used by the LCBD model stage.
#' @param tolerance Maximum allowed absolute difference.
#'
#' @return A tibble reporting the maximum absolute difference for each LCBD
#'   component.

validate_lcbd_beta_matrices <- function(saved_beta_path,
                                        rebuilt_distance_matrices,
                                        ordered_site_ids,
                                        tolerance = 1e-12) {
  # The saved beta bundle is a checkpoint from the preparation stage. It is
  # never used instead of the matrices rebuilt from the aligned response matrix.
  saved_beta_object <- read_rds_checked(
    saved_beta_path,
    "saved_beta_matrices",
    required_class = "list"
  )

  compare_lcbd_beta_matrices(
    saved_beta_object = saved_beta_object,
    rebuilt_distance_matrices = rebuilt_distance_matrices,
    ordered_site_ids = ordered_site_ids,
    tolerance = tolerance
  )
}

#' Compare saved and rebuilt LCBD beta matrices
#'
#' Performs the in-memory structure, identifier, finiteness, and numerical
#' checks used by `validate_lcbd_beta_matrices()`. Keeping comparison separate
#' from file reading lets the scientific validation be tested directly.
#'
#' @param saved_beta_object Saved beta-matrix bundle containing `dist_mats`.
#' @inheritParams validate_lcbd_beta_matrices
#'
#' @return A tibble reporting the maximum absolute difference for each LCBD
#'   component.
compare_lcbd_beta_matrices <- function(saved_beta_object,
                                       rebuilt_distance_matrices,
                                       ordered_site_ids,
                                       tolerance = 1e-12) {
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      is.na(tolerance) || !is.finite(tolerance) || tolerance < 0) {
    stop("`tolerance` must be one non-negative finite number.", call. = FALSE)
  }

  if (anyNA(ordered_site_ids) || any(!nzchar(as.character(ordered_site_ids))) ||
      anyDuplicated(ordered_site_ids)) {
    stop("Ordered LCBD site IDs must be unique and non-missing.", call. = FALSE)
  }
  ordered_site_ids <- as.character(ordered_site_ids)

  if (!"dist_mats" %in% names(saved_beta_object)) {
    stop("The saved beta object does not contain dist_mats.", call. = FALSE)
  }

  saved_distance_matrices <- saved_beta_object$dist_mats
  required_components <- c("total", "replacement", "richness_difference")

  if (!all(required_components %in% names(saved_distance_matrices))) {
    stop("One or more saved dissimilarity components are missing.", call. = FALSE)
  }

  if (!all(required_components %in% names(rebuilt_distance_matrices))) {
    stop("One or more rebuilt dissimilarity components are missing.", call. = FALSE)
  }

  saved_ids_available <- vapply(
    saved_distance_matrices[required_components],
    function(saved_matrix) {
      is.matrix(saved_matrix) &&
        nrow(saved_matrix) == ncol(saved_matrix) &&
        !is.null(rownames(saved_matrix)) &&
        !is.null(colnames(saved_matrix)) &&
        !anyDuplicated(rownames(saved_matrix)) &&
        !anyDuplicated(colnames(saved_matrix)) &&
        identical(rownames(saved_matrix), colnames(saved_matrix)) &&
        all(ordered_site_ids %in% rownames(saved_matrix)) &&
        all(ordered_site_ids %in% colnames(saved_matrix)) &&
        !anyNA(saved_matrix) && all(is.finite(saved_matrix)) &&
        isTRUE(all.equal(saved_matrix, t(saved_matrix), tolerance = tolerance))
    },
    logical(1)
  )

  if (!all(saved_ids_available)) {
    stop(
      "The saved dissimilarities do not cover the full LCBD site set.",
      call. = FALSE
    )
  }


  rebuilt_valid <- vapply(
    rebuilt_distance_matrices[required_components],
    function(rebuilt_matrix) {
      is.matrix(rebuilt_matrix) &&
        nrow(rebuilt_matrix) == ncol(rebuilt_matrix) &&
        identical(rownames(rebuilt_matrix), ordered_site_ids) &&
        identical(colnames(rebuilt_matrix), ordered_site_ids) &&
        !anyNA(rebuilt_matrix) && all(is.finite(rebuilt_matrix)) &&
        isTRUE(all.equal(
          rebuilt_matrix,
          t(rebuilt_matrix),
          tolerance = tolerance
        ))
    },
    logical(1)
  )
  if (!all(rebuilt_valid)) {
    stop("The rebuilt dissimilarities are malformed or misaligned.", call. = FALSE)
  }

  distance_validation <- tibble(
    lcbd_component = required_components,
    maximum_absolute_difference = c(
      max(abs(
        saved_distance_matrices$total[ordered_site_ids, ordered_site_ids] -
          rebuilt_distance_matrices$total
      )),
      max(abs(
        saved_distance_matrices$replacement[
          ordered_site_ids,
          ordered_site_ids
        ] - rebuilt_distance_matrices$replacement
      )),
      max(abs(
        saved_distance_matrices$richness_difference[
          ordered_site_ids,
          ordered_site_ids
        ] - rebuilt_distance_matrices$richness_difference
      ))
    )
  )

  if (any(distance_validation$maximum_absolute_difference > tolerance)) {
    print(distance_validation)
    stop(
      "Rebuilt and saved LCBD dissimilarities are not identical.",
      call. = FALSE
    )
  }

  distance_validation
}
