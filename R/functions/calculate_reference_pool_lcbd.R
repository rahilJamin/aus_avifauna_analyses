import::from(adespatial, LCBD.comp)
import::from(dplyr, bind_rows, filter, left_join, pull, select)
import::from(magrittr, `%>%`)
import::from(tibble, tibble)

#' Calculate LCBD within complete reference pools
#'
#' Calculates raw LCBD values separately within each complete reference pool for
#' total Jaccard, replacement, and richness-difference components. LCBD is a
#' site-level measure of how compositionally unique a site is within its local
#' reference pool. Relative LCBD values are raw LCBD multiplied by the number of
#' sites in the reference pool.
#'
#' @param complete_pool_data Site data for pools meeting the minimum-size rule.
#' @param distance_matrices List with `total`, `replacement`, and
#'   `richness_difference` matrices.
#' @param allocation Allocation index.
#' @param allocation_seed Seed used for the allocation.
#' @param negative_tolerance Threshold for counting meaningful negative LCBD
#'   values.
#' @param lcbd_sum_tolerance Maximum allowed deviation of raw LCBD sums from one
#'   within each pool/component.
#'
#' @return A list with `allocation_lcbd_data` and `pool_checks`.

calculate_reference_pool_lcbd <- function(complete_pool_data,
                                          distance_matrices,
                                          allocation,
                                          allocation_seed,
                                          negative_tolerance,
                                          lcbd_sum_tolerance) {
  if (nrow(complete_pool_data) == 0L) {
    stop("No complete LCBD reference pools are available.", call. = FALSE)
  }
  if (!"reference_pool" %in% names(complete_pool_data) ||
      !"hex_id" %in% names(complete_pool_data)) {
    stop("Complete-pool data must contain `hex_id` and `reference_pool`.", call. = FALSE)
  }
  assert_unique_ids(complete_pool_data, "hex_id", "complete LCBD pool data")

  valid_tolerance <- function(value) {
    is.numeric(value) && length(value) == 1L && !is.na(value) &&
      is.finite(value) && value >= 0
  }
  if (!valid_tolerance(negative_tolerance) ||
      !valid_tolerance(lcbd_sum_tolerance)) {
    stop("LCBD numerical tolerances must be non-negative finite numbers.", call. = FALSE)
  }

  required_components <- c("total", "replacement", "richness_difference")
  if (!all(required_components %in% names(distance_matrices))) {
    stop("The LCBD distance bundle is missing one or more components.", call. = FALSE)
  }

  complete_pool_data$hex_id <- as.character(complete_pool_data$hex_id)
  complete_pool_data$reference_pool <- droplevels(factor(
    complete_pool_data$reference_pool
  ))

  raw_lcbd_rows <- vector(
    "list",
    nlevels(complete_pool_data$reference_pool)
  )
  pool_check_rows <- vector(
    "list",
    nlevels(complete_pool_data$reference_pool)
  )

  # LCBD is a within-pool quantity. A site therefore receives a different raw
  # and relative LCBD value depending on the sites that share its reference pool
  # in this allocation.
  for (pool_index in seq_along(levels(complete_pool_data$reference_pool))) {
    pool_name <- levels(complete_pool_data$reference_pool)[pool_index]

    pool_site_ids <- complete_pool_data %>%
      filter(reference_pool == pool_name) %>%
      pull(hex_id) %>%
      as.character()

    if (length(pool_site_ids) < 2L) {
      stop("Each complete LCBD reference pool must contain at least two sites.", call. = FALSE)
    }

    matrices_cover_pool <- vapply(
      distance_matrices[required_components],
      function(distance_matrix) {
        is.matrix(distance_matrix) &&
          nrow(distance_matrix) == ncol(distance_matrix) &&
          !is.null(rownames(distance_matrix)) &&
          !is.null(colnames(distance_matrix)) &&
          !anyDuplicated(rownames(distance_matrix)) &&
          !anyDuplicated(colnames(distance_matrix)) &&
          identical(rownames(distance_matrix), colnames(distance_matrix)) &&
          !anyNA(distance_matrix) && all(is.finite(distance_matrix)) &&
          isTRUE(all.equal(
            distance_matrix,
            t(distance_matrix),
            tolerance = lcbd_sum_tolerance
          )) &&
          all(pool_site_ids %in% rownames(distance_matrix)) &&
          all(pool_site_ids %in% colnames(distance_matrix))
      },
      logical(1)
    )
    if (!all(matrices_cover_pool)) {
      stop(
        "LCBD distance matrices are malformed or do not cover a complete pool.",
        call. = FALSE
      )
    }

    pool_total <- LCBD.comp(
      distance_matrices$total[pool_site_ids, pool_site_ids, drop = FALSE],
      sqrt.D = TRUE
    )$LCBD

    pool_replacement <- LCBD.comp(
      distance_matrices$replacement[
        pool_site_ids,
        pool_site_ids,
        drop = FALSE
      ],
      sqrt.D = TRUE
    )$LCBD

    pool_richness_difference <- LCBD.comp(
      distance_matrices$richness_difference[
        pool_site_ids,
        pool_site_ids,
        drop = FALSE
      ],
      sqrt.D = TRUE
    )$LCBD

    pool_size <- length(pool_site_ids)

    raw_lcbd_rows[[pool_index]] <- tibble(
      hex_id = pool_site_ids,
      reference_pool = pool_name,
      reference_pool_size = pool_size,
      raw_lcbd_total = as.numeric(pool_total),
      raw_lcbd_replacement = as.numeric(pool_replacement),
      raw_lcbd_richness_difference = as.numeric(pool_richness_difference),
      relative_lcbd_total = as.numeric(pool_total) * pool_size,
      relative_lcbd_replacement = as.numeric(pool_replacement) * pool_size,
      relative_lcbd_richness_difference =
        as.numeric(pool_richness_difference) * pool_size
    )

    pool_check_rows[[pool_index]] <- tibble(
      allocation,
      seed = allocation_seed,
      reference_pool = pool_name,
      reference_pool_size = pool_size,
      lcbd_component = c("total", "replacement", "richness_difference"),
      raw_lcbd_sum = c(
        sum(pool_total),
        sum(pool_replacement),
        sum(pool_richness_difference)
      ),
      minimum_raw_lcbd = c(
        min(pool_total),
        min(pool_replacement),
        min(pool_richness_difference)
      ),
      n_values_below_zero = c(
        sum(pool_total < 0),
        sum(pool_replacement < 0),
        sum(pool_richness_difference < 0)
      ),
      n_values_below_negative_tolerance = c(
        sum(pool_total < -negative_tolerance),
        sum(pool_replacement < -negative_tolerance),
        sum(pool_richness_difference < -negative_tolerance)
      )
    )
  }

  raw_lcbd <- bind_rows(raw_lcbd_rows)
  pool_checks <- bind_rows(pool_check_rows)

  validate_lcbd_pool_checks(
    pool_checks = pool_checks,
    allocation = allocation,
    negative_tolerance = negative_tolerance,
    lcbd_sum_tolerance = lcbd_sum_tolerance
  )

  allocation_lcbd_data <- complete_pool_data %>%
    select(-reference_pool_size) %>%
    left_join(raw_lcbd, by = c("hex_id", "reference_pool"))

  stopifnot(nrow(allocation_lcbd_data) == nrow(complete_pool_data))
  stopifnot(!anyNA(allocation_lcbd_data$relative_lcbd_total))
  stopifnot(!anyNA(allocation_lcbd_data$relative_lcbd_replacement))
  stopifnot(!anyNA(
    allocation_lcbd_data$relative_lcbd_richness_difference
  ))

  list(
    allocation_lcbd_data = allocation_lcbd_data,
    pool_checks = pool_checks
  )
}

#' Validate within-pool LCBD numerical checks
#'
#' Stops before modelling when raw LCBD values fail to sum to one. Negative
#' component LCBD values are retained and counted in `pool_checks` as a
#' diagnostic because component dissimilarity matrices can yield negative site
#' contributions; this matches the original analysis workflow.
#'
#' @param pool_checks Pool/component checks created by
#'   `calculate_reference_pool_lcbd()`.
#' @param allocation Allocation index used in error messages.
#' @param negative_tolerance Threshold used when negative LCBD values were
#'   counted for the diagnostic columns in `pool_checks`.
#' @param lcbd_sum_tolerance Maximum accepted deviation from an LCBD sum of one.
#'
#' @return `pool_checks`, invisibly, when all checks pass.
validate_lcbd_pool_checks <- function(pool_checks,
                                      allocation,
                                      negative_tolerance,
                                      lcbd_sum_tolerance) {
  valid_tolerance <- function(value) {
    is.numeric(value) && length(value) == 1L && !is.na(value) &&
      is.finite(value) && value >= 0
  }
  if (!valid_tolerance(negative_tolerance) ||
      !valid_tolerance(lcbd_sum_tolerance)) {
    stop("LCBD numerical tolerances must be non-negative finite numbers.", call. = FALSE)
  }

  required_columns <- c(
    "raw_lcbd_sum",
    "n_values_below_negative_tolerance"
  )
  if (!all(required_columns %in% names(pool_checks)) || nrow(pool_checks) == 0L) {
    stop("LCBD pool checks are empty or incomplete.", call. = FALSE)
  }

  if (anyNA(pool_checks[required_columns])) {
    stop("LCBD pool checks contain missing values.", call. = FALSE)
  }

  invalid_sums <- abs(pool_checks$raw_lcbd_sum - 1) > lcbd_sum_tolerance
  if (any(invalid_sums)) {
    print(pool_checks[invalid_sums, , drop = FALSE])
    stop(
      "Raw LCBD did not sum to one in allocation ", allocation, ".",
      call. = FALSE
    )
  }

  invisible(pool_checks)
}
