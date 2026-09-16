import::from(dplyr, arrange, count, filter, inner_join, mutate, n, transmute)
import::from(magrittr, `%>%`)
import::from(stats, kmeans)
import::from(tibble, as_tibble)

#' Assign LCBD sites to spatial reference pools
#'
#' Runs one k-means allocation on projected kilometre coordinates, converts
#' arbitrary k-means cluster labels to canonical reference-pool labels, identifies
#' complete pools, and records allocation membership for diagnostics. A reference
#' pool is a group of nearby sites used as the local comparison set for LCBD.
#'
#' @param site_data LCBD site table with `hex_id`, `x_km`, and `y_km` columns.
#' @param allocation Allocation index.
#' @param allocation_seed Random seed used for this allocation.
#' @param reference_pool_k Number of spatial reference pools.
#' @param reference_pool_nstart Number of random starts passed to `kmeans()`.
#' @param minimum_pool_sites Minimum sites required for a pool to be complete.
#'
#' @return A list with the k-means object, site data with pool labels, pool-size
#'   summary, complete-pool site data, and allocation membership table.

assign_lcbd_reference_pools <- function(site_data,
                                        allocation,
                                        allocation_seed,
                                        reference_pool_k,
                                        reference_pool_nstart,
                                        minimum_pool_sites) {
  required_columns <- c("hex_id", "x_km", "y_km")
  missing_columns <- setdiff(required_columns, names(site_data))
  if (length(missing_columns) > 0L) {
    stop(
      "The LCBD site table lacks: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  assert_unique_ids(site_data, "hex_id", "LCBD site table")
  if (!is.numeric(site_data$x_km) || !is.numeric(site_data$y_km) ||
      anyNA(site_data[c("x_km", "y_km")]) ||
      any(!is.finite(as.matrix(site_data[c("x_km", "y_km")]))) ) {
    stop("LCBD reference-pool coordinates must be numeric and finite.", call. = FALSE)
  }

  valid_positive_whole <- function(value) {
    is.numeric(value) && length(value) == 1L && !is.na(value) &&
      is.finite(value) && value >= 1 && value == floor(value)
  }
  if (!valid_positive_whole(allocation) ||
      !valid_positive_whole(reference_pool_k) ||
      !valid_positive_whole(reference_pool_nstart) ||
      !valid_positive_whole(minimum_pool_sites) ||
      !is.numeric(allocation_seed) || length(allocation_seed) != 1L ||
      is.na(allocation_seed) || !is.finite(allocation_seed)) {
    stop("LCBD allocation and pool settings are invalid.", call. = FALSE)
  }
  coordinates <- unique(site_data[c("x_km", "y_km")])
  if (reference_pool_k > nrow(site_data) ||
      reference_pool_k > nrow(coordinates)) {
    stop(
      "The number of reference pools exceeds the number of sites or distinct coordinates.",
      call. = FALSE
    )
  }

  site_data$hex_id <- as.character(site_data$hex_id)
  set.seed(allocation_seed)

  # K-means is run on projected kilometre coordinates. The spatial reference
  # system is chosen upstream in the configuration so that clustering is done
  # consistently across the whole analysis.
  geographical_kmeans <- kmeans(
    dplyr::select(site_data, x_km, y_km),
    centers = reference_pool_k,
    nstart = reference_pool_nstart
  )

  # K-means labels are arbitrary. Ordering pool centroids from west/south to
  # east/north makes labels reproducible enough to compare partitions.
  centre_lookup <- as_tibble(
    geographical_kmeans$centers,
    rownames = "original_pool"
  ) %>%
    mutate(original_pool = as.integer(original_pool)) %>%
    arrange(x_km, y_km) %>%
    mutate(reference_pool = seq_len(n()))

  canonical_pool <- centre_lookup$reference_pool[
    match(geographical_kmeans$cluster, centre_lookup$original_pool)
  ]

  allocation_site_data <- site_data %>%
    mutate(reference_pool = factor(canonical_pool))

  pool_sizes <- allocation_site_data %>%
    count(reference_pool, name = "reference_pool_size") %>%
    mutate(complete_pool = reference_pool_size >= minimum_pool_sites)

  if (!any(pool_sizes$complete_pool)) {
    stop("No reference pool meets the configured minimum site count.", call. = FALSE)
  }

  complete_pool_data <- allocation_site_data %>%
    inner_join(filter(pool_sizes, complete_pool), by = "reference_pool") %>%
    arrange(reference_pool, hex_id) %>%
    mutate(reference_pool = droplevels(factor(reference_pool)))

  allocation_membership <- allocation_site_data %>%
    transmute(
      allocation,
      seed = allocation_seed,
      hex_id,
      reference_pool = as.character(reference_pool),
      complete_pool = as.character(reference_pool) %in%
        as.character(filter(pool_sizes, complete_pool)$reference_pool)
    )

  list(
    geographical_kmeans = geographical_kmeans,
    allocation_site_data = allocation_site_data,
    pool_sizes = pool_sizes,
    complete_pool_data = complete_pool_data,
    allocation_membership = allocation_membership
  )
}
