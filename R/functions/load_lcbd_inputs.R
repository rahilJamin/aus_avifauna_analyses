import::from(dplyr, any_of, arrange, filter, left_join, mutate, select)
import::from(magrittr, `%>%`)
import::from(sf, st_coordinates, st_drop_geometry, st_point_on_surface,
             st_read, st_transform)
import::from(tidyr, drop_na)
import::from(tidyselect, all_of)

#' Load and align inputs for the main LCBD model ensemble
#'
#' Reads the prepared LCBD site table and retained species matrix, adds projected
#' kilometre coordinates for reference-pool allocation and spatial smooths,
#' aligns site order with matrix rows, removes Noisy Miner from the response
#' matrix, and drops response species absent after alignment. A reference pool is
#' a local set of nearby sites within which LCBD is calculated, so sites are
#' compared against spatially comparable neighbours.
#'
#' @param site_data_path Path to the candidate LCBD model dataframe.
#' @param species_matrix_path Path to the retained species matrix.
#' @param hex_path Path to the hexagon spatial layer.
#' @param focal_predictors Named character vector of focal predictor columns.
#' @param noisy_miner_species Scientific name to remove from the LCBD response
#'   matrix.
#' @param reference_pool_crs Projected CRS used for LCBD reference-pool
#'   clustering.
#'
#' @return A list with aligned `site_data`, `lcbd_species_matrix`,
#'   `ordered_site_ids`, `lcbd_site_count`, and `lcbd_species_count`.

load_lcbd_inputs <- function(site_data_path,
                             species_matrix_path,
                             hex_path,
                             focal_predictors,
                             noisy_miner_species,
                             reference_pool_crs) {
  # Load the prepared candidate site table and remove columns that are
  # recalculated inside this stage for each allocation.
  site_data <- read_rds_checked(
    site_data_path,
    "site_data",
    required_class = "data.frame",
    required_names = c("hex_id", unname(focal_predictors), "log_effort")
  )
  assert_unique_ids(site_data, "hex_id", "LCBD site table")
  site_data$hex_id <- as.character(site_data$hex_id)

  if (is.na(sf::st_crs(reference_pool_crs))) {
    stop("`reference_pool_crs` is not a valid coordinate reference system.", call. = FALSE)
  }

  site_data <- site_data %>%
    select(-any_of(c(
      "x", "y", "region", "reference_pool",
      "LCBD_tot", "LCBD_repl", "LCBD_rich"
    )))

  # The full retained species matrix is used for LCBD. GLLVM-specific
  # filtering is intentionally not applied here.
  lcbd_species_matrix <- read_rds_checked(
    species_matrix_path,
    "lcbd_species_matrix",
    required_class = "matrix"
  )

  if (is.null(rownames(lcbd_species_matrix)) ||
      is.null(colnames(lcbd_species_matrix)) ||
      anyDuplicated(rownames(lcbd_species_matrix)) ||
      anyDuplicated(colnames(lcbd_species_matrix)) ||
      anyNA(lcbd_species_matrix) ||
      !all(lcbd_species_matrix %in% c(0, 1))) {
    stop(
      "The LCBD species matrix must be uniquely named binary data without missing values.",
      call. = FALSE
    )
  }

  required_site_columns <- c(
    "hex_id",
    unname(focal_predictors),
    "log_effort"
  )
  missing_site_columns <- setdiff(required_site_columns, names(site_data))

  if (length(missing_site_columns) > 0L) {
    stop(
      "The LCBD site table lacks: ",
      paste(missing_site_columns, collapse = ", "),
      call. = FALSE
    )
  }

  # LCBD reference pools use an Australia-wide projected CRS. Coordinates are
  # stored in kilometres because they enter k-means and the spatial smooth.
  hexagons <- st_read(hex_path, quiet = TRUE)
  assert_unique_ids(hexagons, "hex_id", "hexagon grid")

  hexagons <- hexagons %>%
    mutate(hex_id = as.character(hex_id))

  hexagons <- hexagons %>%
    filter(hex_id %in% site_data$hex_id) %>%
    st_transform(reference_pool_crs)

  hexagon_points <- suppressWarnings(st_point_on_surface(hexagons))
  projected_coordinates <- st_coordinates(hexagon_points)

  coordinate_lookup <- hexagon_points %>%
    mutate(
      x_km = projected_coordinates[, 1] / 1000,
      y_km = projected_coordinates[, 2] / 1000
    ) %>%
    st_drop_geometry() %>%
    select(hex_id, x_km, y_km)
  assert_unique_ids(coordinate_lookup, "hex_id", "LCBD coordinate lookup")

  site_data <- site_data %>%
    left_join(coordinate_lookup, by = "hex_id") %>%
    drop_na(all_of(c(unname(focal_predictors), "log_effort", "x_km", "y_km")))

  common_site_ids <- intersect(
    as.character(site_data$hex_id),
    rownames(lcbd_species_matrix)
  )

  site_data <- site_data %>%
    filter(as.character(hex_id) %in% common_site_ids) %>%
    arrange(match(as.character(hex_id), common_site_ids))

  ordered_site_ids <- as.character(site_data$hex_id)
  lcbd_species_matrix <- lcbd_species_matrix[
    ordered_site_ids,
    ,
    drop = FALSE
  ]

  # Noisy Miner detection is a predictor, so the species is removed from the
  # response matrix exactly once before any dissimilarities are calculated.
  lcbd_species_matrix <- lcbd_species_matrix[
    ,
    colnames(lcbd_species_matrix) != noisy_miner_species,
    drop = FALSE
  ]

  if (nrow(lcbd_species_matrix) < 2L || ncol(lcbd_species_matrix) < 1L ||
      any(rowSums(lcbd_species_matrix) == 0)) {
    stop("Insufficient aligned sites or species for LCBD analysis.", call. = FALSE)
  }
  lcbd_species_matrix <- lcbd_species_matrix[
    ,
    colSums(lcbd_species_matrix) > 0,
    drop = FALSE
  ]

  stopifnot(identical(rownames(lcbd_species_matrix), ordered_site_ids))
  stopifnot(all(lcbd_species_matrix %in% c(0, 1)))

  list(
    site_data = site_data,
    lcbd_species_matrix = lcbd_species_matrix,
    ordered_site_ids = ordered_site_ids,
    lcbd_site_count = nrow(lcbd_species_matrix),
    lcbd_species_count = ncol(lcbd_species_matrix)
  )
}
