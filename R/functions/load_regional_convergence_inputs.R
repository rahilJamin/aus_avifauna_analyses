import::from(dplyr, arrange, filter, left_join, mutate, select,
             transmute)
import::from(magrittr, `%>%`)
import::from(sf, st_drop_geometry, st_join, st_point_on_surface, st_read,
             st_transform)
import::from(tidyr, drop_na)

#' Load and align inputs for regional convergence analysis
#'
#' Reads the candidate site table, unscaled landscape metrics, retained species
#' matrix, hexagon grid, and IBRA7 subregions. IBRA7 subregions are Australian
#' bioregional units used here to compare regional assemblage convergence. The
#' function removes Noisy Miner from the response matrix, attaches each site to
#' an IBRA7 subregion, and aligns species-matrix rows to the site table.
#'
#' @param site_data_path Path to the candidate LCBD site table.
#' @param landscape_data_path Path to processed landscape metrics.
#' @param species_matrix_path Path to the retained species matrix.
#' @param hex_path Path to the hexagon spatial layer.
#' @param ibra_path Path to the IBRA7 subregion spatial layer.
#' @param noisy_miner_species Scientific name to remove from the response matrix.
#' @param regional_crs Projected CRS used for spatial joins.
#'
#' @return A list with aligned `site_data` and `species_matrix`.

load_regional_convergence_inputs <- function(site_data_path,
                                             landscape_data_path,
                                             species_matrix_path,
                                             hex_path,
                                             ibra_path,
                                             noisy_miner_species,
                                             regional_crs) {
  # The candidate-site table supplies Noisy Miner detection. The landscape RDS
  # supplies unscaled values so endpoint definitions remain interpretable.
  site_data <- read_rds_checked(
    site_data_path,
    "site_data",
    required_class = "data.frame",
    required_names = c("hex_id", "miner_detection_rate")
  )
  assert_unique_ids(site_data, "hex_id", "regional site table")

  site_data <- site_data %>%
    transmute(
      hex_id = as.character(hex_id),
      noisy_miner_detection = miner_detection_rate
    )

  landscape_data <- read_rds_checked(
    landscape_data_path,
    "landscape_data",
    required_class = "data.frame",
    required_names = c(
      "hex_id", "pland_agriculture_2500m", "area_mn_woodland_2500m",
      "cohesion_woodland_2500m"
    )
  )
  assert_unique_ids(landscape_data, "hex_id", "regional landscape table")

  landscape_data <- landscape_data %>%
    transmute(
      hex_id = as.character(hex_id),
      agricultural_cover = pland_agriculture_2500m,
      mean_woodland_patch_area = area_mn_woodland_2500m,
      woodland_cohesion = cohesion_woodland_2500m
    )

  site_data <- site_data %>%
    left_join(landscape_data, by = "hex_id") %>%
    drop_na(
      agricultural_cover,
      mean_woodland_patch_area,
      woodland_cohesion,
      noisy_miner_detection
    )

  regional_predictors <- c(
    "agricultural_cover",
    "mean_woodland_patch_area",
    "woodland_cohesion",
    "noisy_miner_detection"
  )
  if (nrow(site_data) < 2L ||
      !all(vapply(site_data[regional_predictors], is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(site_data[regional_predictors])))) {
    stop("Regional predictors must be numeric and finite for at least two sites.", call. = FALSE)
  }

  if (min(site_data$agricultural_cover) < 0 ||
      max(site_data$agricultural_cover) > 100) {
    stop("Agricultural cover must range from 0 to 100%.", call. = FALSE)
  }

  if (min(site_data$noisy_miner_detection) < 0 ||
      max(site_data$noisy_miner_detection) > 1) {
    stop("Noisy Miner detection rate must range from 0 to 1.", call. = FALSE)
  }

  species_matrix <- as.matrix(read_rds_checked(
    species_matrix_path,
    "species_matrix",
    required_class = "matrix"
  ))

  if (is.null(rownames(species_matrix)) || is.null(colnames(species_matrix))) {
    stop("The species matrix must have row and column names.", call. = FALSE)
  }

  if (anyDuplicated(rownames(species_matrix)) ||
      anyDuplicated(colnames(species_matrix)) ||
      anyNA(species_matrix) || !is.numeric(species_matrix) ||
      !all(species_matrix %in% c(0, 1))) {
    stop("The species matrix is not presence-absence data.", call. = FALSE)
  }

  # Noisy Miner is removed from the response matrix because its detection rate
  # defines one of the four explanatory gradients.
  if (noisy_miner_species %in% colnames(species_matrix)) {
    species_matrix <- species_matrix[
      ,
      colnames(species_matrix) != noisy_miner_species,
      drop = FALSE
    ]
  }

  species_matrix <- species_matrix[
    ,
    colSums(species_matrix) > 0,
    drop = FALSE
  ]

  # Attach IBRA7 subregion names to the retained hexagons using point-on-surface
  # joins so every hexagon is represented by one spatial lookup point.
  if (is.na(sf::st_crs(regional_crs))) {
    stop("`regional_crs` is not a valid coordinate reference system.", call. = FALSE)
  }

  ibra <- st_read(ibra_path, quiet = TRUE) %>%
    st_transform(regional_crs)

  if (!"SUB_NAME_7" %in% names(ibra)) {
    stop("The IBRA layer does not contain SUB_NAME_7.", call. = FALSE)
  }

  hexagons <- st_read(hex_path, quiet = TRUE)
  assert_unique_ids(hexagons, "hex_id", "hexagon grid")

  hexagons <- hexagons %>%
    mutate(hex_id = as.character(hex_id))

  hexagons <- hexagons %>%
    filter(hex_id %in% site_data$hex_id) %>%
    st_transform(regional_crs)

  subregion_lookup <- st_join(
    suppressWarnings(st_point_on_surface(hexagons)),
    select(ibra, SUB_NAME_7),
    left = TRUE
  ) %>%
    st_drop_geometry() %>%
    transmute(
      hex_id,
      ibra_subregion = as.character(SUB_NAME_7)
    )
  assert_unique_ids(subregion_lookup, "hex_id", "IBRA subregion lookup")

  site_data <- site_data %>%
    left_join(subregion_lookup, by = "hex_id") %>%
    drop_na(ibra_subregion)

  # Keep only shared sites and put the species matrix in exactly the same order
  # as the site table used by the regional analysis.
  common_site_ids <- intersect(site_data$hex_id, rownames(species_matrix))

  site_data <- site_data %>%
    filter(hex_id %in% common_site_ids) %>%
    arrange(match(hex_id, common_site_ids))

  species_matrix <- species_matrix[
    site_data$hex_id,
    ,
    drop = FALSE
  ]

  species_matrix <- species_matrix[
    ,
    colSums(species_matrix) > 0,
    drop = FALSE
  ]

  if (nrow(species_matrix) < 2L || ncol(species_matrix) < 1L) {
    stop("Insufficient aligned sites or species for regional analysis.", call. = FALSE)
  }
  if (!identical(rownames(species_matrix), site_data$hex_id) ||
      any(rowSums(species_matrix) == 0)) {
    stop("Regional site data and species matrix are not validly aligned.", call. = FALSE)
  }

  site_data$site_row <- seq_len(nrow(site_data))

  list(
    site_data = site_data,
    species_matrix = species_matrix
  )
}
