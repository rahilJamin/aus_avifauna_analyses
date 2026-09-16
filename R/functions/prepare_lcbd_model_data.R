import::from(magrittr, `%>%`)
import::from(dplyr, across, arrange, filter, left_join, mutate,
             n_distinct, select)
import::from(sf, st_coordinates, st_drop_geometry, st_point_on_surface, st_read,
             st_transform)
import::from(tibble, tibble)
import::from(tidyr, drop_na)
import::from(tidyselect, all_of, where)

#' Prepare site data and beta matrices for LCBD modelling
#'
#' Aligns the retained species matrix, Chao2/effort metadata, processed
#' landscape metrics, and 2.5 km grid-cell coordinates. The function applies the
#' urban cover filter, scales numeric landscape predictors after filtering,
#' removes Noisy Miner from the LCBD response matrix, and calculates
#' Podani-Jaccard beta matrices for downstream checks and summaries.
#'
#' @param species_matrix Binary site-by-species matrix with 2.5 km grid-cell
#'   identifiers (`hex_id`) as row names.
#' @param hex_metadata Hexagon metadata from the iNEXT matrix stage.
#' @param landscape_data Processed 2.5 km landscape metrics.
#' @param hex_path Path to the hexagon spatial layer.
#' @param noisy_miner_species Scientific name to remove from the LCBD response
#'   matrix.
#' @param urban_cover_max Maximum permitted urban cover percentage.
#' @param model_coord_crs Projected CRS used to derive model coordinates.
#'
#' @return A list containing `candidate_model_df`, `model_df`,
#'   `lcbd_species_matrix`, `beta_matrices`, and `attrition`.

prepare_lcbd_model_data <- function(species_matrix,
                                    hex_metadata,
                                    landscape_data,
                                    hex_path,
                                    noisy_miner_species,
                                    urban_cover_max,
                                    model_coord_crs) {
  if (!is.numeric(urban_cover_max) || length(urban_cover_max) != 1L ||
      is.na(urban_cover_max) || !is.finite(urban_cover_max)) {
    stop("`urban_cover_max` must be one finite number.", call. = FALSE)
  }
  if (is.na(sf::st_crs(model_coord_crs))) {
    stop("`model_coord_crs` is not a valid coordinate reference system.", call. = FALSE)
  }

  species_matrix <- as.matrix(species_matrix)

  if (is.null(rownames(species_matrix)) || is.null(colnames(species_matrix))) {
    stop("The species matrix must have site row names and species column names.", call. = FALSE)
  }
  if (anyDuplicated(rownames(species_matrix)) ||
      anyDuplicated(colnames(species_matrix))) {
    stop("The species matrix contains duplicated site or species names.", call. = FALSE)
  }
  if (anyNA(species_matrix) || !all(species_matrix %in% c(0, 1))) {
    stop("The species matrix must contain only zeroes and ones.", call. = FALSE)
  }

  required_metadata_columns <- c(
    "hex_id",
    "S_obs",
    "sampling_units",
    "n_checklists",
    "miner_detection_rate"
  )
  missing_metadata_columns <- setdiff(required_metadata_columns, names(hex_metadata))

  if (length(missing_metadata_columns) > 0L) {
    stop(
      "The hex metadata table lacks: ",
      paste(missing_metadata_columns, collapse = ", "),
      call. = FALSE
    )
  }
  assert_unique_ids(hex_metadata, "hex_id", "hex metadata table")

  required_landscape_columns <- c(
    "hex_id",
    "pland_urban_2500m",
    "pland_agriculture_2500m",
    "area_mn_woodland_2500m",
    "cohesion_woodland_2500m"
  )
  missing_landscape_columns <- setdiff(
    required_landscape_columns,
    names(landscape_data)
  )

  if (length(missing_landscape_columns) > 0L) {
    stop(
      "The landscape metrics table lacks: ",
      paste(missing_landscape_columns, collapse = ", "),
      call. = FALSE
    )
  }
  assert_unique_ids(landscape_data, "hex_id", "landscape metrics table")

  landscape_numeric_columns <- setdiff(required_landscape_columns, "hex_id")
  if (!all(vapply(
    landscape_data[landscape_numeric_columns],
    is.numeric,
    logical(1)
  ))) {
    stop("Required landscape metrics must be numeric.", call. = FALSE)
  }

  hex_metadata <- hex_metadata %>%
    mutate(hex_id = as.character(hex_id))

  landscape_scaled <- landscape_data %>%
    mutate(hex_id = as.character(hex_id)) %>%
    drop_na(all_of(required_landscape_columns)) %>%
    filter(pland_urban_2500m <= urban_cover_max) %>%
    mutate(
      across(
        where(is.numeric) & !all_of("hex_id"),
        ~ as.numeric(scale(.x))
      )
    ) %>%
    select(-pland_urban_2500m)

  scaled_focal_columns <- c(
    "pland_agriculture_2500m",
    "area_mn_woodland_2500m",
    "cohesion_woodland_2500m"
  )
  if (nrow(landscape_scaled) < 2L ||
      anyNA(landscape_scaled[scaled_focal_columns]) ||
      any(!is.finite(as.matrix(landscape_scaled[scaled_focal_columns])))) {
    stop(
      "Landscape filtering must leave at least two sites with varying focal metrics.",
      call. = FALSE
    )
  }

  common_ids <- Reduce(
    intersect,
    list(
      rownames(species_matrix),
      hex_metadata$hex_id,
      landscape_scaled$hex_id
    )
  )

  if (length(common_ids) < 2L) {
    stop(
      "Fewer than two sites remain after aligning species, metadata and landscape data.",
      call. = FALSE
    )
  }

  lcbd_species_matrix <- species_matrix[common_ids, , drop = FALSE]

  if (noisy_miner_species %in% colnames(lcbd_species_matrix)) {
    lcbd_species_matrix <- lcbd_species_matrix[
      ,
      colnames(lcbd_species_matrix) != noisy_miner_species,
      drop = FALSE
    ]
  }

  lcbd_species_matrix <- lcbd_species_matrix[
    ,
    colSums(lcbd_species_matrix) > 0,
    drop = FALSE
  ]

  beta_matrices <- calculate_lcbd_beta_matrices(lcbd_species_matrix)
  site_ids <- rownames(lcbd_species_matrix)

  grid <- st_read(hex_path, quiet = TRUE)
  assert_unique_ids(grid, "hex_id", "hexagon grid")

  grid <- grid %>%
    mutate(hex_id = as.character(hex_id))

  grid <- grid %>%
    filter(as.character(hex_id) %in% site_ids) %>%
    st_transform(crs = model_coord_crs)

  # Point-on-surface is intentional for non-rectangular cells. sf warns that
  # attributes are repeated, which is safe here because no attributes change.
  hex_points <- suppressWarnings(st_point_on_surface(grid))
  coordinates <- st_coordinates(hex_points)

  coordinate_lookup <- hex_points %>%
    mutate(
      hex_id = as.character(hex_id),
      x = coordinates[, 1],
      y = coordinates[, 2]
    ) %>%
    st_drop_geometry() %>%
    select(hex_id, x, y)
  assert_unique_ids(coordinate_lookup, "hex_id", "hexagon coordinate lookup")

  model_df <- tibble(
    hex_id = site_ids,
    mean_richness_difference =
      beta_matrices$site_mean_richness_difference[site_ids],
    mean_replacement = beta_matrices$site_mean_replacement[site_ids]
  ) %>%
    mutate(
      replacement_ratio = mean_replacement /
        (mean_replacement + mean_richness_difference)
    ) %>%
    left_join(landscape_scaled, by = "hex_id") %>%
    left_join(hex_metadata, by = "hex_id") %>%
    left_join(coordinate_lookup, by = "hex_id") %>%
    mutate(log_effort = log(pmax(n_checklists, 1))) %>%
    drop_na(x, y) %>%
    arrange(match(hex_id, site_ids))

  final_site_ids <- as.character(model_df$hex_id)
  if (length(final_site_ids) < 2L) {
    stop("Fewer than two sites have usable model coordinates.", call. = FALSE)
  }
  lcbd_species_matrix <- lcbd_species_matrix[final_site_ids, , drop = FALSE]

  beta_matrices <- calculate_lcbd_beta_matrices(lcbd_species_matrix)

  model_df$mean_richness_difference <-
    beta_matrices$site_mean_richness_difference[final_site_ids]
  model_df$mean_replacement <-
    beta_matrices$site_mean_replacement[final_site_ids]
  model_df$replacement_ratio <- model_df$mean_replacement /
    (model_df$mean_replacement + model_df$mean_richness_difference)

  if (anyNA(model_df[c(
    "mean_richness_difference", "mean_replacement", "replacement_ratio"
  )]) || any(!is.finite(as.matrix(model_df[c(
    "mean_richness_difference", "mean_replacement", "replacement_ratio"
  )])))) {
    stop("LCBD site summaries contain undefined values.", call. = FALSE)
  }

  # Preserve identifiers exactly, including non-numeric IDs and leading zeroes.
  model_df$hex_id <- as.character(model_df$hex_id)

  attrition <- tibble(
    step = c(
      "species_matrix_sites",
      "metadata_sites",
      "landscape_sites_after_urban_and_complete_filters",
      "aligned_candidate_sites",
      "final_model_sites",
      "lcbd_response_species"
    ),
    value = c(
      nrow(species_matrix),
      n_distinct(hex_metadata$hex_id),
      nrow(landscape_scaled),
      length(common_ids),
      nrow(model_df),
      ncol(lcbd_species_matrix)
    )
  )

  list(
    candidate_model_df = model_df,
    model_df = model_df,
    lcbd_species_matrix = lcbd_species_matrix,
    beta_matrices = beta_matrices,
    attrition = attrition
  )
}
