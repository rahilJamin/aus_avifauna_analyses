import::from(magrittr, `%>%`)
import::from(dplyr, distinct, filter, mutate)
import::from(tidyr, pivot_wider)

#' Build a binary species-by-hexagon matrix
#'
#' Converts a hex-level occurrence table into a presence-absence matrix for the
#' retained hexagons. Matrix rows are hexagons and columns are species.
#'
#' @param hex_occ Hex-level occurrence table with `hex_id` and `species` columns.
#'   `hex_id` is the unique identifier for a 2.5 km grid cell used as a site.
#' @param valid_hex_ids Hexagon IDs retained after sampling and completeness
#'   filters.
#'
#' @return A binary matrix with row names equal to retained `hex_id` values and
#'   one column per retained species.

make_species_matrix <- function(hex_occ, valid_hex_ids) {
  required_columns <- c("hex_id", "species")
  missing_columns <- setdiff(required_columns, names(hex_occ))
  if (length(missing_columns) > 0L) {
    stop(
      "The occurrence table lacks: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  valid_hex_ids <- as.character(valid_hex_ids)
  if (length(valid_hex_ids) < 1L || anyNA(valid_hex_ids) ||
      any(!nzchar(valid_hex_ids)) || anyDuplicated(valid_hex_ids)) {
    stop("Retained hexagon IDs must be unique, non-empty values.", call. = FALSE)
  }
  if (nrow(hex_occ) == 0L || anyNA(hex_occ$hex_id) || anyNA(hex_occ$species) ||
      any(!nzchar(as.character(hex_occ$hex_id))) ||
      any(!nzchar(as.character(hex_occ$species)))) {
    stop("Species-matrix input requires non-empty site and species values.", call. = FALSE)
  }

  hex_occ$hex_id <- as.character(hex_occ$hex_id)
  missing_valid_ids <- setdiff(valid_hex_ids, unique(hex_occ$hex_id))
  if (length(missing_valid_ids) > 0L) {
    stop(
      "Retained hexagon IDs absent from the occurrence table: ",
      paste(missing_valid_ids, collapse = ", "),
      call. = FALSE
    )
  }

  species_presence <- hex_occ %>%
    filter(hex_id %in% valid_hex_ids) %>%
    distinct(hex_id, species) %>%
    mutate(presence = 1L) %>%
    pivot_wider(names_from = species, values_from = presence,
                values_fill = 0)

  mat <- as.matrix(species_presence[, -1, drop = FALSE])
  rownames(mat) <- as.character(species_presence$hex_id)
  mat <- mat[valid_hex_ids, , drop = FALSE]
  mat[, colSums(mat) > 0, drop = FALSE]
}
