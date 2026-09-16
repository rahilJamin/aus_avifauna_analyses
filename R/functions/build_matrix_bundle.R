import::from(magrittr, `%>%`)
import::from(dplyr, filter, left_join, mutate, n_distinct)
import::from(tibble, tibble)
import::from(tidyr, replace_na)

#' Build the iNEXT-filtered species-matrix bundle
#'
#' Calculates Chao2 richness/completeness by hexagon, keeps sufficiently sampled
#' hexagons, builds the retained binary species matrix, calculates Noisy Miner
#' detection with the same incidence-unit definition, and returns metadata and
#' attrition summaries for downstream models.
#'
#' @param hex_occ Hex-level occurrence table from the hex-assignment stage. Each
#'   `hex_id` is a 2.5 km grid cell treated as one ecological site.
#' @param incidence_unit Incidence-unit definition passed to
#'   `add_incidence_unit()`.
#' @param min_sampling_units Minimum number of incidence units required for a
#'   hexagon.
#' @param min_species Minimum observed species richness required for a hexagon.
#' @param min_completeness Minimum Chao2 sampling completeness required for a
#'   hexagon.
#' @param noisy_miner_species Scientific name used to calculate Noisy Miner
#'   detection.
#' @param show_progress Logical; passed to Chao2 calculation.
#'
#' @return A list containing `chao`, `valid`, `species_matrix`, `hex_metadata`,
#'   and `attrition`.

build_matrix_bundle <- function(hex_occ,
                                incidence_unit,
                                min_sampling_units,
                                min_species,
                                min_completeness,
                                noisy_miner_species,
                                show_progress = TRUE) {
  chao <- calculate_chao2_table(hex_occ, incidence_unit, show_progress)

  valid <- retain_sampled_hexagons(
    chao = chao,
    min_sampling_units = min_sampling_units,
    min_species = min_species,
    min_completeness = min_completeness
  )

  if (nrow(valid) == 0L) {
    stop(
      "No hexagons passed the Chao2 sampling and completeness filters.",
      call. = FALSE
    )
  }

  species_matrix <- make_species_matrix(hex_occ, valid$hex_id)
  miner <- calculate_miner_detection(
    hex_occ = hex_occ,
    noisy_miner_species = noisy_miner_species,
    incidence_unit = incidence_unit
  )

  valid_order <- match(rownames(species_matrix), as.character(valid$hex_id))

  if (anyNA(valid_order)) {
    stop(
      "Could not align species matrix rows with Chao2 metadata.",
      call. = FALSE
    )
  }

  # Keep metadata rows in exactly the same order as the species matrix rows.
  hex_metadata <- valid[valid_order, , drop = FALSE] %>%
    left_join(miner, by = "hex_id") %>%
    mutate(
      n_checklists = replace_na(n_checklists, 0L),
      miner_detection_rate = replace_na(miner_detection_rate, 0)
    )

  if (!identical(as.character(hex_metadata$hex_id), rownames(species_matrix))) {
    stop("Hexagon metadata and species-matrix rows are misaligned.", call. = FALSE)
  }

  attrition <- tibble(
    step = c(
      "occupied_hexagons_before_chao2",
      "hexagons_with_chao2_estimates",
      "hexagons_after_species_and_completeness_filters",
      "species_matrix_sites",
      "species_matrix_species"
    ),
    value = c(
      n_distinct(hex_occ$hex_id),
      sum(chao$chao_status == "ok", na.rm = TRUE),
      nrow(valid),
      nrow(species_matrix),
      ncol(species_matrix)
    )
  )

  list(
    chao = chao,
    valid = valid,
    species_matrix = species_matrix,
    hex_metadata = hex_metadata,
    attrition = attrition
  )
}

#' Retain sufficiently sampled hexagons
#'
#' Applies the sampling-unit, observed-richness, and Chao2-completeness rules to
#' a completed Chao2 table. Values exactly equal to any configured cutoff are
#' retained.
#'
#' @param chao Chao2 summary with one row per hexagon.
#' @param min_sampling_units Minimum retained incidence-unit count.
#' @param min_species Minimum retained observed species richness.
#' @param min_completeness Minimum retained Chao2 completeness in `[0, 1]`.
#'
#' @return The rows with successful Chao2 estimates that meet every cutoff.
retain_sampled_hexagons <- function(chao,
                                    min_sampling_units,
                                    min_species,
                                    min_completeness) {
  required_columns <- c(
    "hex_id",
    "chao_status",
    "sampling_units",
    "S_obs",
    "chao2_completeness"
  )
  missing_columns <- setdiff(required_columns, names(chao))
  if (length(missing_columns) > 0L) {
    stop(
      "The Chao2 table lacks: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  valid_count <- function(value) {
    is.numeric(value) && length(value) == 1L && !is.na(value) &&
      is.finite(value) && value >= 1 && value == floor(value)
  }
  if (!valid_count(min_sampling_units) || !valid_count(min_species)) {
    stop("Minimum sampling units and species must be positive whole numbers.", call. = FALSE)
  }
  if (!is.numeric(min_completeness) || length(min_completeness) != 1L ||
      is.na(min_completeness) || !is.finite(min_completeness) ||
      min_completeness < 0 || min_completeness > 1) {
    stop("Minimum Chao2 completeness must be between zero and one.", call. = FALSE)
  }
  if (anyNA(chao$hex_id) || any(!nzchar(as.character(chao$hex_id))) ||
      anyDuplicated(as.character(chao$hex_id))) {
    stop("The Chao2 table must contain unique, non-missing hexagon IDs.", call. = FALSE)
  }

  numeric_summary_columns <- c(
    "sampling_units", "S_obs", "chao2_completeness"
  )
  if (!all(vapply(chao[numeric_summary_columns], is.numeric, logical(1))) ||
      anyNA(chao$chao_status) ||
      any(!chao$chao_status %in% c("ok", "failed"))) {
    stop("Successful Chao2 summaries must use numeric filtering columns.", call. = FALSE)
  }
  successful <- chao$chao_status == "ok"
  if (anyNA(chao[successful, numeric_summary_columns, drop = FALSE]) ||
      any(!is.finite(as.matrix(
        chao[successful, numeric_summary_columns, drop = FALSE]
      ))) ||
      any(chao$chao2_completeness[successful] < 0) ||
      any(chao$chao2_completeness[successful] > 1)) {
    stop("Successful Chao2 summaries contain invalid filtering values.", call. = FALSE)
  }

  chao %>%
    filter(
      chao_status == "ok",
      sampling_units >= min_sampling_units,
      S_obs >= min_species,
      chao2_completeness >= min_completeness
    )
}
