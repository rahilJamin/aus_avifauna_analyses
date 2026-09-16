import::from(magrittr, `%>%`)
import::from(dplyr, count, distinct, filter, left_join, mutate)
import::from(tidyr, replace_na)

#' Calculate Noisy Miner detection by hexagon
#'
#' Uses the same incidence-unit definition as the Chao2/iNEXT stage. The output
#' retains the historical checklist column names, where a checklist means one
#' unique site (`hex_id`) by incidence-unit combination.
#'
#' @param hex_occ Hex-level occurrence table. Each `hex_id` is a 2.5 km grid cell
#'   treated as one ecological site.
#' @param noisy_miner_species Scientific name identifying Noisy Miner records.
#' @param incidence_unit Incidence-unit definition passed to
#'   `add_incidence_unit()`.
#'
#' @return A data frame with `hex_id`, `n_checklists`, `n_miner_checklists`, and
#'   `miner_detection_rate`.

calculate_miner_detection <- function(hex_occ,
                                      noisy_miner_species,
                                      incidence_unit) {
  if (!is.character(noisy_miner_species) ||
      length(noisy_miner_species) != 1L || is.na(noisy_miner_species) ||
      !nzchar(noisy_miner_species)) {
    stop("`noisy_miner_species` must be one non-empty species name.", call. = FALSE)
  }

  required_columns <- c("hex_id", "species", "date")
  missing_columns <- setdiff(required_columns, names(hex_occ))
  if (length(missing_columns) > 0L) {
    stop(
      "The occurrence table used for Noisy Miner detection lacks: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(hex_occ) == 0L || anyNA(hex_occ$hex_id) ||
      any(!nzchar(as.character(hex_occ$hex_id))) || anyNA(hex_occ$species)) {
    stop("Noisy Miner detection requires non-empty site and species records.", call. = FALSE)
  }

  hex_occ$hex_id <- as.character(hex_occ$hex_id)
  dat <- add_incidence_unit(hex_occ, incidence_unit)

  # Use the same sampling-unit definition as the Chao2/iNEXT step. The output
  # keeps the historical "checklist" names for downstream compatibility.
  checklists <- dat %>%
    distinct(hex_id, incidence_unit)

  effort <- checklists %>%
    count(hex_id, name = "n_checklists")

  miner <- dat %>%
    filter(species == noisy_miner_species) %>%
    distinct(hex_id, incidence_unit) %>%
    count(hex_id, name = "n_miner_checklists")

  effort %>%
    left_join(miner, by = "hex_id") %>%
    mutate(
      n_miner_checklists = replace_na(n_miner_checklists, 0L),
      miner_detection_rate = n_miner_checklists / n_checklists
    )
}
