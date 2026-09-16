#' Validate ecological-group assignments for modelled species
#'
#' Ensures that every modelled species has exactly one ecological-group
#' assignment before group summaries or figures are produced.
#'
#' @param model_species Character vector of species in the fitted model output.
#' @param species_group_key Data frame containing `species` and
#'   `ecological_group` columns.
#'
#' @return `species_group_key`, invisibly, when assignments are complete.
validate_gllvm_species_groups <- function(model_species, species_group_key) {
  required_columns <- c("species", "ecological_group")
  missing_columns <- setdiff(required_columns, names(species_group_key))
  if (length(missing_columns) > 0L) {
    stop(
      "The ecological-group key lacks: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  model_species <- unique(as.character(model_species))
  key_species <- as.character(species_group_key$species)

  if (length(model_species) < 1L || anyNA(model_species) ||
      any(!nzchar(model_species))) {
    stop("Modelled species names must be present and non-empty.", call. = FALSE)
  }

  if (anyNA(key_species) || any(!nzchar(key_species)) ||
      anyNA(species_group_key$ecological_group) ||
      any(!nzchar(as.character(species_group_key$ecological_group)))) {
    stop("The ecological-group key contains missing species or groups.", call. = FALSE)
  }

  duplicated_assignments <- unique(key_species[duplicated(key_species)])
  if (length(duplicated_assignments) > 0L) {
    stop(
      "Species assigned to more than one ecological group: ",
      paste(duplicated_assignments, collapse = ", "),
      call. = FALSE
    )
  }

  missing_species <- setdiff(model_species, key_species)
  if (length(missing_species) > 0L) {
    stop(
      "Modelled species missing from the ecological-group key: ",
      paste(missing_species, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(species_group_key)
}
