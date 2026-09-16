import::from(dplyr, arrange, filter, mutate, select, transmute)
import::from(magrittr, `%>%`)
import::from(tibble, tibble)
import::from(tidyselect, all_of)

#' Prepare aligned matrices and traits for GLLVM fitting
#'
#' Aligns the prepared site table and species matrix, removes sites with
#' incomplete environmental predictors, applies the configured occupancy
#' threshold to that final site set, and aligns valid species traits. This
#' function performs data preparation only and does not fit a model.
#'
#' @param model_df Site-level model data containing `hex_id`, `ibra_sub`, and
#'   the environmental predictors.
#' @param species_matrix Binary site-by-species matrix with `hex_id` row names.
#' @param traits_raw Species trait table containing `species1`, `mass`, and
#'   `hand_wing_index`.
#' @param base_predictors Character vector of environmental predictor columns.
#' @param occupancy_threshold Minimum proportion of final model sites at which
#'   a species must occur. Equality is retained.
#' @param noisy_miner_species Scientific name used to report whether Noisy Miner
#'   remains in the GLLVM response.
#'
#' @return A list containing aligned `X`, `Y`, `gllvm_df`, scaled traits,
#'   `trait_matrix`, species counts, dimensions, and trait-matching diagnostics.
prepare_gllvm_model_data <- function(model_df,
                                     species_matrix,
                                     traits_raw,
                                     base_predictors,
                                     occupancy_threshold,
                                     noisy_miner_species) {
  if (!is.character(base_predictors) || length(base_predictors) < 1L ||
      anyNA(base_predictors) || any(!nzchar(base_predictors)) ||
      anyDuplicated(base_predictors)) {
    stop("`base_predictors` must contain unique, non-empty column names.", call. = FALSE)
  }
  if (!is.character(noisy_miner_species) ||
      length(noisy_miner_species) != 1L || is.na(noisy_miner_species) ||
      !nzchar(noisy_miner_species)) {
    stop("`noisy_miner_species` must be one non-empty species name.", call. = FALSE)
  }

  required_model_columns <- c(
    "hex_id",
    "ibra_sub",
    "x",
    "y",
    base_predictors
  )
  missing_model_columns <- setdiff(required_model_columns, names(model_df))
  if (length(missing_model_columns) > 0L) {
    stop(
      "The GLLVM site table lacks: ",
      paste(missing_model_columns, collapse = ", "),
      call. = FALSE
    )
  }

  assert_unique_ids(model_df, "hex_id", "GLLVM site table")

  if (!is.numeric(occupancy_threshold) || length(occupancy_threshold) != 1L ||
      is.na(occupancy_threshold) || occupancy_threshold <= 0 ||
      occupancy_threshold > 1) {
    stop("`occupancy_threshold` must be greater than zero and at most one.", call. = FALSE)
  }

  species_matrix <- as.matrix(species_matrix)
  if (is.null(rownames(species_matrix)) || is.null(colnames(species_matrix))) {
    stop("The GLLVM species matrix must have row and column names.", call. = FALSE)
  }
  if (anyDuplicated(rownames(species_matrix))) {
    stop("The GLLVM species matrix contains duplicated site IDs.", call. = FALSE)
  }
  if (anyDuplicated(colnames(species_matrix))) {
    stop("The GLLVM species matrix contains duplicated species names.", call. = FALSE)
  }
  if (!is.numeric(species_matrix) || anyNA(species_matrix) ||
      !all(is.finite(species_matrix)) ||
      !all(species_matrix %in% c(0, 1))) {
    stop("The GLLVM species matrix must contain only zeroes and ones.", call. = FALSE)
  }

  model_df <- model_df %>%
    mutate(hex_id = as.character(hex_id)) %>%
    filter(!is.na(ibra_sub)) %>%
    mutate(ibra_sub = droplevels(factor(ibra_sub)))

  common_ids <- intersect(model_df$hex_id, rownames(species_matrix))
  if (length(common_ids) < 2L) {
    stop("Fewer than two sites remain after GLLVM site alignment.", call. = FALSE)
  }

  gllvm_df <- model_df %>%
    filter(hex_id %in% common_ids) %>%
    arrange(match(hex_id, common_ids))
  Y <- species_matrix[common_ids, , drop = FALSE]

  numeric_model_columns <- c("x", "y", base_predictors)
  numeric_columns <- vapply(
    gllvm_df[numeric_model_columns],
    is.numeric,
    logical(1)
  )
  if (!all(numeric_columns)) {
    stop(
      "GLLVM coordinates and environmental predictors must be numeric: ",
      paste(names(numeric_columns)[!numeric_columns], collapse = ", "),
      call. = FALSE
    )
  }
  if (anyNA(gllvm_df[c("x", "y")]) ||
      any(!is.finite(as.matrix(gllvm_df[c("x", "y")]))) ) {
    stop("GLLVM model coordinates must be present and finite.", call. = FALSE)
  }
  if (anyDuplicated(gllvm_df[c("x", "y")])) {
    stop("GLLVM model sites must have unique coordinate pairs.", call. = FALSE)
  }

  X <- gllvm_df %>%
    select(all_of(base_predictors)) %>%
    as.data.frame()

  complete_rows <- stats::complete.cases(X) &
    apply(is.finite(as.matrix(X)), 1L, all)
  X <- X[complete_rows, , drop = FALSE]
  Y <- Y[complete_rows, , drop = FALSE]
  gllvm_df <- gllvm_df[complete_rows, , drop = FALSE]
  gllvm_df$ibra_sub <- droplevels(gllvm_df$ibra_sub)

  if (nrow(Y) < 2L) {
    stop("Fewer than two complete sites remain for GLLVM fitting.", call. = FALSE)
  }
  if (nlevels(gllvm_df$ibra_sub) < 2L) {
    stop("GLLVM fitting requires at least two retained IBRA subregions.", call. = FALSE)
  }
  if (!identical(rownames(Y), gllvm_df$hex_id)) {
    stop("GLLVM site rows are not aligned with the species matrix.", call. = FALSE)
  }

  site_frequency <- colMeans(Y)
  species_keep <- names(site_frequency[site_frequency >= occupancy_threshold])
  if (length(species_keep) == 0L) {
    stop("No species were retained by the GLLVM occupancy filter.", call. = FALSE)
  }
  Y <- Y[, species_keep, drop = FALSE]
  minimum_required_presences <- ceiling(occupancy_threshold * nrow(Y))

  required_trait_columns <- c("species1", "mass", "hand_wing_index")
  missing_trait_columns <- setdiff(required_trait_columns, names(traits_raw))
  if (length(missing_trait_columns) > 0L) {
    stop(
      "The trait table lacks: ",
      paste(missing_trait_columns, collapse = ", "),
      call. = FALSE
    )
  }
  assert_unique_ids(traits_raw, "species1", "trait table")

  missing_traits <- setdiff(colnames(Y), as.character(traits_raw$species1))
  matched_species <- intersect(colnames(Y), as.character(traits_raw$species1))
  if (length(matched_species) < 2L) {
    stop("Fewer than two modelled species have matching trait records.", call. = FALSE)
  }

  matched_traits <- traits_raw %>%
    filter(species1 %in% matched_species) %>%
    arrange(match(species1, matched_species))

  if (!is.numeric(matched_traits$mass) ||
      !is.numeric(matched_traits$hand_wing_index)) {
    stop("Mass and hand-wing index must be numeric.", call. = FALSE)
  }

  invalid_traits <- !is.finite(matched_traits$mass) |
    matched_traits$mass <= 0 |
    !is.finite(matched_traits$hand_wing_index)
  if (any(invalid_traits)) {
    stop(
      "Invalid mass or hand-wing index for trait-matched species: ",
      paste(matched_traits$species1[invalid_traits], collapse = ", "),
      ". Mass must be positive and both traits must be finite.",
      call. = FALSE
    )
  }

  mass_sd <- stats::sd(log(matched_traits$mass))
  hwi_sd <- stats::sd(matched_traits$hand_wing_index)
  if (!is.finite(mass_sd) || !is.finite(hwi_sd) ||
      mass_sd == 0 || hwi_sd == 0) {
    stop("Each GLLVM trait must vary among retained species.", call. = FALSE)
  }

  Y <- Y[, matched_species, drop = FALSE]
  traits <- matched_traits %>%
    transmute(
      species = as.character(species1),
      log_mass = as.numeric(scale(log(mass))),
      HWI = as.numeric(scale(hand_wing_index))
    )

  if (!identical(traits$species, colnames(Y)) ||
      anyNA(traits$log_mass) || anyNA(traits$HWI)) {
    stop("GLLVM traits are not aligned with the response matrix.", call. = FALSE)
  }

  trait_matrix <- traits %>%
    select(-species) %>%
    as.data.frame()
  rownames(trait_matrix) <- traits$species

  species_counts <- tibble(
    species = colnames(Y),
    n_presence = as.integer(colSums(Y)),
    n_absence = as.integer(nrow(Y) - colSums(Y)),
    occupancy = n_presence / nrow(Y)
  )

  if (any(species_counts$n_presence < minimum_required_presences)) {
    stop("A retained species falls below the configured occupancy threshold.", call. = FALSE)
  }

  dimensions <- tibble(
    n_sites = nrow(Y),
    n_species_after_occupancy_and_traits = ncol(Y),
    n_environment_predictors = ncol(X),
    n_traits = ncol(trait_matrix),
    n_ibra_subregions = nlevels(gllvm_df$ibra_sub),
    occupancy_threshold = occupancy_threshold,
    minimum_required_presences = minimum_required_presences,
    noisy_miner_retained = noisy_miner_species %in% colnames(Y)
  )

  list(
    X = X,
    Y = Y,
    gllvm_df = gllvm_df,
    traits = traits,
    trait_matrix = trait_matrix,
    species_counts = species_counts,
    dimensions = dimensions,
    species_order = colnames(Y),
    missing_traits = missing_traits,
    minimum_required_presences = minimum_required_presences
  )
}
