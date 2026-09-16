#' Select spatial dbMEM variables for GLLVM fitting
#'
#' Builds positive spatial eigenvectors from model-site coordinates, runs the
#' global permutation test against the Hellinger-transformed community matrix,
#' and performs forward selection only when the global test is supported.
#'
#' @param species_matrix Binary site-by-species matrix.
#' @param coordinates Numeric matrix or data frame with two coordinate columns
#'   in the same row order as `species_matrix`.
#' @param seed Random seed used by permutation tests.
#' @param permutations Number of permutations.
#' @param global_alpha Significance level for the global dbMEM test.
#' @param selection_alpha Significance level used during forward selection.
#'
#' @return A list containing the global test, adjusted R-squared, forward
#'   selection result, selected names, and selected dbMEM columns.
select_gllvm_dbmem <- function(species_matrix,
                               coordinates,
                               seed,
                               permutations,
                               global_alpha,
                               selection_alpha) {
  species_matrix <- as.matrix(species_matrix)
  coordinates <- as.matrix(coordinates)

  if (nrow(species_matrix) != nrow(coordinates) || ncol(coordinates) != 2L) {
    stop(
      "GLLVM coordinates must contain two columns aligned to the species matrix.",
      call. = FALSE
    )
  }

  if (!is.numeric(coordinates) || anyNA(coordinates) ||
      any(!is.finite(coordinates)) || anyDuplicated(as.data.frame(coordinates))) {
    stop("GLLVM coordinates must be finite and unique.", call. = FALSE)
  }

  valid_probability <- function(value) {
    is.numeric(value) && length(value) == 1L && !is.na(value) &&
      is.finite(value) && value > 0 && value < 1
  }

  if (!valid_probability(global_alpha) || !valid_probability(selection_alpha)) {
    stop("dbMEM alpha values must be between zero and one.", call. = FALSE)
  }

  if (!is.numeric(permutations) || length(permutations) != 1L ||
      is.na(permutations) || permutations < 1L ||
      permutations != floor(permutations)) {
    stop("`permutations` must be a positive whole number.", call. = FALSE)
  }

  set.seed(seed)

  mem_full <- adespatial::dbmem(
    coordinates,
    MEM.autocor = "positive",
    silent = FALSE
  )
  community_hellinger <- vegan::decostand(species_matrix, "hellinger")

  global_model <- vegan::rda(
    community_hellinger ~ .,
    data = as.data.frame(mem_full)
  )
  global_test <- stats::anova(global_model, permutations = permutations)
  adjusted_r2 <- vegan::RsquareAdj(global_model)$adj.r.squared

  selected_mem_names <- character()
  forward_selection <- NULL
  global_p_value <- global_test$`Pr(>F)`[1]

  if (!is.na(global_p_value) && global_p_value < global_alpha) {
    forward_selection <- adespatial::forward.sel(
      Y = community_hellinger,
      X = as.data.frame(mem_full),
      adjR2thresh = adjusted_r2,
      alpha = selection_alpha,
      nperm = permutations
    )

    selected_mem_names <- as.character(forward_selection$variables)
  }

  selected_mems <- as.data.frame(mem_full)[
    ,
    selected_mem_names,
    drop = FALSE
  ]

  list(
    global_test = global_test,
    adjusted_r2 = adjusted_r2,
    forward_selection = forward_selection,
    selected_mem_names = selected_mem_names,
    selected_mems = selected_mems
  )
}
