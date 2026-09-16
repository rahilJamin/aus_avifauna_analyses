import::from(dplyr, bind_rows, filter)
import::from(magrittr, `%>%`)

#' Run regional convergence analysis across all gradients
#'
#' Applies `analyse_regional_gradient()` to every configured gradient and
#' combines the resulting support, endpoint, contrast, bootstrap, and dimension
#' tables into one checkpoint object.
#'
#' @param site_data Aligned regional site table.
#' @param species_matrix Binary site-by-species matrix aligned to `site_data`.
#' @param gradient_definitions Gradient endpoint definitions.
#' @param component_matrices List of Jaccard component matrices.
#' @param component_sum_error Maximum component-additivity error.
#' @param minimum_sites_per_endpoint Minimum endpoint sites required per
#'   subregion.
#' @param minimum_retained_subregions Minimum retained subregions required per
#'   gradient.
#' @param bootstrap_draws Number of site-level bootstrap draws.
#' @param bootstrap_seed Base bootstrap seed.
#' @param noisy_miner_species Species removed from the response matrix.
#' @param species_matrix_source Path to the species matrix checkpoint.
#'
#' @return A list containing settings and all regional convergence result tables.

analyse_regional_convergence <- function(site_data,
                                         species_matrix,
                                         gradient_definitions,
                                         component_matrices,
                                         component_sum_error,
                                         minimum_sites_per_endpoint,
                                         minimum_retained_subregions,
                                         bootstrap_draws,
                                         bootstrap_seed,
                                         noisy_miner_species,
                                         species_matrix_source) {
  regional_results <- vector("list", nrow(gradient_definitions))

  for (gradient_index in seq_len(nrow(gradient_definitions))) {
    regional_results[[gradient_index]] <- analyse_regional_gradient(
      gradient_definition = gradient_definitions[gradient_index, ],
      site_data = site_data,
      species_matrix = species_matrix,
      component_matrices = component_matrices,
      gradient_index = gradient_index,
      minimum_sites_per_endpoint = minimum_sites_per_endpoint,
      minimum_retained_subregions = minimum_retained_subregions,
      bootstrap_draws = bootstrap_draws,
      bootstrap_seed = bootstrap_seed
    )
  }

  # WORKFLOW CHANGE: PCoA removed at the author's request; convergence contrasts
  # are unchanged.
  list(
    settings = list(
      gradient_definitions = gradient_definitions,
      minimum_sites_per_endpoint = minimum_sites_per_endpoint,
      minimum_retained_subregions = minimum_retained_subregions,
      bootstrap_draws = bootstrap_draws,
      bootstrap_seed = bootstrap_seed,
      regional_unit = "IBRA7 subregion (SUB_NAME_7)",
      species_matrix_source = species_matrix_source,
      species_removed = noisy_miner_species,
      component_sum_error = component_sum_error
    ),
    dimensions = bind_rows(lapply(regional_results, `[[`, "dimensions")),
    regional_support = bind_rows(lapply(regional_results, `[[`, "regional_support")),
    subregion_pair_means = bind_rows(lapply(regional_results, `[[`, "subregion_pair_means")),
    subregion_plot_points = bind_rows(lapply(regional_results, `[[`, "subregion_plot_points")),
    subregion_contrasts = bind_rows(lapply(regional_results, `[[`, "subregion_contrasts")),
    endpoint_summary = bind_rows(lapply(regional_results, `[[`, "endpoint_summary")),
    endpoint_contrasts = bind_rows(lapply(regional_results, `[[`, "endpoint_contrasts")),
    bootstrap_results = bind_rows(lapply(regional_results, `[[`, "bootstrap_results"))
  )
}
