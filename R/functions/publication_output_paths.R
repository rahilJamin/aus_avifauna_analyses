#' List publication outputs by analysis stage
#'
#' Defines the files produced from completed model checkpoints. The same paths
#' are returned by the output functions and tracked by the `targets` file
#' targets, so a missing export is detected rather than silently skipped.
#'
#' @param cfg Main workflow configuration with table, figure, and diagnostic
#'   directories.
#' @return A named list of file paths for LCBD, regional convergence, and GLLVM.
publication_output_paths <- function(cfg) {
  regional_gradients <- c(
    "agriculture", "patch_area", "cohesion", "noisy_miner"
  )

  list(
    lcbd = c(
      file.path(cfg$table_dir, "04_local_uniqueness_results.csv"),
      file.path(cfg$table_dir, "04_reference_pool_allocations.csv"),
      file.path(cfg$table_dir, "04_lcbd_curve_summary.csv"),
      file.path(cfg$figure_dir, "lcbd_plots", "04_local_uniqueness.svg")
    ),
    regional = c(
      file.path(cfg$table_dir, "05_regional_support.csv"),
      file.path(cfg$table_dir, "05_regional_results.csv"),
      file.path(cfg$table_dir, "05_subregion_pair_means.csv"),
      file.path(cfg$table_dir, "05_subregion_results.csv"),
      file.path(cfg$table_dir, "05_figure_statistics.csv"),
      file.path(cfg$figure_dir, "regional_plots", "05_regional_convergence_2x2.svg"),
      file.path(
        cfg$figure_dir,
        "regional_plots",
        paste0("05_regional_convergence_", regional_gradients, ".svg")
      )
    ),
    gllvm = c(
      file.path(cfg$table_dir, "06_gllvm_species_missing_traits.txt"),
      file.path(cfg$table_dir, "06_gllvm_model_dimensions_with_dbmem.csv"),
      file.path(cfg$table_dir, "06_gllvm_species_coefficients.csv"),
      file.path(cfg$table_dir, "06_fourth_corner_coefficients.csv"),
      file.path(cfg$diagnostic_dir, "06_gllvm_basic_diagnostics.pdf"),
      file.path(cfg$figure_dir, "gllvm_plots", "beta_heatmap_selected_grouped.svg"),
      file.path(cfg$figure_dir, "gllvm_plots", "heatmap_standalone_legend.svg"),
      file.path(cfg$figure_dir, "gllvm_plots", "beta_heatmap_all.svg"),
      file.path(cfg$figure_dir, "gllvm_plots", "trait_heatmap_selected.svg"),
      file.path(cfg$figure_dir, "gllvm_plots", "gllvm_habitat_direction.svg")
    )
  )
}
