# =============================================================================
# pipeline_outputs.R
# Targets for publication tables, figures, and diagnostics.
# Each stage reads its completed model checkpoint and returns tracked files.
# =============================================================================

output_targets <- list(
  tar_target(
    lcbd_publication_files,
    {
      invisible(lcbd_results_file)
      produce_lcbd_outputs(cfg)
    },
    format = "file"
  ),
  tar_target(
    regional_publication_files,
    {
      invisible(regional_results_file)
      produce_regional_outputs(cfg)
    },
    format = "file"
  ),
  tar_target(
    gllvm_publication_files,
    {
      invisible(gllvm_metadata_file)
      invisible(gllvm_environment_model_file)
      invisible(gllvm_fourth_corner_model_file)
      produce_gllvm_outputs(cfg)
    },
    format = "file"
  ),
  tar_target(
    publication_files,
    c(
      lcbd_publication_files,
      regional_publication_files,
      gllvm_publication_files
    )
  )
)
