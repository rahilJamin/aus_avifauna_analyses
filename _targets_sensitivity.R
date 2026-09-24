# =============================================================================
# _targets_sensitivity.R
# Separate, supplementary robustness analyses for the GEB workflow.
#
# Run only after the main pipeline has completed. This entry point reads the
# main checkpoints and never alters _targets, _targets.R or the main outputs.
# =============================================================================

library(targets)

source(file.path("R", "00_main_config.R"))
source(file.path("R", "01_sensitivity_config.R"))

tar_option_set(
  packages = c(
    "adespatial", "callr", "dplyr", "ggplot2", "gllvm", "iNEXT", "import",
    "magrittr", "mgcv", "readr", "sf", "spdep", "TMB", "tibble",
    "tidyr", "tidyselect", "vegan"
  ),
  format = "rds"
)

main_function_files <- list.files(
  file.path("R", "functions"), pattern = "[.]R$", full.names = TRUE
)
sensitivity_source_files <- c(
  main_function_files,
  list.files(file.path("R", "sensitivity"), pattern = "[.]R$", full.names = TRUE)
)
tar_source(sensitivity_source_files)

source(file.path("pipeline", "pipeline_sensitivity.R"), local = TRUE)

sensitivity_targets
