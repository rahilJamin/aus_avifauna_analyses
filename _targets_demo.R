# =============================================================================
# _targets_demo.R
# Isolated, reduced-data execution of the complete main workflow.
#
# This entry point uses demo/data/inputs, writes only to demo/outputs, and has
# its own _targets_demo/ cache. It never reads from or writes to the main cache
# or output directories. The analysis functions and formulas are shared with
# the main pipeline; only the input subset and LCBD ensemble size are reduced.
# =============================================================================

library(targets)

source(file.path("R", "00_main_config.R"))
source(file.path("R", "02_demo_config.R"))

tar_option_set(
  packages = c(
    "adespatial", "CoordinateCleaner", "cowplot", "dplyr", "forcats",
    "ggh4x", "ggplot2", "gllvm", "iNEXT", "import", "landscapemetrics",
    "magrittr", "mgcv", "patchwork", "rcartocolor", "readr", "scales",
    "sf", "stringr", "terra", "tibble", "tidyr", "tidyselect", "TMB",
    "vegan"
  ),
  format = "rds"
)

function_files <- list.files(
  file.path("R", "functions"), pattern = "[.]R$", full.names = TRUE
)
tar_source(function_files)

source(file.path("pipeline", "pipeline_data_processing.R"), local = TRUE)
source(file.path("pipeline", "pipeline_landscape.R"), local = TRUE)
source(file.path("pipeline", "pipeline_lcbd.R"), local = TRUE)
source(file.path("pipeline", "pipeline_regional_convergence.R"), local = TRUE)
source(file.path("pipeline", "pipeline_gllvm.R"), local = TRUE)
source(file.path("pipeline", "pipeline_outputs.R"), local = TRUE)

c(
  data_processing_targets,
  landscape_targets,
  model_data_targets,
  lcbd_targets,
  regional_convergence_targets,
  gllvm_targets,
  output_targets
)
