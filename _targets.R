# =============================================================================
# _targets.R
# Reproducible pipeline for the 2020-onward GEB main analysis.
#
# The pipeline follows the validated manual workflow:
#   1. Clean occurrence records.
#   2. Assign records to the 2.5 km hexagon grid.
#   3. Estimate Chao2 completeness and build the species matrix.
#   4. Prepare aligned landscape and LCBD model data.
#   5. Fit the Gaussian pool-size-standardised LCBD ensemble.
#   6. Analyse regional convergence among IBRA7 subregions.
#   7. Fit the environment and fourth-corner GLLVM models.
#
# The LCBD GAMMs and the two GLLVM fits are deliberately separate targets. A
# completed heavy target is therefore not rerun when only a downstream summary
# changes. The final Quarto report reads the completed analysis checkpoints.
# =============================================================================

library(targets)


# Configuration and functions ----

source(file.path("R", "00_main_config.R"))
tar_option_set(
  packages = c(
    "adespatial",
    "CoordinateCleaner",
    "dplyr",
    "gllvm",
    "iNEXT",
    "import",
    "magrittr",
    "mgcv",
    "readr",
    "sf",
    "stringr",
    "tibble",
    "tidyr",
    "tidyselect",
    "TMB",
    "vegan"
  ),
  format = "rds"
)
function_files <- list.files(
  file.path("R", "functions"),
  pattern = "[.]R$",
  full.names = TRUE
)
tar_source(function_files)


# Stage definitions ----

source(file.path("pipeline", "pipeline_data_processing.R"), local = TRUE)
source(file.path("pipeline", "pipeline_lcbd.R"), local = TRUE)
source(file.path("pipeline", "pipeline_regional_convergence.R"), local = TRUE)
source(file.path("pipeline", "pipeline_gllvm.R"), local = TRUE)
source(file.path("pipeline", "pipeline_outputs.R"), local = TRUE)


# Complete pipeline ----

c(
  data_processing_targets,
  lcbd_targets,
  regional_convergence_targets,
  gllvm_targets,
  output_targets
)
