# =============================================================================
# _targets.R
# Reproducible pipeline for the 2020-onward GEB main analysis.
#
# The pipeline derives analysis inputs from the frozen external data:
#   1. Create the 2.5 km hexagon grid from the study extent in Australian Albers.
#   2. Clean the 2020-onward occurrence records and assign them to that grid.
#   3. Estimate Chao2 completeness and build the species matrix.
#   4. Reclassify raw rasters and extract metrics at the Chao2-retained sites.
#   5. Prepare aligned landscape and LCBD model data.
#   6. Fit the Gaussian pool-size-standardised LCBD ensemble.
#   7. Analyse regional convergence among IBRA7 subregions.
#   8. Fit the environment and fourth-corner GLLVM models.
#
# The LCBD GAMMs and the two GLLVM fits are deliberately separate targets. A
# completed heavy target is therefore not rerun when only a downstream summary
# changes. Publication-output targets read the completed analysis checkpoints.
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
  format = "rds",
  # Release large raster/model objects between targets. This is configured here
  # because the equivalent tar_make() argument is deprecated.
  garbage_collection = TRUE
)
function_files <- list.files(
  file.path("R", "functions"),
  pattern = "[.]R$",
  full.names = TRUE
)
tar_source(function_files)


# Stage definitions ----

source(file.path("pipeline", "pipeline_data_processing.R"), local = TRUE)
source(file.path("pipeline", "pipeline_landscape.R"), local = TRUE)
source(file.path("pipeline", "pipeline_lcbd.R"), local = TRUE)
source(file.path("pipeline", "pipeline_regional_convergence.R"), local = TRUE)
source(file.path("pipeline", "pipeline_gllvm.R"), local = TRUE)
source(file.path("pipeline", "pipeline_outputs.R"), local = TRUE)


# Complete pipeline ----

c(
  data_processing_targets,
  landscape_targets,
  model_data_targets,
  lcbd_targets,
  regional_convergence_targets,
  gllvm_targets,
  output_targets
)
