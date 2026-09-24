# =============================================================================
# Script: run_sensitivity.R
# Purpose: User-run entry points for the supplementary sensitivity pipeline.
# Inputs: Completed main targets outputs and inputs declared in 00_main_config.R.
# Outputs: outputs/sensitivity/ and _targets_sensitivity/.
# Execution note: Run manually from the repository root after the main pipeline.
# =============================================================================

library(targets)

sensitivity_script <- "_targets_sensitivity.R"
sensitivity_store <- "_targets_sensitivity"

# Inspect planned targets before committing to a long sensitivity run.
targets::tar_manifest(script = sensitivity_script)

# Run the complete sensitivity workflow, including all heavy LCBD scenarios,
# the matched no-dbMEM environmental GLLVM, and appendix-ready SVG exports.
# Remove the leading # when ready.
# targets::tar_make(script = sensitivity_script, store = sensitivity_store)

# Useful staged runs. Each command reuses completed targets in the separate
# store; they do not re-run the main pipeline.
# targets::tar_make(names = sensitivity_lcbd_result,
#   script = sensitivity_script, store = sensitivity_store)
# targets::tar_make(names = miner_lcbd_comparison,
#   script = sensitivity_script, store = sensitivity_store)
# targets::tar_make(names = lcbd_spatial_control,
#   script = sensitivity_script, store = sensitivity_store)
# targets::tar_make(names = gllvm_without_dbmem_file,
#   script = sensitivity_script, store = sensitivity_store)
# targets::tar_make(names = sensitivity_publication_files,
#   script = sensitivity_script, store = sensitivity_store)
