# GEB: Australian bird-community analysis

This repository contains the code used to analyse Australian bird communities
recorded from **1 January 2020 to 28 February 2026**. Using 2.5 km hexagons as
sampling sites, the workflow examines how local community uniqueness and
regional community differences vary with agricultural cover, woodland patch
area, woodland cohesion, and Noisy Miner detection. It also models individual
species responses to habitat and how those responses relate to species traits.

The repository produces the analysis tables and **publication figures** for
three linked analyses:

| Analysis | Main output |
| --- | --- |
| Local contributions to beta diversity (LCBD) | Gaussian model summaries and curves for total, replacement, and richness-difference LCBD across repeated spatial reference pools. |
| Regional convergence | Low-versus-high gradient comparisons of bird-community dissimilarity among IBRA7 subregions. |
| Joint species-distribution models (GLLVM) | Species coefficients, fourth-corner trait effects, habitat-response figures, and model diagnostics. |

The analyses share a common occurrence-cleaning, sampling-completeness, and
landscape-data preparation workflow. Each model then applies the filters needed
for its own analysis. Sensitivity analyses are outside this repository.

## Repository structure

```text
geb_repo/
├── _targets.R                 # Pipeline entry point
├── pipeline/                  # Targets for processing, models, and outputs
├── R/
│   ├── 00_main_config.R       # Input paths and analysis settings
│   └── functions/             # Reusable analysis and output functions
├── tests/testthat/            # Small tests of data and analysis logic
├── reports/
│   └── analysis_outputs.qmd   # Quarto summary of key results
├── renv.lock                  # Recorded R package versions
├── inputs/                    # External inputs and derived data (Git-ignored)
├── outputs/                   # Models, tables, figures, diagnostics (Git-ignored)
└── _targets/                  # Pipeline cache (Git-ignored)
```

`_targets.R` defines the workflow; the files in `pipeline/` keep each stage
readable, and `R/functions/` contains the underlying methods.

## Required data

External data are **not distributed in this repository**. Supply the following
datasets at the paths specified in `R/00_main_config.R` before running the
pipeline. Shapefiles also need their accompanying `.shx`, `.dbf`, and `.prj`
files.

| Dataset | Expected path |
| --- | --- |
| Atlas of Living Australia occurrence records | `inputs/data_raw/species/ala_df_raw.csv` |
| Species traits | `inputs/data_raw/species/trait_set.csv` |
| 2.5 km hexagon grid | `inputs/data_raw/spatial_grids/hex_2_5.shp` |
| IBRA7 subregions | `inputs/data_raw/spatial_grids/IBRA7_subregions/IBRA7_subregions.shp` |
| Processed landscape metrics | `inputs/data_raw/lsm_processed/selected_metrics_2_5km.rds` |

To reproduce the reported results, use the same versions of these datasets.
Their provenance and access conditions should be supplied alongside the public
release.

## Run the workflow

Run commands from the repository root with **R 4.6.0** and the Quarto command
line tool available. The package versions used for this workflow are recorded
in `renv.lock`. Install `renv` if necessary, restore the packages, and load the
project library for this R session:

```r
install.packages("renv")
renv::restore(project = ".")
renv::load(project = ".")
```

The repository does not currently auto-activate `renv`, so call `renv::load()`
when starting a new R session in this project. Then run the tests and pipeline:

```r
testthat::test_dir(file.path("tests", "testthat"))
targets::tar_outdated()  # Review what will run
targets::tar_make()      # Build the complete analysis and report
```

The LCBD ensemble and GLLVM fits are computationally intensive and can take
hours. `targets` keeps completed steps in `_targets/` and only rebuilds a step
when its inputs, code, or settings change. A named target can be built on its
own, together with any missing upstream steps; for example,
`targets::tar_make(names = lcbd_results_file)`.

## Find the results

| Location | Contents |
| --- | --- |
| `inputs/data_processed/` | Cleaned occurrences, completeness estimates, species matrix, and aligned model data. |
| `outputs/models/` | Saved LCBD, regional, and GLLVM model checkpoints. |
| `outputs/tables/` | Attrition summaries, analysis results, and model coefficients. |
| `outputs/figures/lcbd_plots/` | LCBD response curves. |
| `outputs/figures/regional_plots/` | Regional convergence panels. |
| `outputs/figures/gllvm_plots/` | Species and trait response figures. |
| `outputs/diagnostics/` | GLLVM diagnostic PDF. |
| `reports/analysis_outputs.html` | Concise rendered overview of key results. |

Generated data, models, figures, and the `_targets/` cache are Git-ignored.
The figures in `outputs/figures/` are the publication exports produced by this
workflow; the Quarto report displays a selection of them. The tracked source
code and `renv.lock` allow the workflow to be run with the required external
data.
