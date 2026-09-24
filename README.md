# GEB: Australian bird-community analysis

This repository contains the code used to analyse Australian bird communities
recorded from **1 January 2020 to 28 February 2026**. Using hexagons with 2.5 km
opposite-edge spacing in Australian Albers as sampling sites, the workflow
examines how local community uniqueness and
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
for its own analysis.

## Repository structure

```text
geb_repo/
├── _targets.R                 # Pipeline entry point
├── _targets_sensitivity.R     # Separate supplementary sensitivity pipeline
├── _targets_demo.R            # Isolated quick-run pipeline
├── DESCRIPTION                # Project metadata and direct R dependencies
├── pipeline/                  # Targets for processing, models, and outputs
├── R/
│   ├── 00_main_config.R       # Input paths and analysis settings
│   ├── 01_sensitivity_config.R # Sensitivity settings and isolated paths
│   ├── 02_demo_config.R        # Demo-only input/output paths and run settings
│   ├── functions/              # Reusable analysis and output functions
│   └── sensitivity/            # Supplementary-analysis functions
├── scripts/                   # Data-acquisition reference and sensitivity runner
├── tests/testthat/            # Small tests of data and analysis logic
├── demo/
│   └── data/                  # Tracked inputs for the quick-run demonstration
├── metadata/
│   ├── data_manifest.json      # Input/output inventory and workflow references
│   └── inputs/                 # Field dictionaries and source class legends
├── renv.lock                  # Recorded R package versions
├── inputs/
│   ├── data_raw/              # External data (Git-ignored)
│   └── data_processed/        # Generated data (Git-ignored)
├── outputs/                   # Models, tables, figures, diagnostics (Git-ignored)
├── _targets/                  # Main pipeline cache (Git-ignored)
└── _targets_demo/             # Demo pipeline cache (Git-ignored)
```

`_targets.R` defines the workflow; the files in `pipeline/` keep each stage
readable, and `R/functions/` contains the underlying methods.

## Required data

The full analysis inputs are **not distributed in this repository**. Supply the
following datasets at the paths specified in `R/00_main_config.R` before running
the main pipeline. Shapefiles also need their accompanying `.shx`, `.dbf`, and
`.prj` files. A reduced input set for the quick demo is included under
`demo/data/inputs/`.

The machine-readable [`metadata/data_manifest.json`](metadata/data_manifest.json)
indexes each input, its pipeline role, and the corresponding detailed metadata
file. The input-specific JSON files define the columns and units, spatial
attributes, source raster codes, and the exact reclassification used by the
workflow. Fields that still need author confirmation are labelled in those
files; they must be resolved before the Zenodo release.

The workflow starts from these frozen external inputs. It creates the hexagon
grid, cleaned records, retained site identifiers and landscape metrics during
execution; none of those derived objects needs to be supplied.

| Dataset | Expected path | Field and source metadata |
| --- | --- | --- |
| Atlas of Living Australia occurrence records | `inputs/data_raw/species/ala_df_raw.csv` | [ALA occurrence fields](metadata/inputs/ala_occurrences.json) |
| AVONET-derived species traits | `inputs/data_raw/species/trait_set.csv` | [Trait fields and model transformations](metadata/inputs/avonet_traits.json) |
| IBRA7 subregions | `inputs/data_raw/spatial_grids/IBRA7_subregions/ibra7_subregions.shp` | [IBRA attributes](metadata/inputs/ibra7_subregions.json) |
| Analysis study area | `inputs/data_raw/spatial_grids/study_area.shp` | [Study-area metadata](metadata/inputs/study_area.json) |
| CLUM land-use raster and sidecars | `inputs/data_raw/raw_rasters/clum_50m_2023_v2.tif` | [CLUM raster fields and SIMPN legend](metadata/inputs/clum_landuse.json) |
| NVIS major-vegetation-group raster and sidecars | `inputs/data_raw/raw_rasters/nvis_mvg.tif` | [NVIS raster fields and reclassification](metadata/inputs/nvis_vegetation.json) |

To reproduce the reported results, use the same versions of these datasets.
Their provenance and access conditions should be supplied alongside the public
release.

Use the frozen ALA snapshot for reproduction. Data acquisition is separate from
the analysis: the targets pipeline does not submit a new ALA query or download
rasters. The query reference under `scripts/data_acquisition/` documents the
available acquisition code and is not sourced by `_targets.R`. It uses the
optional `galah` package, which can be installed separately when a new ALA
download is required; `galah` is not needed to reproduce the frozen analysis.

## Run the workflow

Run commands from the repository root with **R 4.6.0**. The package versions
used for this workflow are recorded in `renv.lock`. Install `renv` if necessary,
restore the packages, and restart R from the repository root:

```r
install.packages("renv")
renv::restore(project = ".")
```

The committed `.Rprofile` automatically activates the project library when R
starts from the repository root. Then run the tests and pipeline:

```r
testthat::test_dir(file.path("tests", "testthat"))
targets::tar_outdated()  # Review what will run
targets::tar_make()      # Build the complete analysis and publication outputs
```

The LCBD ensemble and GLLVM fits are computationally intensive and can take
hours. `targets` keeps completed steps in `_targets/` and only rebuilds a step
when its inputs, code, or settings change. A named target can be built on its
own, together with any missing upstream steps; for example,
`targets::tar_make(names = lcbd_results_file)`.

## Run the quick demo

The demo runs the shared cleaning, Chao2 filtering, raster reclassification,
landscape-metric extraction, LCBD, regional-convergence and GLLVM code on a
small spatially coherent sample. It includes cropped source-code CLUM and NVIS
rasters, so collaborators can exercise raster processing without downloading
or processing the full-country rasters. Demo inputs and outputs have their own
paths and cache; they cannot replace the main workflow's files.

The tracked demo inputs contain 150 selected sites in a compact extent spanning
New South Wales and Victoria and occurrence records for 50 common,
trait-matched species, including Noisy Miner. The selection and provenance are
recorded in `demo/data/demo_input_manifest.json`. These subsets retain the
source licences and attribution requirements documented for the corresponding
full inputs under `metadata/inputs/`.

After restoring the packages recorded in `renv.lock`, run the demo from the
repository root:

```r
targets::tar_manifest(script = "_targets_demo.R")
targets::tar_outdated(script = "_targets_demo.R", store = "_targets_demo")
targets::tar_make(
  script = "_targets_demo.R",
  store = "_targets_demo"
)
```

The main occurrence date, protocol, coordinate, Chao2, completeness and
occupancy filters are preserved. GLLVM formulas, spatial selection and
initialisation settings are preserved. Regional endpoint support is reduced to
one site in at least two subregions for the small dataset. The demo LCBD ensemble uses three
allocations of three reference pools rather than the main analysis's 100
allocations of 25 pools, with smaller pool-support minimums appropriate for the
107-site demo; its outputs are illustrative and not publication
estimates. The demo pipeline tests the main analysis and plotting code, stores
model checkpoints under `demo/outputs/models/`, and writes figures beneath
`demo/outputs/figures/`. Demo outputs and the `_targets_demo/` cache are
generated locally and are not committed.

## Run the sensitivity analyses

The sensitivity analyses deliberately use a separate pipeline. They read the
completed main checkpoints but use their own `_targets_sensitivity/` cache and
write only to `outputs/sensitivity/`, so they never overwrite main results.
Keeping this entry point separate also lets collaborators run, clear, or update
supplementary analyses without placing the completed main models at risk.

They assess completeness thresholds of 75%, 80%, 85%, and 90%; a 12-month
minimum sampling rule; daily rather than monthly incidence units; structured
protocols only; 20 and 30 LCBD reference pools; Noisy Miner detection after
landscape residualisation; LCBD spatial control; and GLLVM diagnostics with a
matched environmental model without dbMEM predictors. The supplementary dbMEM
comparison uses only the two matched environmental GLLVMs; fourth-corner
diagnostics remain separate from that comparison.

The LCBD sensitivity scenarios use 30 matched allocation seeds. The main 70%
setting is refitted with those same 30 seeds, so comparisons do not mix a
reduced sensitivity ensemble with the 100-allocation publication ensemble.
The latter remains unchanged in the main pipeline. From the repository root,
inspect the sensitivity plan and then run a named section or the full
supplementary workflow:

```r
source(file.path("scripts", "run_sensitivity.R"))

targets::tar_make(
  names = sensitivity_lcbd_result,
  script = "_targets_sensitivity.R",
  store = "_targets_sensitivity"
)

targets::tar_make(
  names = sensitivity_publication_files,
  script = "_targets_sensitivity.R",
  store = "_targets_sensitivity"
)
```

The second command builds any missing sensitivity prerequisites, including the
matched no-dbMEM environmental GLLVM, before exporting numbered supplementary
figures and tables to `outputs/sensitivity/appendix/`.
If all fitted targets are already cached, only this lightweight output target
runs. To run all sensitivity targets in one step, omit `names` from
`targets::tar_make()`. The full command list is also in
`scripts/run_sensitivity.R`.

## How the pipeline is organised

The main dependency chain is:

1. Create the hexagon grid from the study extent in Australian Albers.
2. Clean the 2020-onward ALA records and assign them to that grid.
3. Estimate iNEXT/Chao2 completeness, retain eligible sites and build their
   binary species matrix. Extract landscape metrics at those same sites from
   the reclassified CLUM and NVIS rasters.
4. Align species, sampling, spatial, and landscape data into one model dataset.
5. Fit the LCBD, regional convergence, and GLLVM analyses.
6. Export publication tables, figures and diagnostics.

Grid creation transforms the extent to **GDA94 / Australian Albers (EPSG:3577)**
before applying the 2500 m cell spacing. Cells are point-topped; the grid origin
is the southwest corner of the projected extent's bounding box. Whole cells
intersecting the extent are retained. They have equal projected area of about
5.413 km², including at the boundary. `hex_id` values are assigned before spatial
filtering and remain tied to the generated geometry. Grid settings are declared
under `cfg$hex_grid`; the saved layer is
`inputs/data_processed/spatial_grids/hex_2_5.shp`.

Occurrence longitude/latitude coordinates are transformed to the generated
grid CRS before spatial assignment. The grid file is shared by occurrence
assignment, landscape extraction, LCBD and regional input preparation; GLLVM
uses the resulting aligned model data. Records intersecting more than one
hexagon cause an explicit error, requiring boundary assignments to be reviewed.

Raster preprocessing uses the same study extent, transformed to each source
raster's native CRS before cropping and masking. Classified rasters retain the
configured landscape projection and nearest-neighbour resampling. The same
generated hexagons are transformed to that projection for centroid calculation
and sampling in circular windows of **2500 m radius**. This radius is a separate
setting from hexagon spacing. Classification maps, CRS, buffer size, requested
metrics and output columns are declared under `cfg$landscape`.

The `landscape_site_ids` target comes directly from the current iNEXT/Chao2
retained-site table. Raster reclassification can run independently of occurrence
processing, while metric extraction depends on that table and the generated
grid. There is one date window and one completeness-filtering stage.

The selected landscape table retains complete values across ten predictors.
Model-data preparation applies the urban filter and scales landscape predictors
over that filtered landscape table before aligning its sites with the
Chao2-retained species matrix. Changing the retained sites therefore changes
both site availability and scaling. GLLVM and regional analyses also
require IBRA assignments; their final site counts can differ from LCBD.

The generated landscape table is
`inputs/data_processed/lsm_processed/selected_metrics_2_5km.rds`. It is the only
landscape input used by model-data preparation. To build just this stage and its
required upstream targets, run:

```r
targets::tar_make(names = landscape_metrics_file)
```

To create the grid alone, run:

```r
targets::tar_make(names = hex_grid_files)
```

To build and inspect preprocessing before fitting the heavy models, run:

```r
targets::tar_make(names = c(model_data_files, gllvm_prepared_data))
targets::tar_read(inext_matrix_bundle)$attrition
targets::tar_read(landscape_metric_bundle)$attrition
targets::tar_read(model_data_bundle)$attrition
targets::tar_read(gllvm_prepared_data)$dimensions
```

These targets rebuild missing upstream stages but do not perform dbMEM selection
or fit the LCBD, regional or GLLVM models. After reviewing them, use
`targets::tar_make()` to build the complete analysis and publication outputs.

Changing the extent, grid settings, occurrence filters, Chao2 settings, raster
inputs or landscape settings invalidates the affected downstream targets
automatically. `targets::tar_outdated()` shows what would be rebuilt before a run.

## Find the results

| Location | Contents |
| --- | --- |
| `inputs/data_processed/` | Generated grid, cleaned occurrences, completeness estimates, species matrix, reclassified rasters, landscape metrics and aligned model data. |
| `outputs/models/` | Saved LCBD, regional, and GLLVM model checkpoints. |
| `outputs/tables/` | Attrition summaries, analysis results, and model coefficients. |
| `outputs/figures/lcbd_plots/` | LCBD response curves. |
| `outputs/figures/regional_plots/` | Regional convergence panels. |
| `outputs/figures/gllvm_plots/` | Species and trait response figures. |
| `outputs/diagnostics/` | GLLVM diagnostic PDF. |
| `outputs/sensitivity/tables/` | Machine-readable sensitivity results and diagnostics. |
| `outputs/sensitivity/diagnostics/` | Package-native GLLVM diagnostic files and compact diagnostic objects. |
| `outputs/sensitivity/appendix/` | Numbered, appendix-ready SVG figures and tables. |
| `demo/data/inputs/` | Tracked input subsets for the quick-run demonstration. |
| `demo/outputs/` | Locally generated demo data, models, tables and figures. |

Generated data, models, figures, and the `_targets/` and `_targets_sensitivity/`
caches are Git-ignored.
The figures in `outputs/figures/` are the main publication exports produced by
this workflow. Supplementary figures and tables are numbered in methods order
under `outputs/sensitivity/appendix/` for transfer into the final Word file.
The tracked source code and `renv.lock` allow the workflow to be run with the
required external data.

## Citation

If you use this workflow, code, or demonstration dataset, cite the archived
software release using the metadata in [`CITATION.cff`](CITATION.cff). GitHub's
**Cite this repository** control reads that file. A version-specific Zenodo DOI
will be added after the repository release is archived.

## Licence

The original software in this repository is available under the
[MIT License](LICENSE). The licence permits use, copying, modification, and
redistribution provided that the copyright and licence notice are retained.
Please also cite the archived software release as described above.

External datasets and the demonstration subsets derived from them retain their
source licences and attribution requirements. The MIT License does not
relicense ALA, AVONET, CLUM, NVIS, IBRA7, or other third-party material. Their
provenance and reuse conditions are recorded under [`metadata/inputs/`](metadata/inputs/).
