# =============================================================================
# Sensitivity settings. Source after R/00_main_config.R.
# Defines settings only; the main configuration and its outputs are unchanged.
# =============================================================================

sensitivity_cfg <- local({
  output_dir <- file.path("outputs", "sensitivity")
  n_lcbd_allocations <- 30L
  list(
    # A separate targets store prevents sensitivity runs from replacing main fits.
    script = "_targets_sensitivity.R",
    store = "_targets_sensitivity",
    output_dir = output_dir,
    model_dir = file.path(output_dir, "models"),
    table_dir = file.path(output_dir, "tables"),
    figure_dir = file.path(output_dir, "figures"),
    diagnostic_dir = file.path(output_dir, "diagnostics"),
    # Change completeness alone. Monthly incidence, minimum effort/richness,
    # urban filtering, pool rules and model formulae remain those of the main run.
    completeness = c(0.75, 0.80, 0.85, 0.90),
    minimum_months = 12L,
    alternative_pool_counts = c(20L, 30L),
    structured_protocol = "community_complete_or_structured",
    # Read only needed columns from the large stage-02 CSV in bounded chunks.
    occurrence_chunk_size = 250000L,

    # Use a matched 30-allocation ensemble for every LCBD sensitivity scenario,
    # including a refit of the main 70% setting. The publication analysis keeps
    # its full 100-allocation ensemble in the separate main pipeline.
    n_lcbd_allocations = n_lcbd_allocations,
    allocation_seeds = cfg$lcbd_allocation_seed_start +
      seq_len(n_lcbd_allocations) - 1L,

    # Residualisation uses broad IBRA regions, as in the original 05a script.
    miner_region_field = "REG_NAME_7",
    miner_region_crs = cfg$lcbd_reference_pool_crs,

    # Reproducible Dunn-Smyth residuals and native GLLVM Q-Q envelopes.
    residual_seed = 4208L,
    qq_envelope_simulations = 150L,
    distance_breaks_m = seq(0, 25000, by = 2500),
    # Preserve the original GLLVM correlogram: row-standardised weights and
    # a one-sided test for positive residual autocorrelation, BH within each band.
    gllvm_moran_style = "W",
    gllvm_moran_alternative = "greater",
    # Original LCBD diagnostic: symmetric binary weights and two-sided tests.
    lcbd_moran_style = "B",
    lcbd_moran_alternative = "two.sided"
  )
})
