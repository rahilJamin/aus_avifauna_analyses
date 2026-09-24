#' Write compact sensitivity tables and figures from completed targets
#'
#' This function reads no raw inputs and fits no model. It collects the compact
#' branch results written by the sensitivity pipeline and makes the files used by
#' the appendix. The main publication-output functions are not called or changed.
#' @param baseline Main LCBD baseline bundle.
#' @param lcbd_results List of scenario results.
#' @param regional_results List of regional scenario results.
#' @param miner_results Noisy Miner residualisation comparison.
#' @param lcbd_spatial Spatial-control diagnostic.
#' @param gllvm_diagnostic_files Diagnostic RDS/SVG paths.
#' @param settings Sensitivity settings.
#' @return Paths to all generated compact exports.
produce_sensitivity_outputs <- function(baseline, lcbd_results, regional_results,
                                        miner_results, lcbd_spatial,
                                        gllvm_diagnostic_files, settings) {
  dir.create(settings$table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(settings$figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(settings$model_dir, recursive = TRUE, showWarnings = FALSE)
  # Derive this output-only path locally. Adding it to the shared sensitivity
  # configuration would unnecessarily invalidate fitted sensitivity targets.
  appendix_dir <- file.path(settings$output_dir, "appendix")
  dir.create(appendix_dir, recursive = TRUE, showWarnings = FALSE)
  write_table <- function(data, filename) {
    path <- sensitivity_path(settings$table_dir, filename)
    readr::write_csv(data, path)
    path
  }
  bind_part <- function(values, name) {
    dplyr::bind_rows(lapply(values, function(x) x[[name]]))
  }

  # `main_70` is one of the fitted sensitivity branches, so every LCBD row in
  # these outputs represents the same reduced allocation design.
  manifest <- bind_part(lcbd_results, "status")
  support <- bind_part(lcbd_results, "support")
  endpoint_results <- bind_part(lcbd_results, "results")
  curves <- bind_part(lcbd_results, "curves")
  main_regional_contrasts <- transform(
    baseline$regional$endpoint_contrasts,
    scenario_id = "main_70"
  )
  main_regional_status <- data.frame(
    gradient = unique(main_regional_contrasts$gradient),
    status = "complete",
    reason = "Completed main analysis",
    scenario_id = "main_70"
  )
  main_regional_definitions <- transform(
    baseline$regional$settings$gradient_definitions,
    scenario_id = "main_70"
  )
  regional_status <- dplyr::bind_rows(
    main_regional_status, bind_part(regional_results, "status")
  )
  regional_contrasts <- dplyr::bind_rows(
    main_regional_contrasts, bind_part(regional_results, "contrasts")
  )
  regional_definitions <- dplyr::bind_rows(
    main_regional_definitions, bind_part(regional_results, "definitions")
  )
  if (!nrow(regional_contrasts)) {
    regional_contrasts <- data.frame(
      scenario_id = character(), gradient = character(), gradient_label = character(),
      component = character(), high_minus_low = numeric(),
      difference_lower_95 = numeric(), difference_upper_95 = numeric(),
      bootstrap_two_sided_p_value = numeric()
    )
  }

  # GLLVM diagnostics are already compact RDS files. Avoid reopening their
  # source models while producing the manuscript tables.
  diag_rds <- unlist(gllvm_diagnostic_files, use.names = FALSE)
  diag_rds <- diag_rds[grepl("[.]rds$", diag_rds)]
  diagnostics <- lapply(diag_rds, readRDS)
  gllvm_summary <- dplyr::bind_rows(lapply(diagnostics, `[[`, "summary"))
  gllvm_moran <- dplyr::bind_rows(lapply(diagnostics, `[[`, "moran"))
  gllvm_coefficients <- dplyr::bind_rows(lapply(diagnostics, `[[`, "coefficients"))
  # Pair the environmental coefficients by term only after verifying the same
  # response in the model diagnostic. This describes coefficient movement, not
  # a formal species-level difference test.
  coefficient_change <- data.frame(
    term = character(), estimate_with_dbmem = numeric(), se_with_dbmem = numeric(),
    estimate_without_dbmem = numeric(), se_without_dbmem = numeric(),
    estimate_change = numeric()
  )
  if (all(c("environment_with_dbmem", "environment_without_dbmem") %in% gllvm_coefficients$model)) {
    with_mem <- gllvm_coefficients[gllvm_coefficients$model == "environment_with_dbmem", c("term", "estimate", "se")]
    without_mem <- gllvm_coefficients[gllvm_coefficients$model == "environment_without_dbmem", c("term", "estimate", "se")]
    names(with_mem)[-1] <- c("estimate_with_dbmem", "se_with_dbmem")
    names(without_mem)[-1] <- c("estimate_without_dbmem", "se_without_dbmem")
    coefficient_change <- dplyr::mutate(dplyr::inner_join(with_mem, without_mem, by = "term"),
      estimate_change = estimate_with_dbmem - estimate_without_dbmem)
  }

  # Build the two compact tables that are cited in the supplementary text.
  # These are exported as SVGs so they can be inserted directly into Word.
  scenario_table <- dplyr::transmute(
    manifest,
    Scenario = scenario_label,
    `Chao2 sites` = n_chao_sites,
    `LCBD sites` = n_model_sites,
    Species = n_response_species,
    `LCBD fit` = ifelse(status == "complete", "Completed", "Insufficient pool support")
  )
  gllvm_counts <- dplyr::summarise(
    dplyr::group_by(gllvm_moran, model),
    n_tests = sum(status == "ok"),
    n_detected = sum(status == "ok" & moran_i > 0 & p_bh < 0.05, na.rm = TRUE),
    n_species_detected = dplyr::n_distinct(
      series[status == "ok" & moran_i > 0 & !is.na(p_bh) & p_bh < 0.05]
    ),
    .groups = "drop"
  )
  gllvm_overview <- dplyr::left_join(gllvm_summary, gllvm_counts, by = "model")
  model_levels <- c("environment_with_dbmem", "environment_without_dbmem")
  model_labels <- c(
    environment_with_dbmem = "Environmental, with dbMEMs",
    environment_without_dbmem = "Environmental, without dbMEMs"
  )
  gllvm_overview <- gllvm_overview[match(model_levels, gllvm_overview$model), ]
  gllvm_overview <- gllvm_overview[!is.na(gllvm_overview$model), ]
  gllvm_table <- dplyr::transmute(
    gllvm_overview,
    Model = unname(model_labels[model]),
    Sites = n_sites,
    Species = n_species,
    dbMEMs = n_dbmem,
    AIC = round(aic),
    `Species with detected band` = paste0(n_species_detected, "/", n_species),
    `Detected species-band tests` = paste0(n_detected, "/", n_tests)
  )

  files <- c(
    write_table(manifest, "sensitivity_scenario_manifest.csv"),
    write_table(support, "sensitivity_lcbd_pool_support.csv"),
    write_table(endpoint_results, "sensitivity_lcbd_endpoint_results.csv"),
    write_table(curves, "sensitivity_lcbd_curves.csv"),
    write_table(regional_status, "sensitivity_regional_status.csv"),
    write_table(regional_definitions, "sensitivity_regional_thresholds.csv"),
    write_table(regional_contrasts, "sensitivity_regional_contrasts.csv"),
    write_table(miner_results$residualisation$association, "sensitivity_miner_residualisation.csv"),
    write_table(miner_results$residualisation$diagnostics,
                "sensitivity_miner_residualisation_diagnostics.csv"),
    write_table(miner_results$comparison$summary, "sensitivity_miner_lcbd_results.csv"),
    write_table(lcbd_spatial$summary, "sensitivity_lcbd_spatial_models.csv"),
    write_table(lcbd_spatial$effects, "sensitivity_lcbd_spatial_effects.csv"),
    write_table(lcbd_spatial$moran, "sensitivity_lcbd_residual_moran.csv"),
    write_table(gllvm_summary, "sensitivity_gllvm_diagnostics.csv"),
    write_table(gllvm_moran, "sensitivity_gllvm_residual_moran.csv"),
    write_table(coefficient_change, "sensitivity_gllvm_dbmem_coefficient_change.csv")
  )

  table_s1 <- sensitivity_path(
    appendix_dir, "01_table_s1_sampling_and_reference_pool_scenarios.svg"
  )
  save_sensitivity_table_svg(
    scenario_table,
    table_s1,
    title = "Table S1. Sampling and reference-pool sensitivity scenarios",
    note = paste(
      "Chao2 sites satisfy each scenario's sampling and completeness criteria;",
      "LCBD sites are candidates before gradient-specific pool filtering."
    ),
    width = 11
  )
  files <- c(files, table_s1)

  # Report complete smooths: endpoint differences cannot describe a U-shaped
  # response. Keep endpoint CSVs for inspection, without using them as a test
  # of whether the whole smooth is associated with its predictor.
  curve_plots <- plot_sensitivity_lcbd_curves(curves, manifest)
  lcbd_filenames <- c(
    agriculture = "02_figure_s1_lcbd_agriculture.svg",
    patch_area = "03_figure_s2_lcbd_patch_area.svg",
    cohesion = "04_figure_s3_lcbd_cohesion.svg",
    noisy_miner = "05_figure_s4_lcbd_noisy_miner.svg"
  )
  for (gradient in names(curve_plots)) {
    if (!gradient %in% names(lcbd_filenames)) {
      stop("No supplementary figure label is defined for gradient: ", gradient,
           call. = FALSE)
    }
    lcbd_figure <- sensitivity_path(appendix_dir, lcbd_filenames[[gradient]])
    ggplot2::ggsave(
      lcbd_figure,
      curve_plots[[gradient]],
      device = grDevices::svg,
      width = 11,
      height = 8,
      units = "in",
      bg = "white"
    )
    files <- c(files, lcbd_figure)
  }
  if (nrow(miner_results$comparison$curves)) {
    curve_mean <- dplyr::summarise(dplyr::group_by(miner_results$comparison$curves,
      scenario, lcbd_component, gradient_position_percent),
      mean_partial_effect = mean(partial_effect), .groups = "drop")
    miner_plot <- ggplot2::ggplot(curve_mean,
      ggplot2::aes(gradient_position_percent, mean_partial_effect, colour = scenario)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.35) +
      ggplot2::geom_line(linewidth = 0.8) + ggplot2::facet_wrap(~ lcbd_component, scales = "free_y") +
      ggplot2::scale_colour_manual(values = c("Raw detection rate" = "#C44E52",
        "Landscape-residualised detection" = "#2A7F9E")) +
      ggplot2::labs(x = "Position in each central predictor range (%)",
        y = "Partial effect on pool-size-standardised LCBD", colour = NULL) +
      ggplot2::theme_classic(base_size = 10)
    miner_figure <- sensitivity_path(
      appendix_dir, "06_figure_s5_noisy_miner_residualisation.svg"
    )
    ggplot2::ggsave(
      miner_figure,
      miner_plot,
      device = grDevices::svg,
      width = 8.5,
      height = 3.5,
      units = "in",
      bg = "white"
    )
    files <- c(files, miner_figure)
  }
  if (nrow(regional_contrasts)) {
    regional_plot <- ggplot2::ggplot(regional_contrasts,
      ggplot2::aes(scenario_id, high_minus_low, colour = component)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.35) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = difference_lower_95,
        ymax = difference_upper_95), width = 0.15) + ggplot2::geom_point() +
      ggplot2::facet_wrap(~ gradient_label, scales = "free_y") +
      ggplot2::labs(x = NULL, y = "High-minus-low between-subregion dissimilarity",
                    colour = "Component") + ggplot2::theme_classic(base_size = 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    regional_figure <- sensitivity_path(
      appendix_dir, "07_figure_s6_regional_convergence.svg"
    )
    ggplot2::ggsave(
      regional_figure,
      regional_plot,
      device = grDevices::svg,
      width = 11,
      height = 7,
      units = "in",
      bg = "white"
    )
    files <- c(files, regional_figure)
  }
  # Plot the existing Moran diagnostics; no residuals or tests are recomputed.
  # Keep these alongside the package-native diagnostic figures.
  dir.create(settings$diagnostic_dir, recursive = TRUE, showWarnings = FALSE)
  spatial_plots <- plot_sensitivity_spatial_diagnostics(gllvm_moran, lcbd_spatial$moran)
  spatial_filenames <- c(
    lcbd = "08_figure_s7_lcbd_spatial_autocorrelation.svg",
    gllvm = "09_figure_s8_gllvm_spatial_autocorrelation.svg"
  )
  for (name in names(spatial_plots)) {
    spatial_figure <- sensitivity_path(appendix_dir, spatial_filenames[[name]])
    ggplot2::ggsave(spatial_figure, spatial_plots[[name]],
      device = grDevices::svg, width = 11,
      height = if (name == "lcbd") 7 else 4.5, units = "in", bg = "white")
    files <- c(files, spatial_figure)
  }

  table_s2 <- sensitivity_path(
    appendix_dir, "10_table_s2_gllvm_spatial_diagnostics.svg"
  )
  save_sensitivity_table_svg(
    gllvm_table,
    table_s2,
    title = "Table S2. Environmental GLLVM spatial diagnostics",
    note = paste(
      "Detection is positive Moran's I with within-band BH-adjusted p < 0.05.",
      "AIC is comparable only between the matched environmental models."
    ),
    width = 12
  )
  files <- c(files, table_s2)

  expected_appendix_files <- vapply(
    c(
      "01_table_s1_sampling_and_reference_pool_scenarios.svg",
      "02_figure_s1_lcbd_agriculture.svg",
      "03_figure_s2_lcbd_patch_area.svg",
      "04_figure_s3_lcbd_cohesion.svg",
      "05_figure_s4_lcbd_noisy_miner.svg",
      "06_figure_s5_noisy_miner_residualisation.svg",
      "07_figure_s6_regional_convergence.svg",
      "08_figure_s7_lcbd_spatial_autocorrelation.svg",
      "09_figure_s8_gllvm_spatial_autocorrelation.svg",
      "10_table_s2_gllvm_spatial_diagnostics.svg"
    ),
    function(filename) sensitivity_path(appendix_dir, filename),
    character(1)
  )
  missing_appendix_files <- expected_appendix_files[!file.exists(expected_appendix_files)]
  if (length(missing_appendix_files)) {
    stop(
      "The following supplementary exports were not created: ",
      paste(basename(missing_appendix_files), collapse = ", "),
      call. = FALSE
    )
  }
  # Store compact results outside the targets cache so a cleared cache does not
  # discard the evidence used by the appendix. The fitted residualisation model
  # stays in the cache: the exported checkpoint contains only its diagnostics
  # and LCBD summaries, never a duplicate fitted model.
  checkpoint <- sensitivity_path(settings$model_dir, "sensitivity_summary.rds")
  saveRDS(list(manifest = manifest, lcbd_results = endpoint_results,
    regional_contrasts = regional_contrasts,
    miner = list(
      association = miner_results$residualisation$association,
      diagnostics = miner_results$residualisation$diagnostics,
      lcbd_summary = miner_results$comparison$summary,
      pool_support = miner_results$comparison$support
    ),
    lcbd_spatial = lcbd_spatial, gllvm = list(summary = gllvm_summary,
      moran = gllvm_moran, coefficient_change = coefficient_change)), checkpoint)
  c(files, checkpoint)
}

#' Save a compact supplementary table as an editable SVG
#'
#' The table uses only base grid graphics, which keeps the export independent
#' of HTML/Quarto rendering and avoids adding a package dependency.
#' @param data Data frame to display.
#' @param path Output SVG path.
#' @param title Standalone table label and title.
#' @param note Brief explanatory note shown below the table.
#' @param width SVG width in inches.
#' @return The output path, invisibly.
save_sensitivity_table_svg <- function(data, path, title, note, width = 11) {
  if (!is.data.frame(data) || !nrow(data) || !ncol(data)) {
    stop("Supplementary SVG tables require a non-empty data frame.", call. = FALSE)
  }
  display <- as.data.frame(lapply(data, as.character), check.names = FALSE)
  row_height <- 0.36
  height <- max(3.2, 1.45 + row_height * (nrow(display) + 1L))
  character_widths <- vapply(seq_along(display), function(index) {
    max(nchar(c(names(display)[index], display[[index]])), na.rm = TRUE)
  }, numeric(1))
  character_widths <- pmax(8, pmin(character_widths, 34))

  grDevices::svg(path, width = width, height = height, pointsize = 10)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
  grid::grid.text(
    title,
    x = grid::unit(0.025, "npc"),
    y = grid::unit(0.95, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontface = "bold", fontsize = 12)
  )

  table_top <- 0.84
  table_bottom <- 0.18
  table_height <- table_top - table_bottom
  layout <- grid::grid.layout(
    nrow = nrow(display) + 1L,
    ncol = ncol(display),
    widths = grid::unit(character_widths, "null"),
    heights = grid::unit(rep(1, nrow(display) + 1L), "null")
  )
  grid::pushViewport(grid::viewport(
    x = 0.5,
    y = (table_top + table_bottom) / 2,
    width = 0.95,
    height = table_height,
    layout = layout
  ))
  for (column in seq_len(ncol(display))) {
    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = column))
    grid::grid.rect(gp = grid::gpar(fill = "#E8E8E8", col = "white"))
    grid::grid.text(
      names(display)[column],
      x = grid::unit(0.04, "npc"),
      just = "left",
      gp = grid::gpar(fontface = "bold", fontsize = 9)
    )
    grid::popViewport()
  }
  for (row in seq_len(nrow(display))) {
    fill <- if (row %% 2L) "white" else "#F7F7F7"
    for (column in seq_len(ncol(display))) {
      grid::pushViewport(grid::viewport(
        layout.pos.row = row + 1L,
        layout.pos.col = column
      ))
      grid::grid.rect(gp = grid::gpar(fill = fill, col = "#DDDDDD", lwd = 0.35))
      grid::grid.text(
        display[row, column],
        x = grid::unit(0.04, "npc"),
        just = "left",
        gp = grid::gpar(fontsize = 8.5)
      )
      grid::popViewport()
    }
  }
  grid::popViewport()
  grid::grid.text(
    note,
    x = grid::unit(0.025, "npc"),
    y = grid::unit(0.08, "npc"),
    just = c("left", "bottom"),
    gp = grid::gpar(fontsize = 8.5, col = "#333333")
  )
  invisible(path)
}

#' Plot complete LCBD sensitivity curves from saved ensemble summaries
#'
#' Each figure represents one environmental gradient. Scenario panels retain
#' every fitted grid point and its existing pointwise uncertainty interval.
#' No models are fitted and no effects or intervals are re-estimated.
#' @param curves Combined scenario-level ensemble curves.
#' @param manifest Scenario labels and completion status.
#' @return Named list of ggplot figures, one per focal gradient.
plot_sensitivity_lcbd_curves <- function(curves, manifest) {
  if (!nrow(curves)) return(list())
  supported <- manifest[manifest$status == "complete", ]
  if (any(!curves$scenario_id %in% supported$scenario_id)) {
    stop("LCBD curves must belong to completed sensitivity scenarios.", call. = FALSE)
  }
  curves$scenario_label <- factor(curves$scenario_id,
    levels = supported$scenario_id, labels = supported$scenario_label)
  curves$lcbd_component <- factor(curves$lcbd_component,
    levels = c("total", "replacement", "richness_difference"),
    labels = c("Total", "Replacement", "Richness difference"))
  colours <- c(Total = "#333333", Replacement = "#94346E",
    `Richness difference` = "#2A7F9E")

  # Standardisation and the central predictor range can differ by scenario.
  # Retain the saved x values and show each scenario separately; do not imply
  # that an x position denotes identical raw cover/rate across these panels.
  plots <- lapply(split(curves, curves$focal_gradient), function(dat) {
    ggplot2::ggplot(dat, ggplot2::aes(gradient_value, mean_partial_effect,
      colour = lcbd_component, fill = lcbd_component)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower_95, ymax = upper_95),
        alpha = 0.10, colour = NA) +
      ggplot2::geom_line(linewidth = 0.65) +
      ggplot2::facet_wrap(~ scenario_label, ncol = 3, scales = "free") +
      ggplot2::scale_colour_manual(values = colours, drop = FALSE) +
      ggplot2::scale_fill_manual(values = colours, drop = FALSE) +
      ggplot2::labs(title = unique(dat$focal_label),
        x = "Predictor (scenario-specific standardised units)",
        y = "Partial effect on pool-size-standardised LCBD",
        colour = NULL, fill = NULL) +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(legend.position = "bottom",
        panel.grid.minor = ggplot2::element_blank(),
        strip.background = ggplot2::element_rect(fill = "grey98"))
  })
  plots
}

#' Compare saved residual spatial diagnostics across model specifications
#'
#' GLLVM boxes show variation among species, not uncertainty of a pooled mean.
#' LCBD lines retain signed Moran statistics and display the saved adjusted
#' tests for the representative allocation. Failed tests are gaps, never zeros.
#' @param gllvm_moran Saved species-by-distance-band GLLVM diagnostics.
#' @param lcbd_moran Saved LCBD component-by-distance-band diagnostics.
#' @return Named list of GLLVM and LCBD ggplot figures.
plot_sensitivity_spatial_diagnostics <- function(gllvm_moran, lcbd_moran) {
  plots <- list()
  if (nrow(gllvm_moran)) {
    # The dbMEM sensitivity is a matched comparison of the two environmental
    # models. Fourth-corner diagnostics are retained in the diagnostic files,
    # but are not mixed into this comparison figure or its summary table.
    dat <- gllvm_moran[gllvm_moran$model %in% c(
      "environment_without_dbmem", "environment_with_dbmem"
    ), , drop = FALSE]
    dat$moran_i[dat$status != "ok"] <- NA_real_
    dat$model <- factor(dat$model,
      levels = c("environment_without_dbmem", "environment_with_dbmem"),
      labels = c("Environmental: without dbMEMs", "Environmental: with dbMEMs"))
    # Plot the stored expectation as well as the observed statistic. The
    # expectation depends on the number of connected sites in each band.
    expected <- unique(dat[c("model", "midpoint_km", "expected_i")])
    plots$gllvm <- ggplot2::ggplot(dat,
      ggplot2::aes(midpoint_km, moran_i, group = midpoint_km)) +
      ggplot2::geom_boxplot(width = 1.7, fill = "#C9DFDE", colour = "#2F6F73",
        outlier.size = 0.6, linewidth = 0.35, na.rm = TRUE) +
      ggplot2::geom_line(data = expected,
        ggplot2::aes(midpoint_km, expected_i, group = 1), inherit.aes = FALSE,
        linetype = "dashed", colour = "grey35", linewidth = 0.45) +
      ggplot2::facet_wrap(~ model, nrow = 1) +
      ggplot2::scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25)) +
      ggplot2::labs(x = "Distance-band midpoint (km)", y = "Residual Moran's I") +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  }
  if (nrow(lcbd_moran)) {
    dat <- lcbd_moran
    dat$moran_i[dat$status != "ok"] <- NA_real_
    dat$model <- factor(dat$model, levels = c("without_gp", "with_gp"),
      labels = c("Without spatial smooth", "With spatial smooth"))
    dat$component <- factor(dat$series,
      levels = c("total", "replacement", "richness_difference"),
      labels = c("Total", "Replacement", "Richness difference"))
    dat$gradient <- factor(dat$focal_gradient,
      levels = c("agriculture", "patch_area", "cohesion", "noisy_miner"),
      labels = c("Agricultural cover", "Mean patch area", "Woodland cohesion", "Noisy Miner"))
    dat$detection <- ifelse(is.na(dat$p_bh), NA_character_,
      ifelse(dat$p_bh < 0.05, "BH p < 0.05", "BH p >= 0.05"))
    expected <- unique(dat[c("gradient", "component", "midpoint_km", "expected_i")])
    plots$lcbd <- ggplot2::ggplot(dat,
      ggplot2::aes(midpoint_km, moran_i, colour = model, group = model)) +
      ggplot2::geom_line(data = expected,
        ggplot2::aes(midpoint_km, expected_i, group = 1), inherit.aes = FALSE,
        linetype = "dashed", colour = "grey50", linewidth = 0.35) +
      ggplot2::geom_line(linewidth = 0.6, na.rm = TRUE) +
      ggplot2::geom_point(ggplot2::aes(shape = detection), size = 1.8, na.rm = TRUE) +
      ggplot2::facet_grid(component ~ gradient, scales = "free_y") +
      ggplot2::scale_colour_manual(values = c("Without spatial smooth" = "#C44E52",
        "With spatial smooth" = "#2A7F9E")) +
      ggplot2::scale_shape_manual(values = c("BH p < 0.05" = 16, "BH p >= 0.05" = 1)) +
      ggplot2::scale_x_continuous(breaks = c(0, 10, 20)) +
      ggplot2::labs(x = "Distance-band midpoint (km)", y = "Residual Moran's I",
        colour = NULL, shape = NULL) +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(legend.position = "bottom",
        panel.grid.minor = ggplot2::element_blank())
  }
  plots
}
