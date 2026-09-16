#' Produce LCBD tables and publication figure
#'
#' Reads the completed model checkpoint and reproduces the validated
#' publication tables and plots. Plot data, scales, colours, and export
#' settings are carried over from the original output script unchanged.
#'
#' @param cfg Main workflow configuration with input and output paths.
#' @return Paths to the generated files, for a `targets` file target.

import::from(dplyr, arrange, case_when, count, desc, first, group_by, last,
             left_join, mutate, n, recode, slice, summarise, transmute,
             ungroup)
import::from(ggh4x, facet_grid2)
import::from(ggplot2, aes, element_blank, element_line, element_rect,
             element_text, geom_hline, geom_line, ggplot, ggsave, labs, theme,
             theme_bw)
import::from(grDevices, svg)
import::from(magrittr, `%>%`)
import::from(readr, write_csv)
import::from(stats, median, quantile)

produce_lcbd_outputs <- function(cfg) {
  # Produce LCBD tables ---------------------------------------------------------
  # Paths -----------------------------------------------------------------------

  input_path <- cfg$lcbd_results_path

  table_dir <- cfg$table_dir

  main_result_path <- file.path(
    table_dir,
    "04_local_uniqueness_results.csv"
  )
  allocation_validation_path <- file.path(
    table_dir,
    "04_reference_pool_allocations.csv"
  )
  curve_summary_path <- file.path(
    table_dir,
    "04_lcbd_curve_summary.csv"
  )

  if (!file.exists(input_path)) {
    stop(
      "Missing stage-04 result object: ",
      input_path,
      call. = FALSE
    )
  }

  # Validate output destinations before creating folders or exporting.
  invisible(lapply(c(table_dir), assert_not_raw_output_path))
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  # Load the curves and model summaries ----------------------------------------

  lcbd_results <- read_rds_checked(input_path, "saved analysis result", required_class = "list")

  # Require the confirmed Gaussian checkpoint, including when a path is overridden.
  if (!identical(lcbd_results$settings$model_family,
                 "Gaussian identity for all three LCBD components")) {
    stop("Expected the Gaussian LCBD main result checkpoint.", call. = FALSE)
  }

  required_objects <- c(
    "curves_by_allocation",
    "model_results_by_allocation",
    "main_results",
    "allocation_validation"
  )

  missing_objects <- setdiff(required_objects, names(lcbd_results))
  if (length(missing_objects) > 0L) {
    stop(
      "The stage-04 result object lacks: ",
      paste(missing_objects, collapse = ", "),
      call. = FALSE
    )
  }

  allocation_curves <- lcbd_results$curves_by_allocation
  allocation_model_results <- lcbd_results$model_results_by_allocation
  main_results <- lcbd_results$main_results
  allocation_validation <- lcbd_results$allocation_validation

  # Export the primary 12-row result table -------------------------------------

  local_uniqueness_table <- main_results %>%
    transmute(
      focal_gradient = focal_label,
      lcbd_component = recode(
        lcbd_component,
        total = "Total",
        replacement = "Replacement",
        richness_difference = "Richness difference"
      ),
      response = "Raw LCBD multiplied by reference-pool site count",
      model_family = "Gaussian identity",
      n_allocations,
      n_unique_partitions,
      minimum_n_sites,
      mean_n_sites,
      maximum_n_sites,
      minimum_n_informative_pools,
      mean_n_informative_pools,
      maximum_n_informative_pools,
      endpoint_low,
      endpoint_high,
      mean_endpoint_difference,
      lower_95,
      upper_95,
      between_allocation_sd,
      direction_agreement,
      median_edf,
      all_models_converged,
      interpretation
    ) %>%
    arrange(focal_gradient, lcbd_component)

  write_csv(local_uniqueness_table, main_result_path)

  # Export the one-row-per-allocation validation table -------------------------

  write_csv(allocation_validation, allocation_validation_path)

  required_curve_columns <- c(
    "allocation",
    "seed",
    "focal_gradient",
    "focal_label",
    "focal_variable",
    "lcbd_component",
    "gradient_value",
    "partial_effect"
  )

  missing_curve_columns <- setdiff(
    required_curve_columns,
    names(allocation_curves)
  )
  if (length(missing_curve_columns) > 0L) {
    stop(
      "The saved LCBD curves lack: ",
      paste(missing_curve_columns, collapse = ", "),
      call. = FALSE
    )
  }

  # Calculate endpoint change and curve effect range ----------------------------

  # Endpoint change is the partial effect at the high end minus the partial
  # effect at the low end of the common plotted range. It describes direction.
  #
  # Effect range is the maximum minus minimum partial effect anywhere across that
  # range. It describes magnitude and therefore remains informative for U-shaped
  # or hump-shaped curves whose two endpoints may be similar.
  allocation_curve_effects <- allocation_curves %>%
    arrange(
      allocation,
      focal_gradient,
      lcbd_component,
      gradient_value
    ) %>%
    group_by(
      allocation,
      seed,
      focal_gradient,
      focal_label,
      focal_variable,
      lcbd_component
    ) %>%
    summarise(
      gradient_low = first(gradient_value),
      gradient_midpoint = median(gradient_value),
      gradient_high = last(gradient_value),
      partial_effect_low = first(partial_effect),
      partial_effect_midpoint = partial_effect[
        which.min(abs(gradient_value - median(gradient_value)))
      ],
      partial_effect_high = last(partial_effect),
      endpoint_change = partial_effect_high - partial_effect_low,
      curve_effect_range = max(partial_effect) - min(partial_effect),
      gradient_at_curve_minimum = gradient_value[which.min(partial_effect)],
      gradient_at_curve_maximum = gradient_value[which.max(partial_effect)],
      .groups = "drop"
    ) %>%
    mutate(
      curve_shape = case_when(
        partial_effect_low < partial_effect_midpoint &
          partial_effect_midpoint < partial_effect_high ~ "Increasing",
        partial_effect_low > partial_effect_midpoint &
          partial_effect_midpoint > partial_effect_high ~ "Decreasing",
        partial_effect_midpoint < partial_effect_low &
          partial_effect_midpoint < partial_effect_high ~ "U-shaped",
        partial_effect_midpoint > partial_effect_low &
          partial_effect_midpoint > partial_effect_high ~ "Hump-shaped",
        TRUE ~ "Nonlinear"
      )
    )

  # Summarise variation among the 100 geographical allocations -----------------

  curve_summary <- allocation_curve_effects %>%
    group_by(
      focal_gradient,
      focal_label,
      focal_variable,
      lcbd_component
    ) %>%
    summarise(
      n_allocations = n(),
      gradient_low = first(gradient_low),
      gradient_high = first(gradient_high),
      endpoint_change_median = median(endpoint_change),
      endpoint_change_2_5_percent = quantile(
        endpoint_change,
        0.025,
        names = FALSE
      ),
      endpoint_change_97_5_percent = quantile(
        endpoint_change,
        0.975,
        names = FALSE
      ),
      curve_effect_range_median = median(curve_effect_range),
      curve_effect_range_2_5_percent = quantile(
        curve_effect_range,
        0.025,
        names = FALSE
      ),
      curve_effect_range_97_5_percent = quantile(
        curve_effect_range,
        0.975,
        names = FALSE
      ),
      median_gradient_at_curve_minimum = median(
        gradient_at_curve_minimum
      ),
      median_gradient_at_curve_maximum = median(
        gradient_at_curve_maximum
      ),
      .groups = "drop"
    )

  # Identify the most frequent broad curve shape without a custom helper.
  shape_summary <- allocation_curve_effects %>%
    count(
      focal_gradient,
      focal_label,
      focal_variable,
      lcbd_component,
      curve_shape,
      name = "n_allocations_with_shape"
    ) %>%
    group_by(
      focal_gradient,
      focal_label,
      focal_variable,
      lcbd_component
    ) %>%
    mutate(
      shape_agreement = n_allocations_with_shape /
        sum(n_allocations_with_shape)
    ) %>%
    arrange(desc(n_allocations_with_shape), curve_shape) %>%
    slice(1L) %>%
    ungroup() %>%
    transmute(
      focal_gradient,
      focal_label,
      focal_variable,
      lcbd_component,
      dominant_curve_shape = curve_shape,
      shape_agreement
    )

  # EDF is included only to keep this table comparable with the previous curve
  # summary. P-values are already reported by stage 04 and are not duplicated.
  edf_summary <- allocation_model_results %>%
    group_by(
      focal_gradient,
      focal_label,
      focal_variable,
      lcbd_component
    ) %>%
    summarise(
      edf_median = median(edf),
      edf_2_5_percent = quantile(edf, 0.025, names = FALSE),
      edf_97_5_percent = quantile(edf, 0.975, names = FALSE),
      .groups = "drop"
    )

  # Create one readable 12-row table -------------------------------------------

  lcbd_curve_summary <- curve_summary %>%
    left_join(
      shape_summary,
      by = c(
        "focal_gradient",
        "focal_label",
        "focal_variable",
        "lcbd_component"
      )
    ) %>%
    left_join(
      edf_summary,
      by = c(
        "focal_gradient",
        "focal_label",
        "focal_variable",
        "lcbd_component"
      )
    ) %>%
    transmute(
      focal_gradient = focal_label,
      lcbd_component = recode(
        lcbd_component,
        total = "Total",
        replacement = "Replacement",
        richness_difference = "Richness difference"
      ),
      n_allocations,
      gradient_low,
      gradient_high,
      dominant_curve_shape,
      shape_agreement,
      edf_median,
      edf_2_5_percent,
      edf_97_5_percent,
      endpoint_change_median,
      endpoint_change_2_5_percent,
      endpoint_change_97_5_percent,
      curve_effect_range_median,
      curve_effect_range_2_5_percent,
      curve_effect_range_97_5_percent,
      median_gradient_at_curve_minimum,
      median_gradient_at_curve_maximum,
      units = "Pool-size-standardised LCBD partial-effect units",
      interval_definition =
        "2.5th and 97.5th percentiles among geographical allocations"
    ) %>%
    arrange(focal_gradient, lcbd_component)

  write_csv(lcbd_curve_summary, curve_summary_path)

  print(lcbd_curve_summary)
  message("Finished LCBD publication tables")
  message("Wrote: ", main_result_path)
  message("Wrote: ", allocation_validation_path)
  message("Wrote: ", curve_summary_path)



  # Produce LCBD curve figure ---------------------------------------------------
  # Paths -----------------------------------------------------------------------

  figure_dir <- file.path(cfg$figure_dir, "lcbd_plots")
  svg_path <- file.path(figure_dir, "04_local_uniqueness.svg")

  # Validate output destinations before creating folders or exporting.
  invisible(lapply(c(figure_dir), assert_not_raw_output_path))
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  # Reuse the checkpoint loaded for the tables and check its plotting data.
  required_objects <- c("curves_by_allocation", "main_curves")
  missing_objects <- setdiff(required_objects, names(lcbd_results))

  if (length(missing_objects) > 0L) {
    stop(
      "The stage-04 result object lacks: ",
      paste(missing_objects, collapse = ", "),
      call. = FALSE
    )
  }

  allocation_curves <- lcbd_results$curves_by_allocation
  mean_curves <- lcbd_results$main_curves

  focal_label_levels <- c(
    "Agricultural cover",
    "Mean woodland patch area",
    "Woodland cohesion",
    "Noisy Miner detection rate"
  )

  component_levels <- c(
    "Total",
    "Replacement",
    "Richness difference"
  )

  # Prepare the allocation-specific and mean curves ----------------------------

  figure_allocation_data <- allocation_curves %>%
    mutate(
      focal_gradient = factor(
        focal_label,
        levels = focal_label_levels
      ),
      lcbd_component = factor(
        lcbd_component,
        levels = c("total", "replacement", "richness_difference"),
        labels = component_levels
      ),
      allocation = factor(allocation)
    )

  figure_mean_data <- mean_curves %>%
    mutate(
      focal_gradient = factor(
        focal_label,
        levels = focal_label_levels
      ),
      lcbd_component = factor(
        lcbd_component,
        levels = c("total", "replacement", "richness_difference"),
        labels = component_levels
      )
    )

  # Reproduce the accepted LCBD curve figure -----------------------------------

  # Each faint line is one geographical allocation, not a bootstrap resample of
  # sites. The black line is the mean partial effect across the 100 allocations.
  local_uniqueness_figure <- ggplot() +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "grey60",
      linewidth = 0.5
    ) +
    geom_line(
      data = figure_allocation_data,
      aes(
        x = gradient_value,
        y = partial_effect,
        group = allocation
      ),
      colour = "#94346E",
      alpha = 0.09,
      linewidth = 0.2
    ) +
    geom_line(
      data = figure_mean_data,
      aes(
        x = gradient_value,
        y = mean_partial_effect
      ),
      colour = "black",
      linewidth = 0.6
    ) +
    facet_grid2(
      lcbd_component ~ focal_gradient,
      scales = "free",
      independent = "y"
    ) +
    labs(
      x = NULL,
      y = "Partial effect on pool-size-standardised LCBD"
    ) +
    theme_bw(base_size = 12) +
    theme(
      strip.background = element_rect(
        fill = "grey98",
        colour = "black"
      ),
      strip.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey95"),
      axis.text = element_text(colour = "black"),
      legend.position = "none"
    )



  ggsave(
    filename = svg_path,
    plot = local_uniqueness_figure,
    device = svg,
    width = 11,
    height = 7,
    units = "in",
    bg = "white"
  )

  message("Finished LCBD publication figure")
  message("Wrote: ", svg_path)


  message("Completed LCBD publication outputs")

  output_files <- publication_output_paths(cfg)$lcbd
  missing_files <- output_files[!file.exists(output_files)]
  if (length(missing_files) > 0L) {
    stop("Publication files were not created: ",
         paste(missing_files, collapse = ", "), call. = FALSE)
  }
  output_files
}
