#' Produce regional convergence tables and publication figures
#'
#' Reads the completed model checkpoint and reproduces the validated
#' publication tables and plots. Plot data, scales, colours, and export
#' settings are carried over from the original output script unchanged.
#'
#' @param cfg Main workflow configuration with input and output paths.
#' @return Paths to the generated files, for a `targets` file target.

import::from(dplyr, arrange, desc, filter, left_join, mutate, rename, select)
import::from(ggplot2, aes, element_blank, element_text, expansion, facet_wrap,
             geom_errorbar, geom_point, ggplot, ggsave, labs, margin,
             position_jitter, scale_colour_manual, scale_fill_manual,
             scale_x_discrete, scale_y_continuous, theme, theme_classic, vars)
import::from(grDevices, svg)
import::from(grid, unit)
import::from(magrittr, `%>%`)
import::from(readr, write_csv)
import::from(tidyr, pivot_wider)

produce_regional_outputs <- function(cfg) {
  # Produce regional tables -----------------------------------------------------
  # Paths ----------------------------------------------------------------------

  input_path <- cfg$regional_results_path

  table_dir <- cfg$table_dir

  support_table_path <- file.path(table_dir, "05_regional_support.csv")
  regional_result_path <- file.path(table_dir, "05_regional_results.csv")
  pair_result_path <- file.path(table_dir, "05_subregion_pair_means.csv")
  subregion_result_path <- file.path(table_dir, "05_subregion_results.csv")
  figure_statistics_path <- file.path(table_dir, "05_figure_statistics.csv")

  if (!file.exists(input_path)) {
    stop("Missing stage-05 result object: ", input_path, call. = FALSE)
  }

  # Validate output destinations before creating folders or exporting.
  invisible(lapply(c(table_dir), assert_not_raw_output_path))
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  # Load and check the saved analysis ------------------------------------------

  regional_objects <- read_rds_checked(input_path, "saved analysis result", required_class = "list")

  required_objects <- c(
    "regional_support",
    "subregion_pair_means",
    "subregion_plot_points",
    "subregion_contrasts",
    "endpoint_summary",
    "endpoint_contrasts"
  )

  missing_objects <- setdiff(required_objects, names(regional_objects))
  if (length(missing_objects) > 0L) {
    stop(
      "The stage-05 result object lacks: ",
      paste(missing_objects, collapse = ", "),
      call. = FALSE
    )
  }

  regional_support <- regional_objects$regional_support
  subregion_pair_means <- regional_objects$subregion_pair_means
  subregion_plot_points <- regional_objects$subregion_plot_points
  subregion_contrasts <- regional_objects$subregion_contrasts
  endpoint_summary <- regional_objects$endpoint_summary
  endpoint_contrasts <- regional_objects$endpoint_contrasts

  # Export endpoint support ----------------------------------------------------

  regional_support_table <- regional_support %>%
    select(
      gradient,
      gradient_label,
      predictor,
      ibra_subregion,
      retained,
      n_sites_all,
      n_sites_low,
      n_sites_high,
      low_threshold,
      high_threshold,
      threshold_rule,
      minimum_predictor_value,
      maximum_predictor_value
    ) %>%
    arrange(gradient, desc(retained), ibra_subregion)

  write_csv(regional_support_table, support_table_path)

  # Export the directly observed means and SDs for each subregion pair ---------

  write_csv(subregion_pair_means, pair_result_path)

  # Export the three components for each of the four gradients ----------------

  regional_result_table <- endpoint_contrasts %>%
    select(
      gradient,
      gradient_label,
      low_threshold,
      high_threshold,
      component,
      low_end_mean,
      high_end_mean,
      high_minus_low,
      difference_bootstrap_se,
      difference_lower_95,
      difference_upper_95,
      bootstrap_two_sided_p_value,
      interpretation
    )

  write_csv(regional_result_table, regional_result_path)

  # Export one low-high comparison for every retained subregion ---------------

  subregion_result_table <- subregion_contrasts %>%
    left_join(
      subregion_plot_points %>%
        select(
          gradient,
          component,
          ibra_subregion,
          endpoint,
          sd_among_subregion_pairs
        ) %>%
        pivot_wider(
          names_from = endpoint,
          values_from = sd_among_subregion_pairs,
          names_prefix = "sd_"
        ),
      by = c("gradient", "component", "ibra_subregion")
    ) %>%
    rename(
      low_end_sd_among_subregion_pairs = `sd_Low end`,
      high_end_sd_among_subregion_pairs = `sd_High end`
    ) %>%
    arrange(gradient, component, ibra_subregion)

  write_csv(subregion_result_table, subregion_result_path)

  # Export compact values for manual figure annotation ------------------------

  figure_statistics <- endpoint_summary %>%
    left_join(
      endpoint_contrasts %>%
        select(
          gradient,
          component,
          high_minus_low,
          difference_bootstrap_se,
          difference_lower_95,
          difference_upper_95,
          bootstrap_two_sided_p_value
        ),
      by = c("gradient", "component")
    ) %>%
    select(
      gradient,
      gradient_label,
      component,
      endpoint,
      endpoint_definition,
      low_threshold,
      high_threshold,
      n_subregions,
      n_subregion_pairs,
      endpoint_mean = estimate,
      sd_among_subregion_pair_means,
      endpoint_bootstrap_se = bootstrap_se,
      endpoint_lower_95 = lower_95,
      endpoint_upper_95 = upper_95,
      high_minus_low,
      difference_bootstrap_se,
      difference_lower_95,
      difference_upper_95,
      bootstrap_two_sided_p_value
    )

  write_csv(figure_statistics, figure_statistics_path)

  print(regional_result_table, n = Inf)

  message("Finished regional convergence tables")



  # Produce regional figures ----------------------------------------------------
  # Figure settings ------------------------------------------------------------

  # Change these values without rerunning the analysis.
  component_to_plot <- "Total Jaccard"

  # Cool petrol teal and warm burnt coral distinguish the two endpoints.
  endpoint_colours <- c(
    "Low end" = "#2F6F73",
    "High end" = "#D06B4F"
  )

  background_fill <- "#BDBDBD"
  background_outline <- "#6F6F6F"
  background_point_size <- 2.0
  background_point_alpha <- 0.72
  summary_point_size <- 3.2
  base_text_size <- 10
  figure_width <- 3.6
  figure_height <- 3.6
  combined_figure_width <- 7.2
  combined_figure_height <- 6.2

  # Paths ----------------------------------------------------------------------

  figure_dir <- file.path(cfg$figure_dir, "regional_plots")
  combined_figure_path <- file.path(
    figure_dir,
    "05_regional_convergence_2x2.svg"
  )

  # Validate output destinations before creating folders or exporting.
  invisible(lapply(c(figure_dir), assert_not_raw_output_path))
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  # Reuse the checkpoint loaded for the tables and check its plotting data.
  required_objects <- c("subregion_plot_points", "endpoint_summary", "settings")
  missing_objects <- setdiff(required_objects, names(regional_objects))

  if (length(missing_objects) > 0L) {
    stop(
      "The stage-05 result object lacks: ",
      paste(missing_objects, collapse = ", "),
      call. = FALSE
    )
  }

  gradient_definitions <- regional_objects$settings$gradient_definitions

  if (is.null(gradient_definitions)) {
    stop("The stage-05 result object lacks gradient definitions.", call. = FALSE)
  }

  subregion_points <- regional_objects$subregion_plot_points %>%
    filter(component == component_to_plot)

  endpoint_summary <- regional_objects$endpoint_summary %>%
    filter(component == component_to_plot)

  endpoint_levels <- c("Low end", "High end")
  endpoint_axis_labels <- c(
    "Low end" = "Low endpoint",
    "High end" = "High endpoint"
  )

  # Use one data-derived y-axis across all four SVGs. This keeps both the canvas
  # and the visual scale identical when the panels are assembled in Inkscape.
  all_y_values <- c(
    subregion_points$mean_dissimilarity_to_other_subregions,
    endpoint_summary$lower_95,
    endpoint_summary$upper_95
  )

  shared_y_limits <- c(
    floor(min(all_y_values, na.rm = TRUE) * 20) / 20,
    ceiling(max(all_y_values, na.rm = TRUE) * 20) / 20
  )

  shared_y_breaks <- seq(
    shared_y_limits[1],
    shared_y_limits[2],
    by = 0.05
  )

  # Export one separate SVG per gradient --------------------------------------

  for (gradient_index in seq_len(nrow(gradient_definitions))) {
    current_definition <- gradient_definitions[gradient_index, ]
    current_gradient <- current_definition$gradient

    current_subregion_points <- subregion_points %>%
      filter(gradient == current_gradient) %>%
      mutate(
        endpoint = factor(endpoint, levels = endpoint_levels)
      )

    current_endpoint_summary <- endpoint_summary %>%
      filter(gradient == current_gradient) %>%
      mutate(endpoint = factor(endpoint, levels = endpoint_levels))

    if (
      nrow(current_subregion_points) == 0L ||
        nrow(current_endpoint_summary) != 2L
    ) {
      stop(
        "Incomplete plotting data for gradient: ",
        current_gradient,
        call. = FALSE
      )
    }

    regional_convergence_figure <- ggplot() +
      geom_point(
        data = current_subregion_points,
        aes(
          x = endpoint,
          y = mean_dissimilarity_to_other_subregions
        ),
        position = position_jitter(
          width = 0.055,
          height = 0,
          seed = 4202
        ),
        shape = 21,
        fill = background_fill,
        colour = background_outline,
        stroke = 0.35,
        size = background_point_size,
        alpha = background_point_alpha
      ) +
      geom_errorbar(
        data = current_endpoint_summary,
        aes(
          x = endpoint,
          ymin = lower_95,
          ymax = upper_95,
          colour = endpoint
        ),
        width = 0.08,
        linewidth = 0.75
      ) +
      geom_point(
        data = current_endpoint_summary,
        aes(x = endpoint, y = estimate, fill = endpoint),
        shape = 21,
        colour = "black",
        stroke = 0.45,
        size = summary_point_size
      ) +
      scale_x_discrete(labels = endpoint_axis_labels) +
      scale_colour_manual(values = endpoint_colours, guide = "none") +
      scale_fill_manual(values = endpoint_colours, guide = "none") +
      scale_y_continuous(
        limits = shared_y_limits,
        breaks = shared_y_breaks,
        expand = expansion(mult = c(0.02, 0.02))
      ) +
      labs(
        x = NULL,
        y = "Mean total β-diversity\nbetween subregions"
      ) +
      theme_classic(base_size = base_text_size) +
      theme(
        axis.text = element_text(colour = "black"),
        axis.title.y = element_text(margin = margin(r = 7)),
        legend.position = "none",
        plot.margin = margin(5.5, 5.5, 5.5, 5.5, unit = "pt")
      )

    figure_path <- file.path(
      figure_dir,
      paste0("05_regional_convergence_", current_gradient, ".svg")
    )

    ggsave(
      filename = figure_path,
      plot = regional_convergence_figure,
      device = svg,
      width = figure_width,
      height = figure_height,
      units = "in",
      bg = "white"
    )

    message("Exported: ", figure_path)
  }

  # Export the publication-ready 2 x 2 layout ---------------------------------

  panel_labels <- c(
    agriculture = "a) Agricultural cover",
    patch_area = "b) Mean woodland patch area",
    cohesion = "c) Woodland cohesion",
    noisy_miner = "d) Noisy Miner detection rate"
  )

  combined_subregion_points <- subregion_points %>%
    mutate(
      endpoint = factor(endpoint, levels = endpoint_levels),
      panel = factor(
        panel_labels[gradient],
        levels = unname(panel_labels)
      )
    )

  combined_endpoint_summary <- endpoint_summary %>%
    mutate(
      endpoint = factor(endpoint, levels = endpoint_levels),
      panel = factor(
        panel_labels[gradient],
        levels = unname(panel_labels)
      )
    )

  regional_convergence_2x2 <- ggplot() +
    geom_point(
      data = combined_subregion_points,
      aes(
        x = endpoint,
        y = mean_dissimilarity_to_other_subregions
      ),
      position = position_jitter(
        width = 0.055,
        height = 0,
        seed = 4202
      ),
      shape = 21,
      fill = background_fill,
      colour = background_outline,
      stroke = 0.35,
      size = background_point_size,
      alpha = background_point_alpha
    ) +
    geom_errorbar(
      data = combined_endpoint_summary,
      aes(
        x = endpoint,
        ymin = lower_95,
        ymax = upper_95,
        colour = endpoint
      ),
      width = 0.08,
      linewidth = 0.75
    ) +
    geom_point(
      data = combined_endpoint_summary,
      aes(x = endpoint, y = estimate, fill = endpoint),
      shape = 21,
      colour = "black",
      stroke = 0.45,
      size = summary_point_size
    ) +
    facet_wrap(
      vars(panel),
      ncol = 2,
      scales = "free_x"
    ) +
    scale_x_discrete(labels = endpoint_axis_labels) +
    scale_colour_manual(values = endpoint_colours, guide = "none") +
    scale_fill_manual(values = endpoint_colours, guide = "none") +
    scale_y_continuous(
      limits = shared_y_limits,
      breaks = shared_y_breaks,
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    labs(
      x = NULL,
      y = "Mean total β-diversity between subregions"
    ) +
    theme_classic(base_size = base_text_size) +
    theme(
      axis.text = element_text(colour = "black", size = 9.5),
      axis.title.y = element_text(
        size = 10.5,
        margin = margin(r = 9)
      ),
      strip.background = element_blank(),
      strip.text = element_text(
        colour = "black",
        face = "plain",
        size = 10.5,
        hjust = 0
      ),
      panel.spacing = unit(1.2, "lines"),
      legend.position = "none",
      plot.margin = margin(8, 8, 8, 10, unit = "pt")
    )

  ggsave(
    filename = combined_figure_path,
    plot = regional_convergence_2x2,
    device = svg,
    width = combined_figure_width,
    height = combined_figure_height,
    units = "in",
    bg = "white"
  )

  message("Exported: ", combined_figure_path)
  message("Finished regional convergence figures")


  message("Completed regional convergence publication outputs")

  output_files <- publication_output_paths(cfg)$regional
  missing_files <- output_files[!file.exists(output_files)]
  if (length(missing_files) > 0L) {
    stop("Publication files were not created: ",
         paste(missing_files, collapse = ", "), call. = FALSE)
  }
  output_files
}
