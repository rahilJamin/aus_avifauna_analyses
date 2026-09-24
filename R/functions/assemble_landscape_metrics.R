#' Assemble the landscape predictors from three extraction jobs
#'
#' Creates the wide tables, fills absent PLAND classes with zero where coverage
#' exists, derives habitat measures, selects ten predictors, and applies the
#' complete-case filter. A hex_id is the
#' grid polygon identifier used to align rows across rasters and bird matrices.
#'
#' @param vegetation_metrics Long sample_lsm output for vegetation classes.
#' @param woodland_metrics Long sample_lsm output for binary woodland.
#' @param landuse_metrics Long sample_lsm output for land-use classes.
#' @param settings Extraction settings from cfg$landscape$extraction.
#' @return List with all_metrics, selected_metrics (the main input schema), and
#'   a small attrition table. Raw coverage statistics remain in the three input
#'   targets and are not converted into additional filters here.
assemble_landscape_metrics <- function(vegetation_metrics, woodland_metrics,
                                       landuse_metrics, settings) {
  suffix <- paste0(settings$buffer_m, "m")

  # Each source yields at most one value for a site/class/metric combination.
  # Duplicate keys indicate an extraction problem; do not average them away.
  reshape_job <- function(raw, labels, metric_functions) {
    required <- c("plot_id", "class", "metric", "value")
    if (!is.data.frame(raw) || !all(required %in% names(raw)) ||
        !is.numeric(raw$value)) {
      stop("Extraction tables require plot_id, class, metric and numeric value.",
           call. = FALSE)
    }
    if (anyNA(raw$plot_id) || any(!nzchar(as.character(raw$plot_id)))) {
      stop("Extraction table contains missing site identifiers.", call. = FALSE)
    }
    if (anyDuplicated(raw[c("plot_id", "class", "metric")])) {
      stop("Duplicate site/class/metric values in landscape extraction.", call. = FALSE)
    }
    metric_names <- sub("^lsm_c_", "", metric_functions)
    if (anyNA(raw$metric) || any(!raw$metric %in% metric_names)) {
      stop("Unexpected metric in landscape extraction.", call. = FALSE)
    }
    raw <- raw[as.character(raw$class) %in% names(labels), , drop = FALSE]
    wide <- data.frame(hex_id = unique(raw$plot_id))

    # Declare columns even if a class is absent from the whole subset. This keeps
    # the output schema stable for small or geographically restricted datasets.
    for (metric in metric_names) {
      for (class_code in names(labels)) {
        column <- paste(metric, labels[[class_code]], suffix, sep = "_")
        entries <- raw[raw$metric == metric &
                         as.character(raw$class) == class_code, , drop = FALSE]
        wide[[column]] <- entries$value[match(wide$hex_id, entries$plot_id)]
      }
    }
    wide
  }

  vegetation_labels <- settings$vegetation_labels[
    !names(settings$vegetation_labels) %in% as.character(settings$vegetation_drop_classes)
  ]
  wide <- reshape_job(vegetation_metrics, vegetation_labels,
                      settings$vegetation_metrics)
  woodland <- reshape_job(
    woodland_metrics, stats::setNames("woodland", as.character(settings$woodland_class)),
    settings$woodland_metrics
  )
  landuse <- reshape_job(landuse_metrics, settings$landuse_labels,
                         settings$landuse_metrics)

  # Use the vegetation site population, then align binary woodland and land-use
  # values by hex_id rather than relying on row position.
  for (right in list(woodland, landuse)) {
    for (column in setdiff(names(right), "hex_id")) {
      if (column %in% names(wide)) {
        stop("Landscape column appears in more than one job: ", column, call. = FALSE)
      }
      wide[[column]] <- right[[column]][match(wide$hex_id, right$hex_id)]
    }
  }

  # Any observed PLAND value across
  # vegetation OR land use allows missing PLAND cells in that row to become 0.
  # Entirely missing rows stay missing. Patch area/cohesion NAs never become 0.
  pland_columns <- grep("^pland_", names(wide), value = TRUE)
  has_any_data <- rowSums(!is.na(wide[pland_columns])) > 0
  for (column in pland_columns) {
    wide[[column]][has_any_data & is.na(wide[[column]])] <- 0
  }

  metric_column <- function(name) wide[[paste(name, suffix, sep = "_")]]
  habitat <- metric_column("pland_wet_dense") + metric_column("pland_sclerophyll") +
    metric_column("pland_mallee_arid")
  native <- metric_column("pland_native_prod")
  agriculture <- metric_column("pland_agriculture")
  wide[[paste0("PLAND_habitat_", suffix)]] <- habitat
  wide[[paste0("prop_wet_of_habitat_", suffix)]] <- ifelse(
    habitat > 0, metric_column("pland_wet_dense") / habitat, NA_real_
  )
  wide[[paste0("matrix_permeability_", suffix)]] <- ifelse(
    native + agriculture > 0, native / (native + agriculture), NA_real_
  )

  # Dropping incomplete rows across all ten selected predictors is intentional,
  # even though the downstream final models use fewer predictors.
  selected_columns <- c("hex_id", paste(settings$selected_metrics, suffix, sep = "_"))
  missing_columns <- setdiff(selected_columns, names(wide))
  if (length(missing_columns) > 0L) {
    stop("Missing selected landscape columns: ",
         paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  selected <- wide[selected_columns]
  selected <- selected[stats::complete.cases(selected), , drop = FALSE]
  rownames(selected) <- NULL

  list(
    all_metrics = wide,
    selected_metrics = selected,
    attrition = data.frame(
      step = c("sites_in_vegetation_table", "sites_with_complete_selected_metrics"),
      n_sites = c(nrow(wide), nrow(selected))
    )
  )
}
