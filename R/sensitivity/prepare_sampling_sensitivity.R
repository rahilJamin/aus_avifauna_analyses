#' Declare one-change-at-a-time sensitivity scenarios
#' @param config,settings Main and sensitivity settings.
#' @return Named scenario lists for dynamic targets branches.
define_sensitivity_scenarios <- function(config, settings) {
  make <- function(id, label, type, completeness = config$min_completeness,
                   minimum_units = config$min_sampling_units,
                   pool_k = config$lcbd_reference_pool_k) {
    list(id = id, label = label, type = type, completeness = completeness,
         minimum_units = minimum_units, pool_k = pool_k)
  }
  scenarios <- list(
    make("main_70", "Main settings (70% completeness)", "main")
  )
  scenarios <- c(scenarios, lapply(settings$completeness, function(x) make(
    paste0("completeness_", round(100 * x)), paste0(round(100 * x), "% completeness"),
    "completeness", completeness = x)))
  scenarios <- c(scenarios, list(
    make("effort_12", "At least 12 sampled months", "effort", minimum_units = settings$minimum_months),
    make("daily_incidence", "Daily incidence units", "daily"),
    make("structured_protocol", "Structured protocols only", "structured")
  ), lapply(settings$alternative_pool_counts, function(k) make(
    paste0("pools_", k), paste(k, "reference pools"), "pools", pool_k = k)))
  stats::setNames(scenarios, vapply(scenarios, `[[`, "", "id"))
}

#' Read the sampling-sensitivity occurrence subset in bounded chunks
#' @param path Main stage-02 occurrence CSV, read only.
#' @param settings Sensitivity settings.
#' @return Distinct site/species/date/protocol records; no coordinate deduplication.
read_sampling_occurrences <- function(path, settings) {
  callback <- readr::DataFrameCallback$new(function(x, pos) {
    unique(x)
  })
  out <- readr::read_csv_chunked(path, callback,
    chunk_size = settings$occurrence_chunk_size, progress = TRUE,
    col_types = readr::cols_only(hex_id = readr::col_character(),
      species = readr::col_character(), date = readr::col_date(),
      protocol_class = readr::col_character()))
  if (!all(c("hex_id", "species", "date", "protocol_class") %in% names(out)) ||
      !nrow(out) || anyNA(out[c("hex_id", "species", "date")])) {
    stop("Stage-02 records lack valid sampling-sensitivity columns.", call. = FALSE)
  }
  unique(out)
}

#' Label daily sampling units without modifying the main incidence function
#' @param data Occurrence data with complete dates.
#' @param incidence_unit Must be "day".
#' @return Input data with calendar-date incidence labels.
add_daily_incidence_unit <- function(data, incidence_unit) {
  if (!identical(incidence_unit, "day")) stop("This sensitivity labels daily units only.", call. = FALSE)
  dates <- as.Date(data$date)
  if (length(dates) != nrow(data) || anyNA(dates)) stop("Daily incidence requires complete dates.", call. = FALSE)
  data$incidence_unit <- as.character(dates)
  data
}

#' Build daily incidence with the unchanged main Chao2 and detection algorithms
#' @param occurrences Main-footprint occurrences.
#' @param config Main settings.
#' @return The standard matrix bundle, with both effort and detection counted daily.
build_daily_matrix_bundle <- function(occurrences, config) {
  # The main function supports monthly/event units. Give LOCAL COPIES of its
  # four callers the extra daily labeller; main functions, their environments,
  # source files and targets hashes remain untouched. Estimation itself is shared.
  daily <- new.env(parent = environment(build_matrix_bundle))
  daily$add_incidence_unit <- add_daily_incidence_unit
  for (name in c("build_matrix_bundle", "calculate_chao2_table",
                 "estimate_chao2_for_hex", "calculate_miner_detection")) {
    fun <- get(name, envir = environment(build_matrix_bundle), inherits = TRUE)
    environment(fun) <- daily
    assign(name, fun, envir = daily)
  }
  daily$build_matrix_bundle(occurrences, "day", config$min_sampling_units,
    config$min_species, config$min_completeness, config$noisy_miner_species,
    config$show_progress)
}

#' Prepare a full workflow alternative after occurrence assignment
#' @param scenario One declared scenario.
#' @param baseline Validated main inputs.
#' @param occurrences Stage-02 subset, needed only for daily/protocol alternatives.
#' @param config,settings Main and sensitivity settings.
#' @return Prepared LCBD data plus scenario identity, settings and sample counts.
prepare_sensitivity_scenario <- function(scenario, baseline, occurrences, config, settings) {
  if (scenario$type == "completeness") {
    prepared <- prepare_completeness_data(scenario$completeness, baseline, config)
  } else if (scenario$type %in% c("main", "pools")) {
    dat <- baseline$main$analysis_site_data
    y <- baseline$species_matrix[dat$hex_id,
      colnames(baseline$species_matrix) != config$noisy_miner_species, drop = FALSE]
    y <- y[, colSums(y) > 0, drop = FALSE]
    prepared <- list(site_data = dat, species_matrix = y,
      distance_matrices = calculate_lcbd_beta_matrices(y)$dist_mats,
      n_chao_sites = nrow(baseline$species_matrix), attrition = data.frame())
  } else {
    if (scenario$type == "effort") {
      valid <- retain_sampled_hexagons(baseline$chao, scenario$minimum_units,
                                      config$min_species, scenario$completeness)
      ids <- as.character(valid$hex_id)
      matrix <- baseline$species_matrix[ids, , drop = FALSE]
      metadata <- baseline$metadata[match(ids, baseline$metadata$hex_id), , drop = FALSE]
    } else {
      if (scenario$type == "daily") {
        bundle <- build_daily_matrix_bundle(occurrences, config)
      } else {
        selected <- occurrences[!is.na(occurrences$protocol_class) &
          occurrences$protocol_class == settings$structured_protocol, , drop = FALSE]
        bundle <- build_matrix_bundle(selected, config$incidence_unit,
          config$min_sampling_units, config$min_species, config$min_completeness,
          config$noisy_miner_species, config$show_progress)
      }
      matrix <- bundle$species_matrix
      metadata <- bundle$hex_metadata
      ids <- rownames(matrix)
    }
    model <- prepare_lcbd_model_data(matrix, metadata,
      baseline$landscape[baseline$landscape$hex_id %in% ids, , drop = FALSE],
      config$hex_2_5_path, config$noisy_miner_species, config$urban_cover_max, config$model_coord_crs)
    dat <- add_lcbd_pool_coordinates(
      model$candidate_model_df, config$hex_2_5_path, config$lcbd_reference_pool_crs
    )
    prepared <- list(site_data = dat, species_matrix = model$lcbd_species_matrix,
      distance_matrices = model$beta_matrices$dist_mats, n_chao_sites = length(ids), attrition = model$attrition)
  }
  prepared$scenario <- scenario
  prepared$threshold <- scenario$completeness
  prepared
}
