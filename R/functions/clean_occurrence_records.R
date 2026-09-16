import::from(CoordinateCleaner, clean_coordinates)
import::from(magrittr, `%>%`)
import::from(dplyr, any_of, bind_rows, distinct, filter, mutate, n_distinct,
             select)
import::from(readr, read_csv)
import::from(stringr, regex, str_detect)
import::from(tibble, tibble)
import::from(tidyr, replace_na)

#' Clean raw ALA occurrence records
#'
#' Filters raw occurrence records to presence records in the configured date
#' window, removes records with missing taxon/date/coordinate fields, screens
#' survey protocols, applies coordinate-cleaning checks, deduplicates
#' species-date-coordinate records, and applies broad range filters for selected
#' Tasmanian endemic and mainland obligate species.
#'
#' @param raw_path Path to the raw ALA occurrence CSV.
#' @param start_date First date retained.
#' @param end_date Last date retained.
#' @param coord_precision Maximum accepted coordinate uncertainty in metres.
#' @param cleaning_rules Named list containing trusted resources, protocol
#'   exclusions, Tasmanian endemics, and mainland-obligate species.
#' @param protocol_filter Protocol filter to apply. Use `"all_retained"` for the
#'   main analysis or `"community_complete_or_structured"` for the stricter
#'   structured-protocol subset.
#' @param show_progress Logical; if `TRUE`, print progress messages after each
#'   cleaning checkpoint.
#' @param coordinate_cleaner Function used to add CoordinateCleaner validity
#'   flags. The default is `CoordinateCleaner::clean_coordinates()`; accepting
#'   the function as an argument allows deterministic unit testing.
#'
#' @return A list with `records`, the cleaned occurrence table, and `attrition`,
#'   a table of record/species counts after each filtering step.

clean_occurrence_records <- function(raw_path,
                                     start_date,
                                     end_date,
                                     coord_precision,
                                     cleaning_rules,
                                     protocol_filter = "all_retained",
                                     show_progress = TRUE,
                                     coordinate_cleaner = clean_coordinates) {
  valid_protocol_filters <- c(
    "all_retained",
    "community_complete_or_structured"
  )
  
  if (!protocol_filter %in% valid_protocol_filters) {
    stop(
      "Unknown protocol_filter: ", protocol_filter,
      ". Expected one of: ", paste(valid_protocol_filters, collapse = ", "),
      call. = FALSE
    )
  }

  if (!inherits(start_date, "Date") || length(start_date) != 1L ||
      is.na(start_date) || !inherits(end_date, "Date") ||
      length(end_date) != 1L || is.na(end_date) || start_date > end_date) {
    stop("`start_date` and `end_date` must be valid ordered Date values.", call. = FALSE)
  }

  if (!is.numeric(coord_precision) || length(coord_precision) != 1L ||
      is.na(coord_precision) || !is.finite(coord_precision) ||
      coord_precision < 0) {
    stop("`coord_precision` must be one non-negative finite number.", call. = FALSE)
  }

  required_rule_names <- c(
    "trusted_resources",
    "protocol_exclusions",
    "tasmanian_endemics",
    "mainland_obligates"
  )
  missing_rule_names <- setdiff(required_rule_names, names(cleaning_rules))
  if (length(missing_rule_names) > 0L) {
    stop(
      "The cleaning rule set lacks: ",
      paste(missing_rule_names, collapse = ", "),
      call. = FALSE
    )
  }
  valid_rule_vectors <- vapply(
    cleaning_rules[required_rule_names],
    function(rule) is.character(rule) && !anyNA(rule),
    logical(1)
  )
  if (!all(valid_rule_vectors) ||
      length(cleaning_rules$protocol_exclusions) < 1L ||
      any(!nzchar(cleaning_rules$protocol_exclusions))) {
    stop(
      "Cleaning rules must be character vectors, and protocol exclusions cannot be empty.",
      call. = FALSE
    )
  }
  if (!is.function(coordinate_cleaner)) {
    stop("`coordinate_cleaner` must be a function.", call. = FALSE)
  }
  
  started_at <- Sys.time()
  step_number <- 0L
  total_steps <- 9L
  
  format_elapsed <- function(start_time) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    
    if (elapsed < 60) {
      return(sprintf("%.1f sec", elapsed))
    }
    
    sprintf("%d min %.1f sec", floor(elapsed / 60), elapsed %% 60)
  }
  
  format_count <- function(x) {
    format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
  }
  
  # Retain ALA basis-of-record classes that represent in-situ observations.
  filter_bor <- c(
    "HUMAN_OBSERVATION",
    "MACHINE_OBSERVATION",
    "OBSERVATION",
    "OCCURRENCE"
  )
  
  # Track how many records and species remain after each major cleaning step.
  # The same checkpoint also reports timed progress during manual runs.
  attrition <- list()
  
  add_attrition <- function(data, step, label = step) {
    step_number <<- step_number + 1L
    
    attrition[[length(attrition) + 1]] <<- tibble(
      step = step,
      n_records = nrow(data),
      n_species = n_distinct(data$species, na.rm = TRUE)
    )
    
    if (isTRUE(show_progress)) {
      message(
        sprintf(
          "[%02d/%02d] %s: %s records, %s species; elapsed %s",
          step_number,
          total_steps,
          label,
          format_count(nrow(data)),
          format_count(n_distinct(data$species, na.rm = TRUE)),
          format_elapsed(started_at)
        )
      )
    }
    
    invisible(data)
  }
  
  # Read records and validate the input schema before applying filters. A clear
  # schema error is more useful than a later failure inside a dplyr expression.
  raw <- read_csv(raw_path, show_col_types = FALSE)
  required_columns <- c(
    "occurrenceStatus", "date", "species", "longitude", "latitude",
    "basisOfRecord", "samplingProtocol", "dataResourceName",
    "coordinateUncertaintyInMeters"
  )
  missing_columns <- setdiff(required_columns, names(raw))
  if (length(missing_columns) > 0L) {
    stop(
      "The raw occurrence table lacks: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  raw <- raw %>%
    mutate(date = as.Date(date))

  numeric_coordinate_columns <- c(
    "longitude", "latitude", "coordinateUncertaintyInMeters"
  )
  if (!all(vapply(raw[numeric_coordinate_columns], is.numeric, logical(1)))) {
    stop(
      "Longitude, latitude, and coordinate uncertainty must be numeric.",
      call. = FALSE
    )
  }
  
  add_attrition(raw, "raw_input", "Loaded raw records")
  
  # Keep presence records inside the configured analysis window.
  dat <- raw %>%
    filter(
      occurrenceStatus == "PRESENT",
      date >= start_date,
      date <= end_date
    )
  
  add_attrition(
    dat,
    "present_and_date_window",
    "Filtered presence records and analysis dates"
  )
  
  # Records without species names, dates, or coordinates cannot be assigned to
  # species-level occurrence histories or spatial hexagons.
  dat <- dat %>%
    filter(
      !is.na(species),
      !is.na(date),
      !is.na(longitude),
      !is.na(latitude)
    )
  
  add_attrition(
    dat,
    "required_taxon_date_coordinates",
    "Removed records missing species, date, or coordinates"
  )
  
  # Remove specimen and non-observation records.
  dat <- dat %>%
    filter(basisOfRecord %in% filter_bor)
  
  add_attrition(
    dat,
    "basis_of_record_filter",
    "Filtered basis-of-record classes"
  )
  
  # Remove protocols known to be unsuitable for community occurrence modelling.
  # Missing protocols are retained only for trusted data resources.
  exclusion_regex <- regex(
    paste(cleaning_rules$protocol_exclusions, collapse = "|"),
    ignore_case = TRUE
  )
  
  dat <- dat %>%
    filter(
      !str_detect(replace_na(samplingProtocol, ""), exclusion_regex),
      !is.na(samplingProtocol) |
        (is.na(samplingProtocol) &
           dataResourceName %in% cleaning_rules$trusted_resources)
    ) %>%
    mutate(
      protocol_class = classify_protocol(dataResourceName, samplingProtocol)
    )
  
  if (protocol_filter == "community_complete_or_structured") {
    dat <- dat %>%
      filter(protocol_class == "community_complete_or_structured")
  }
  
  add_attrition(
    dat,
    "protocol_screen_and_classification",
    "Screened and classified survey protocols"
  )
  
  # CoordinateCleaner flags spatially impossible or suspicious coordinates. The
  # original workflow used equal-coordinate, sea, and zero-coordinate checks.
  flags <- coordinate_cleaner(
    x = dat,
    lon = "longitude",
    lat = "latitude",
    species = "species",
    tests = c("equal", "seas", "zeros"),
    value = "spatialvalid"
  )

  required_flag_columns <- c(".sea", ".zer", ".equ")
  if (!is.data.frame(flags) || nrow(flags) != nrow(dat) ||
      !all(required_flag_columns %in% names(flags))) {
    stop(
      "The coordinate cleaner did not return one flagged row per input record.",
      call. = FALSE
    )
  }
  
  add_attrition(
    flags,
    "coordinate_cleaner_flags_added",
    "Added CoordinateCleaner spatial flags"
  )
  
  # Keep records with coordinate uncertainty within the configured precision.
  # When uncertainty is missing, retain records only if the protocol text
  # suggests a standardised, structured, camera, or spatially bounded survey.
  dat <- flags %>%
    mutate(
      is_standardised = str_detect(
        samplingProtocol,
        regex(
          "Birds Australia|BirdLife Australia|Garden Bird Survey",
          ignore_case = TRUE
        )
      ),
      is_camera = str_detect(
        samplingProtocol,
        regex("Camera|camera trap", ignore_case = TRUE)
      ),
      is_structured = str_detect(
        samplingProtocol,
        regex(
          "Stationary|transect|Point spot|Bird count|Area search",
          ignore_case = TRUE
        )
      ),
      has_dimensions = str_detect(
        samplingProtocol,
        regex("\\d+\\s?(m|km|ha)", ignore_case = TRUE)
      )
    ) %>%
    filter(
      coordinateUncertaintyInMeters <= coord_precision |
        (
          is.na(coordinateUncertaintyInMeters) &
            (is_standardised | is_camera | is_structured | has_dimensions)
        ),
      .sea == TRUE,
      .zer == TRUE,
      .equ == TRUE
    ) %>%
    select(-is_standardised, -is_camera, -is_structured, -has_dimensions)
  
  add_attrition(
    dat,
    "coordinate_uncertainty_and_spatial_validity",
    "Filtered coordinate uncertainty and spatial validity"
  )
  
  # Retain one record per species per date at each exact coordinate.
  dat <- dat %>%
    distinct(latitude, longitude, date, species, .keep_all = TRUE)
  
  add_attrition(
    dat,
    "deduplicated_lat_lon_date_species",
    "Removed duplicate species-date-coordinate records"
  )
  
  # Remove records outside accepted broad ranges for Tasmanian endemics and
  # mainland obligates while leaving all other species unchanged.
  dat <- dat %>%
    filter(
      (species %in% cleaning_rules$tasmanian_endemics & latitude <= -39.5) |
        (species %in% cleaning_rules$mainland_obligates & latitude > -39.5) |
        (!species %in% cleaning_rules$tasmanian_endemics &
           !species %in% cleaning_rules$mainland_obligates)
    ) %>%
    select(-any_of(c(".val", ".equ", ".zer", ".sea", ".summary")))
  
  add_attrition(
    dat,
    "tasmanian_endemic_and_mainland_obligate_range_filter",
    "Applied broad species-range filters"
  )
  
  list(records = dat, attrition = bind_rows(attrition))
}
