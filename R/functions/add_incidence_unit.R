import::from(magrittr, `%>%`)
import::from(dplyr, mutate)

#' Add an incidence unit column
#'
#' Creates the incidence-unit identifier used by iNEXT/Chao2 calculations and
#' Noisy Miner detection. An incidence unit is the repeat-survey unit used to
#' count how often a species was detected at a site. The main workflow uses
#' month-year incidence units.
#'
#' @param data Occurrence table containing at least `date`; the
#'   `checklist_event` option also uses `hex_id`, `dataResourceName`, and
#'   `samplingProtocol`.
#' @param incidence_unit Incidence-unit definition. Supported values are
#'   `"month_year"` and `"checklist_event"`.
#'
#' @return The input data with an `incidence_unit` column.

add_incidence_unit <- function(data, incidence_unit) {
  if (!is.character(incidence_unit) || length(incidence_unit) != 1L ||
      is.na(incidence_unit) || !nzchar(incidence_unit)) {
    stop("`incidence_unit` must be one non-empty character value.", call. = FALSE)
  }

  if (!"date" %in% names(data)) {
    stop("The occurrence table lacks the required `date` column.", call. = FALSE)
  }

  parsed_dates <- as.Date(data$date)
  if (anyNA(parsed_dates)) {
    stop("Incidence-unit dates must be complete and valid.", call. = FALSE)
  }

  if (incidence_unit == "month_year") {
    data %>% mutate(incidence_unit = format(parsed_dates, "%Y-%m"))
  } else if (incidence_unit == "checklist_event") {
    required_columns <- c("hex_id", "dataResourceName", "samplingProtocol")
    missing_columns <- setdiff(required_columns, names(data))
    if (length(missing_columns) > 0L) {
      stop(
        "Checklist-event incidence requires: ",
        paste(missing_columns, collapse = ", "),
        call. = FALSE
      )
    }

    data %>%
      mutate(incidence_unit = paste(hex_id, parsed_dates, dataResourceName,
                                    samplingProtocol, sep = "__"))
  } else {
    stop("Unknown incidence_unit: ", incidence_unit, call. = FALSE)
  }
}
