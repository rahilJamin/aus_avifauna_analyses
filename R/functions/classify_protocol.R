import::from(dplyr, case_when, coalesce)
import::from(stringr, regex, str_detect)

#' Classify survey protocol type
#'
#' Classifies retained records into broad protocol classes used by the occurrence
#' cleaning stage. The classification is based on protocol text and, for
#' iNaturalist, data-resource name.
#'
#' @param data_resource Character vector of data-resource names.
#' @param sampling_protocol Character vector of sampling-protocol descriptions.
#'
#' @return A character vector with protocol classes:
#'   `"community_complete_or_structured"`, `"opportunistic"`, or
#'   `"semi_structured_or_unclear"`.

classify_protocol <- function(data_resource, sampling_protocol) {
  resource <- coalesce(data_resource, "")
  protocol <- coalesce(sampling_protocol, "")

  case_when(
    str_detect(
      protocol,
      regex(
        paste(c(
          "Stationary", "Birds Australia", "BirdLife Australia",
          "Garden Bird Survey", "transect", "Point spot", "Bird count",
          "area search", "2ha search", "500m area search",
          "5km area search", "fixed route", "defined area"
        ), collapse = "|"),
        ignore_case = TRUE
      )
    ) ~ "community_complete_or_structured",
    str_detect(
      protocol,
      regex("Incidental|Opportunistic|iNaturalist", ignore_case = TRUE)
    ) | resource == "iNaturalist Australia" ~ "opportunistic",
    TRUE ~ "semi_structured_or_unclear"
  )
}
