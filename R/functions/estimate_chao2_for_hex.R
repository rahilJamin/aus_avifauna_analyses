import::from(magrittr, `%>%`)
import::from(dplyr, count, distinct, n_distinct)
import::from(iNEXT, ChaoRichness, DataInfo)
import::from(tibble, tibble)

#' Estimate Chao2 richness and completeness for one hexagon
#'
#' Converts occurrences from one 2.5 km grid cell (`hex_id`) into
#' incidence-frequency input for `iNEXT::ChaoRichness()` and returns observed
#' richness, estimated richness, standard error, completeness, sample coverage,
#' and singleton/doubleton incidence counts. Failed iNEXT calls return a row with
#' `chao_status = "failed"` rather than stopping the whole workflow.
#'
#' @param hex_data Occurrence records for one hexagon.
#' @param incidence_unit Incidence-unit definition passed to
#'   `add_incidence_unit()`.
#'
#' @return A one-row tibble of Chao2 and sampling-completeness statistics.

estimate_chao2_for_hex <- function(hex_data, incidence_unit = "month_year") {
  if (nrow(hex_data) == 0L || !"species" %in% names(hex_data) ||
      anyNA(hex_data$species) ||
      any(!nzchar(as.character(hex_data$species)))) {
    stop("Chao2 estimation requires non-empty species records.", call. = FALSE)
  }

  dat <- add_incidence_unit(hex_data, incidence_unit)

  incidence <- dat %>%
    distinct(species, incidence_unit) %>%
    count(species, name = "freq")

  n_units <- n_distinct(dat$incidence_unit)
  if (n_units < 1L || nrow(incidence) < 1L ||
      any(incidence$freq < 1L) || any(incidence$freq > n_units)) {
    stop("The incidence-frequency vector is internally inconsistent.", call. = FALSE)
  }
  input <- c(n_units, incidence$freq)

  richness <- tryCatch(
    ChaoRichness(input, datatype = "incidence_freq", conf = 0.95),
    error = function(e) NULL
  )
  info <- tryCatch(
    DataInfo(input, datatype = "incidence_freq"),
    error = function(e) NULL
  )

  if (is.null(richness)) {
    return(tibble(
      S_obs = n_distinct(dat$species),
      sampling_units = n_units,
      chao2_estimate = NA_real_,
      chao2_se = NA_real_,
      chao2_completeness = NA_real_,
      sample_coverage = NA_real_,
      Q1 = sum(incidence$freq == 1),
      Q2 = sum(incidence$freq == 2),
      chao_status = "failed"
    ))
  }

  observed <- as.numeric(richness[["Observed"]][1])
  estimate <- as.numeric(richness[["Estimator"]][1])
  se_col <- intersect(c("Est_s.e.", "Est_s.e", "s.e.", "SE"), names(richness))
  se <- if (length(se_col) > 0) {
    as.numeric(richness[[se_col[1]]][1])
  } else {
    NA_real_
  }
  sample_coverage <- if (!is.null(info) && "SC" %in% names(info)) {
    as.numeric(info$SC[1])
  } else {
    NA_real_
  }
  if (!is.finite(observed) || observed < 1 || !is.finite(estimate) ||
      estimate < observed) {
    return(tibble(
      S_obs = n_distinct(dat$species),
      sampling_units = n_units,
      chao2_estimate = NA_real_,
      chao2_se = NA_real_,
      chao2_completeness = NA_real_,
      sample_coverage = NA_real_,
      Q1 = sum(incidence$freq == 1),
      Q2 = sum(incidence$freq == 2),
      chao_status = "failed"
    ))
  }
  if (!is.na(sample_coverage) &&
      (!is.finite(sample_coverage) || sample_coverage < 0 ||
       sample_coverage > 1)) {
    sample_coverage <- NA_real_
  }
  completeness <- if (is.finite(estimate) && estimate > 0) {
    pmin(observed / estimate, 1)
  } else {
    NA_real_
  }

  tibble(
    S_obs = observed,
    sampling_units = n_units,
    chao2_estimate = estimate,
    chao2_se = se,
    chao2_completeness = completeness,
    sample_coverage = sample_coverage,
    Q1 = sum(incidence$freq == 1),
    Q2 = sum(incidence$freq == 2),
    chao_status = "ok"
  )
}
