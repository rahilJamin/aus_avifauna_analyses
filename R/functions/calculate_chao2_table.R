import::from(magrittr, `%>%`)
import::from(dplyr, bind_rows, mutate)
import::from(utils, setTxtProgressBar, txtProgressBar)

#' Estimate Chao2 richness and completeness for all hexagons
#'
#' Splits a hex-level occurrence table by 2.5 km grid cell (`hex_id`) and applies
#' `estimate_chao2_for_hex()` to each occupied grid cell.
#'
#' @param hex_occ Hex-level occurrence table.
#' @param incidence_unit Incidence-unit definition passed to
#'   `estimate_chao2_for_hex()`.
#' @param show_progress Logical; if `TRUE`, print a progress bar.
#'
#' @return A tibble with one row per occupied hexagon and Chao2/completeness
#'   statistics.

calculate_chao2_table <- function(hex_occ, incidence_unit, show_progress = TRUE) {
  if (!is.logical(show_progress) || length(show_progress) != 1L ||
      is.na(show_progress)) {
    stop("`show_progress` must be TRUE or FALSE.", call. = FALSE)
  }
  required_columns <- c("hex_id", "species", "date")
  missing_columns <- setdiff(required_columns, names(hex_occ))
  if (length(missing_columns) > 0L) {
    stop(
      "The hexagon occurrence table lacks: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(hex_occ) == 0L || anyNA(hex_occ$hex_id) ||
      any(!nzchar(as.character(hex_occ$hex_id)))) {
    stop("Chao2 calculation requires at least one valid hexagon ID.", call. = FALSE)
  }

  hex_occ$hex_id <- as.character(hex_occ$hex_id)
  hex_ids <- sort(unique(hex_occ$hex_id))
  hex_split <- split(hex_occ, hex_occ$hex_id)

  pb <- NULL
  if (isTRUE(show_progress) && length(hex_ids) > 0) {
    message("Estimating Chao2/iNEXT completeness for ", length(hex_ids), " hexagons")
    pb <- txtProgressBar(min = 0, max = length(hex_ids), style = 3)
    on.exit(close(pb), add = TRUE)
  }

  out <- vector("list", length(hex_ids))
  for (i in seq_along(hex_ids)) {
    this_hex <- as.character(hex_ids[[i]])
    out[[i]] <- estimate_chao2_for_hex(hex_split[[this_hex]], incidence_unit) %>%
      mutate(hex_id = hex_ids[[i]], .before = 1)

    if (!is.null(pb)) setTxtProgressBar(pb, i)
  }

  bind_rows(out)
}
