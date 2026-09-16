import::from(magrittr, `%>%`)
import::from(dplyr, arrange, count, filter, mutate, row_number, select)
import::from(sf, st_as_sf, st_crs, st_drop_geometry, st_intersects, st_join,
             st_read, st_transform)

#' Assign cleaned occurrence records to hexagons
#'
#' Converts cleaned longitude/latitude records to points, projects them to the
#' hexagon grid CRS, joins them to intersecting hexagons, and keeps one hexagon
#' assignment per input record.
#'
#' @param cleaned_records Data frame with `longitude` and `latitude` columns.
#' @param hex_path Path to the hexagon grid spatial layer.
#' @param input_lonlat_crs Coordinate reference system of the occurrence
#'   longitude and latitude values.
#'
#' @return A data frame of cleaned records with joined hexagon attributes,
#'   including `hex_id`, the unique 2.5 km grid-cell identifier used as the site
#'   ID in downstream matrices.

assign_records_to_hex <- function(cleaned_records, hex_path, input_lonlat_crs) {
  required_columns <- c("longitude", "latitude")
  missing_columns <- setdiff(required_columns, names(cleaned_records))
  if (length(missing_columns) > 0L) {
    stop(
      "The cleaned occurrence table lacks: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  coordinate_columns_are_numeric <- vapply(
    cleaned_records[required_columns],
    is.numeric,
    logical(1)
  )
  if (!all(coordinate_columns_are_numeric)) {
    stop("Occurrence longitude and latitude must be numeric.", call. = FALSE)
  }
  if (anyNA(cleaned_records[required_columns]) ||
      any(!is.finite(as.matrix(cleaned_records[required_columns])))) {
    stop("Occurrence coordinates must be present and finite.", call. = FALSE)
  }

  input_crs <- st_crs(input_lonlat_crs)
  if (is.na(input_crs)) {
    stop("`input_lonlat_crs` is not a valid coordinate reference system.", call. = FALSE)
  }

  hex_grid <- st_read(hex_path, quiet = TRUE) %>%
    mutate(hex_id = as.character(hex_id))
  assert_unique_ids(hex_grid, "hex_id", "hexagon grid")
  if (is.na(st_crs(hex_grid))) {
    stop("The hexagon grid has no coordinate reference system.", call. = FALSE)
  }

  records_sf <- st_as_sf(
    cleaned_records %>% mutate(record_id = row_number()),
    coords = c("longitude", "latitude"),
    crs = input_crs,
    remove = FALSE
  ) %>%
    st_transform(st_crs(hex_grid))

  joined <- st_join(records_sf, hex_grid, join = st_intersects, left = FALSE)

  ambiguous_records <- joined %>%
    st_drop_geometry() %>%
    count(record_id, name = "n_hexagons") %>%
    filter(n_hexagons > 1L)
  if (nrow(ambiguous_records) > 0L) {
    stop(
      nrow(ambiguous_records),
      " occurrence record(s) intersect more than one hexagon. Resolve boundary assignments before continuing.",
      call. = FALSE
    )
  }

  joined %>%
    st_drop_geometry() %>%
    arrange(record_id) %>%
    select(-record_id)
}
