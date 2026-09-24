#' Prepare projected hexagon centroids for landscape sampling
#'
#' A hex_id identifies one grid polygon. Retains the requested polygons in
#' grid-file order, transforms them to the landscape raster CRS, and calculates
#' the centroid used as the centre of each sampling window.
#'
#' @param hex_path Hexagon polygon file containing unique hex_id values.
#' @param site_ids Identifiers of the hexagons to sample.
#' @param raster_crs Projected CRS of the raster used for sampling.
#' @return An sf point table containing hex_id and projected centroid geometry.
prepare_landscape_sites <- function(hex_path, site_ids, raster_crs) {
  requested <- data.frame(hex_id = as.character(site_ids))
  assert_unique_ids(requested, label = "requested landscape sites")
  if (nrow(requested) == 0L) {
    stop("At least one landscape site is required.", call. = FALSE)
  }
  grid <- sf::st_read(hex_path, quiet = TRUE)
  assert_unique_ids(grid, label = "hexagon grid")
  missing_ids <- setdiff(requested$hex_id, as.character(grid$hex_id))
  if (length(missing_ids) > 0L) {
    stop("Requested landscape sites are absent from the grid: ",
         paste(utils::head(missing_ids, 5L), collapse = ", "), call. = FALSE)
  }
  grid <- grid[as.character(grid$hex_id) %in% requested$hex_id, "hex_id"]
  grid <- sf::st_transform(grid, raster_crs)
  if (isTRUE(sf::st_is_longlat(grid))) {
    stop("Landscape sampling requires a projected CRS in metres.", call. = FALSE)
  }
  sf::st_geometry(grid) <- sf::st_centroid(sf::st_geometry(grid))
  grid
}
