#' Sample class-level landscape metrics around hexagon centroids
#'
#' Runs one of three extraction jobs: vegetation classes, binary woodland, or
#' land-use classes. Returns the full sample_lsm table, including its coverage
#' information. targets caches each job separately.
#'
#' @param raster_path Reclassified GeoTIFF file.
#' @param sites Projected sf point table with unique hex_id values.
#' @param buffer_m Circular sampling radius in metres.
#' @param metrics Character vector of landscapemetrics function names.
#' @param show_progress Whether sample_lsm should display progress.
#' @return Long-format metric table with plot_id, class, metric, and value columns.
extract_landscape_metrics <- function(raster_path, sites, buffer_m, metrics,
                                      show_progress) {
  assert_unique_ids(sites, label = "landscape sampling sites")
  if (!is.numeric(buffer_m) || length(buffer_m) != 1L ||
      !is.finite(buffer_m) || buffer_m <= 0) {
    stop("`buffer_m` must be one positive distance in metres.", call. = FALSE)
  }
  raster <- terra::rast(raster_path)
  if (!inherits(sites, "sf") || nrow(sites) == 0L ||
      !all(sf::st_geometry_type(sites) == "POINT") ||
      isTRUE(sf::st_is_longlat(sites)) ||
      !isTRUE(sf::st_crs(sites) == sf::st_crs(terra::crs(raster)))) {
    stop("Sampling sites must be points in the same projected CRS as the raster.",
         call. = FALSE)
  }
  message("Extracting ", paste(metrics, collapse = ", "),
          " for ", nrow(sites), " hexagons")
  started <- proc.time()[["elapsed"]]

  # Use sample_lsm defaults for patch connectivity and raster-edge handling.
  result <- landscapemetrics::sample_lsm(
    landscape = raster, y = sites, plot_id = sites$hex_id,
    shape = "circle", size = buffer_m, what = metrics, progress = show_progress
  )
  message("Extraction finished in ",
          round((proc.time()[["elapsed"]] - started) / 60, 1), " minutes")
  result
}
