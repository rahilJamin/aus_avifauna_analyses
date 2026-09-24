#' Track a GeoTIFF and its category and georeferencing sidecars
#'
#' CLUM and NVIS may store category tables outside the TIFF, for example in
#' .tif.aux.xml or .tif.vat.dbf. These files affect reclassification and must be
#' tracked with the raster so a changed category table invalidates the target.
#'
#' @param raster_path Path to the source .tif or .tiff file.
#' @return Character vector containing the TIFF and its existing sidecars.
track_raster_files <- function(raster_path) {
  if (!is.character(raster_path) || length(raster_path) != 1L ||
      is.na(raster_path) || !file.exists(raster_path) ||
      !tolower(tools::file_ext(raster_path)) %in% c("tif", "tiff")) {
    stop("Expected an existing GeoTIFF: ", raster_path, call. = FALSE)
  }
  candidates <- list.files(dirname(raster_path), full.names = TRUE)
  raster_name <- tolower(basename(raster_path))
  stem <- tolower(tools::file_path_sans_ext(basename(raster_path)))
  candidate_names <- tolower(basename(candidates))
  keep <- candidate_names == raster_name |
    startsWith(candidate_names, paste0(raster_name, ".")) |
    candidate_names %in% paste0(stem, c(".tfw", ".tfwx", ".wld", ".prj"))
  candidates[keep & !dir.exists(candidates)]
}
