#' Rebuild the vegetation and land-use rasters used for landscape metrics
#'
#' Crops and masks each source in its own CRS, applies the configured category
#' mappings, and projects each result separately using nearest-neighbour
#' resampling. The binary vegetation raster is derived before projection.
#'
#' @param landuse_path Source CLUM GeoTIFF, including its category sidecars.
#' @param vegetation_path Source NVIS GeoTIFF, including its category sidecars.
#' @param study_area_path Study-area polygon file used to mask both rasters.
#' @param classification Classification rules and CRS from cfg$landscape.
#' @param output_paths Named paths: vegetation_classes, vegetation_binary,
#'   landuse_classes. Outputs must be outside raw-data directories.
#' @return Paths of the three generated GeoTIFFs. Raster objects themselves are
#'   not returned because terra objects contain pointers unsuitable for RDS caches.
reclassify_landscape_rasters <- function(landuse_path, vegetation_path,
                                        study_area_path, classification,
                                        output_paths) {
  output_names <- c("vegetation_classes", "vegetation_binary", "landuse_classes")
  if (!all(output_names %in% names(output_paths))) {
    stop("Missing a reclassified raster output path.", call. = FALSE)
  }
  output_paths <- unlist(output_paths[output_names], use.names = TRUE)
  invisible(lapply(output_paths, assert_not_raw_output_path))
  if (anyDuplicated(vapply(output_paths, normalise_path_for_check, character(1)))) {
    stop("Reclassified rasters need distinct output paths.", call. = FALSE)
  }
  input_paths <- c(landuse_path, vegetation_path, study_area_path)
  if (any(vapply(output_paths, normalise_path_for_check, character(1)) %in%
          vapply(input_paths, normalise_path_for_check, character(1)))) {
    stop("Reclassified raster outputs must not overwrite an input.", call. = FALSE)
  }

  message("Cropping and masking CLUM and NVIS rasters")
  study_area <- sf::st_read(study_area_path, quiet = TRUE)
  landuse <- terra::rast(landuse_path)
  vegetation <- terra::rast(vegetation_path)
  if (terra::nlyr(landuse) != 1L || terra::nlyr(vegetation) != 1L) {
    stop("Each source raster must have exactly one categorical layer.", call. = FALSE)
  }

  # Transform the study area into each source CRS before cropping and masking.
  study_area <- sf::st_transform(study_area, terra::crs(landuse))
  landuse <- terra::mask(terra::crop(landuse, study_area), study_area)
  study_area <- sf::st_transform(study_area, terra::crs(vegetation))
  vegetation <- terra::mask(terra::crop(vegetation, study_area), study_area)

  landuse_matrix <- make_landscape_reclassification(
    terra::cats(landuse)[[1]], "SIMPN", classification$landuse_groups
  )
  vegetation_matrix <- make_landscape_reclassification(
    terra::cats(vegetation)[[1]], "Value", classification$vegetation_groups
  )
  landuse <- terra::classify(landuse, landuse_matrix,
                            others = classification$landuse_other)

  # Preserve the source missing-cell mask explicitly. terra applies `others = 5`
  # to source NA cells as well as unknown numeric codes. Without this guard,
  # missing cells become nonhabitat and alter PLAND and cohesion denominators.
  vegetation_missing <- is.na(vegetation)
  vegetation <- terra::classify(vegetation, vegetation_matrix,
                               others = classification$vegetation_other)
  vegetation <- terra::mask(
    vegetation, vegetation_missing, maskvalues = 1, updatevalue = NA
  )

  # Derive woodland/non-woodland after restoring the source NA footprint.
  # Explicitly mapped water/sea categories are already NA after classification.
  binary_matrix <- cbind(
    as.numeric(names(classification$vegetation_binary)),
    unname(classification$vegetation_binary)
  )
  vegetation_binary <- terra::classify(vegetation, binary_matrix)

  # Project each classified raster independently. terra chooses the output
  # resolution and extent for each source; nearest-neighbour resampling retains
  # categorical values.
  rasters <- list(vegetation, vegetation_binary, landuse)
  for (i in seq_along(rasters)) {
    message("Projecting and saving ", output_names[i])
    projected <- terra::project(rasters[[i]], classification$crs,
                                method = classification$resampling)
    dir.create(dirname(output_paths[i]), recursive = TRUE, showWarnings = FALSE)
    terra::writeRaster(projected, output_paths[i], overwrite = TRUE)
  }
  unname(output_paths)
}
