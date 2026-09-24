#' Create the analysis hexagon grid from the study extent
#'
#' Transforms the supplied extent to a projected CRS in metres, creates a regular
#' hexagonal grid, and retains whole cells intersecting the extent. Boundary
#' cells are not clipped, so every site has the same area in the grid CRS.
#'
#' @param study_area An sf or sfc polygon object with a declared CRS.
#' @param cellsize_m Distance between opposite hexagon edges in projected metres.
#'   This is neither the hexagon radius nor its area.
#' @param grid_crs Projected CRS in metres; the main analysis uses EPSG:3577
#'   (GDA94 / Australian Albers).
#' @param flat_topped Logical; FALSE creates point-topped hexagons.
#'
#' @return An sf polygon table in grid_crs with character hex_id values. IDs are
#'   assigned in full-grid order before spatial filtering. They are reproducible
#'   for a fixed extent, CRS, cell size, orientation and sf environment.
create_hex_grid <- function(study_area, cellsize_m, grid_crs, flat_topped) {
  if (!inherits(study_area, c("sf", "sfc"))) {
    stop("The study extent must be an sf or sfc polygon object.", call. = FALSE)
  }
  area <- sf::st_geometry(study_area)
  if (length(area) == 0L || any(sf::st_is_empty(area))) {
    stop("The study extent must contain non-empty polygons.", call. = FALSE)
  }
  if (is.na(sf::st_crs(area))) {
    stop("The study extent has no CRS; provide its correct source CRS.", call. = FALSE)
  }
  if (any(!as.character(sf::st_geometry_type(area)) %in% c("POLYGON", "MULTIPOLYGON")) ||
      !all(sf::st_is_valid(area) %in% TRUE)) {
    stop("The study extent must contain valid polygon geometries.", call. = FALSE)
  }
  if (!is.numeric(cellsize_m) || length(cellsize_m) != 1L ||
      is.na(cellsize_m) || !is.finite(cellsize_m) || cellsize_m <= 0) {
    stop("`cellsize_m` must be one positive finite number of metres.", call. = FALSE)
  }
  if (!is.logical(flat_topped) || length(flat_topped) != 1L || is.na(flat_topped)) {
    stop("`flat_topped` must be TRUE or FALSE.", call. = FALSE)
  }
  projected_crs <- sf::st_crs(grid_crs)
  if (is.na(projected_crs) || isTRUE(sf::st_is_longlat(projected_crs)) ||
      !isTRUE(tolower(projected_crs$units_gdal) %in%
                c("metre", "meter", "metres", "meters", "m"))) {
    stop("Grid construction requires a projected CRS in metres.", call. = FALSE)
  }

  # Transform coordinates before applying a size in metres. Assigning a CRS to
  # longitude/latitude coordinates would relabel them without projecting them.
  area <- sf::st_transform(area, projected_crs)
  bounds <- sf::st_bbox(area)
  if (any(!is.finite(bounds)) || bounds["xmax"] <= bounds["xmin"] ||
      bounds["ymax"] <= bounds["ymin"]) {
    stop("The projected study extent has no usable spatial extent.", call. = FALSE)
  }

  # Anchor the grid to the southwest corner of the projected study bounds.
  # A different extent can change both the grid origin and its identifiers.
  cells <- sf::st_make_grid(
    area,
    cellsize = cellsize_m,
    offset = unname(bounds[c("xmin", "ymin")]),
    square = FALSE,
    flat_topped = flat_topped
  )
  grid <- sf::st_sf(
    hex_id = sprintf("hex_%07d", seq_along(cells)),
    geometry = cells
  )

  # Spatial subsetting retains full cells and one row per cell, including when
  # the supplied study extent contains several polygons.
  grid <- sf::st_filter(grid, area, .predicate = sf::st_intersects)
  if (nrow(grid) == 0L) {
    stop("No generated hexagons intersect the study extent.", call. = FALSE)
  }
  rownames(grid) <- NULL
  grid
}
