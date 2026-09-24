make_grid_test_extent <- function() {
  sf::st_as_sfc(sf::st_bbox(
    c(xmin = 0, ymin = 0, xmax = 6000, ymax = 4000),
    crs = sf::st_crs(3577)
  ))
}

test_that("hexagons have the specified metric size and full boundary-cell areas", {
  skip_if_not_installed("sf")
  source_function("create_hex_grid.R")
  area <- make_grid_test_extent()
  grid <- create_hex_grid(area, 2500, 3577, FALSE)

  expect_equal(sf::st_crs(grid)$epsg, 3577)
  expect_false(sf::st_is_longlat(grid))
  expect_gt(nrow(grid), 0)
  expect_type(grid$hex_id, "character")
  expect_identical(anyDuplicated(grid$hex_id), 0L)
  expect_true(all(lengths(sf::st_intersects(grid, area)) > 0L))

  # For opposite-edge spacing d, the area is sqrt(3)/2 * d^2. This detects
  # radius/diameter confusion, degree-based sizes, and clipped boundary cells.
  expected_area <- sqrt(3) / 2 * 2500^2
  expect_equal(as.numeric(sf::st_area(grid)), rep(expected_area, nrow(grid)),
               tolerance = 1e-8)
  expect_gt(sum(as.numeric(sf::st_area(grid))), as.numeric(sf::st_area(area)))
  expect_equal(sum(as.numeric(sf::st_area(grid))),
               as.numeric(sf::st_area(sf::st_union(grid))), tolerance = 1e-8)

  box <- sf::st_bbox(grid[1, ])
  expect_equal(as.numeric(box["xmax"] - box["xmin"]), 2500, tolerance = 1e-8)
  expect_equal(as.numeric(box["ymax"] - box["ymin"]), 2 * 2500 / sqrt(3),
               tolerance = 1e-8)
  flat <- create_hex_grid(area, 2500, 3577, TRUE)
  flat_box <- sf::st_bbox(flat[1, ])
  expect_equal(as.numeric(flat_box["ymax"] - flat_box["ymin"]), 2500,
               tolerance = 1e-8)
  again <- create_hex_grid(area, 2500, 3577, FALSE)
  expect_identical(grid$hex_id, again$hex_id)
  expect_identical(sf::st_geometry(grid), sf::st_geometry(again))
})

test_that("longitude-latitude extents are transformed before making the grid", {
  skip_if_not_installed("sf")
  source_function("create_hex_grid.R")
  area <- sf::st_as_sfc(sf::st_bbox(
    c(xmin = 145, ymin = -35.05, xmax = 145.05, ymax = -35),
    crs = sf::st_crs(4326)
  ))
  grid <- create_hex_grid(area, 2500, 3577, FALSE)

  expect_equal(as.numeric(sf::st_area(grid)),
               rep(sqrt(3) / 2 * 2500^2, nrow(grid)), tolerance = 1e-8)
  centres <- sf::st_transform(sf::st_centroid(sf::st_geometry(grid)), 4326)
  xy <- sf::st_coordinates(centres)
  expect_true(all(xy[, 1] > 144.9 & xy[, 1] < 145.15))
  expect_true(all(xy[, 2] > -35.15 & xy[, 2] < -34.9))
})

test_that("disconnected extent polygons share one grid and never duplicate sites", {
  skip_if_not_installed("sf")
  source_function("create_hex_grid.R")
  left <- sf::st_as_sfc(sf::st_bbox(
    c(xmin = 0, ymin = 0, xmax = 1000, ymax = 1000), crs = sf::st_crs(3577)
  ))
  right <- sf::st_as_sfc(sf::st_bbox(
    c(xmin = 9000, ymin = 0, xmax = 10000, ymax = 1000), crs = sf::st_crs(3577)
  ))
  area <- c(left, right)
  grid <- create_hex_grid(area, 500, 3577, FALSE)
  expect_true(all(lengths(sf::st_intersects(grid, area)) > 0L))
  expect_identical(anyDuplicated(grid$hex_id), 0L)
  expect_true(any(diff(as.integer(sub("^hex_", "", grid$hex_id))) > 1L))

  # Repeating an extent polygon must not create repeated sampling sites.
  repeated <- create_hex_grid(c(area, left), 500, 3577, FALSE)
  expect_identical(repeated$hex_id, grid$hex_id)
  expect_identical(sf::st_geometry(repeated), sf::st_geometry(grid))
})

test_that("grid construction rejects missing CRSs, angular sizes and invalid extents", {
  skip_if_not_installed("sf")
  source_function("create_hex_grid.R")
  area <- make_grid_test_extent()
  expect_error(create_hex_grid(data.frame(x = 1), 2500, 3577, FALSE), "sf or sfc")
  expect_error(create_hex_grid(sf::st_set_crs(area, NA), 2500, 3577, FALSE), "no CRS")
  expect_error(create_hex_grid(area, 2500, 4326, FALSE), "projected CRS in metres")
  expect_error(create_hex_grid(area, 2500, 2263, FALSE), "projected CRS in metres")
  for (size in list(0, -1, NA_real_, Inf, "2500", c(500, 1000))) {
    expect_error(create_hex_grid(area, size, 3577, FALSE), "positive finite")
  }
  expect_error(create_hex_grid(area, 2500, 3577, NA), "TRUE or FALSE")
  expect_error(create_hex_grid(sf::st_sfc(crs = 3577), 2500, 3577, FALSE), "non-empty")
  expect_error(create_hex_grid(sf::st_centroid(area), 2500, 3577, FALSE), "polygon")
  crossed <- sf::st_sfc(sf::st_polygon(list(matrix(
    c(0, 0, 100, 100, 0, 100, 100, 0, 0, 0), ncol = 2, byrow = TRUE
  ))), crs = 3577)
  expect_error(create_hex_grid(crossed, 2500, 3577, FALSE), "valid polygon")
})

test_that("occurrence assignment uses the generated projected grid", {
  skip_if_not_installed("sf")
  source_function("create_hex_grid.R")
  source_function("assert_unique_ids.R")
  source_function("assign_records_to_hex.R")
  source_function("prepare_landscape_sites.R")
  area <- sf::st_as_sfc(sf::st_bbox(
    c(xmin = 145, ymin = -35.05, xmax = 145.05, ymax = -35),
    crs = sf::st_crs(4326)
  ))
  grid <- create_hex_grid(area, 2500, 3577, FALSE)
  path <- tempfile(fileext = ".gpkg")
  on.exit(unlink(path), add = TRUE)
  sf::st_write(grid, path, quiet = TRUE)

  selected <- grid[c(nrow(grid), 1L), ]
  points <- sf::st_transform(sf::st_centroid(sf::st_geometry(selected)), 4326)
  xy <- sf::st_coordinates(points)
  records <- data.frame(longitude = xy[, 1], latitude = xy[, 2], species = "toy_bird")
  assigned <- assign_records_to_hex(records, path, 4326)
  expect_identical(assigned$hex_id, selected$hex_id)
  expect_equal(assigned$longitude, records$longitude)
  expect_equal(assigned$latitude, records$latitude)

  # Sampling in the raster projection must still place each centre inside its
  # own original grid cell. This detects relabelled rather than transformed CRS.
  sampled <- prepare_landscape_sites(path, selected$hex_id, 32655)
  expect_equal(sf::st_crs(sampled)$epsg, 32655)
  containing_cells <- sf::st_within(sf::st_transform(sampled, 3577), grid)
  expect_true(all(lengths(containing_cells) == 1L))
  expect_identical(grid$hex_id[unlist(containing_cells)], sampled$hex_id)
})
