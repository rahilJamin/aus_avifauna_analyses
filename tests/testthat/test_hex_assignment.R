make_test_hex_grid <- function(hex_ids = c("h01", "h02")) {
  sf::st_sf(
    hex_id = hex_ids,
    geometry = sf::st_sfc(
      sf::st_polygon(list(matrix(
        c(0, 0, 1, 0, 1, 1, 0, 1, 0, 0),
        ncol = 2,
        byrow = TRUE
      ))),
      sf::st_polygon(list(matrix(
        c(1, 0, 2, 0, 2, 1, 1, 1, 1, 0),
        ncol = 2,
        byrow = TRUE
      )))
    ),
    crs = 4326
  )
}

make_overlapping_test_hex_grid <- function() {
  sf::st_sf(
    hex_id = c("h01", "h02"),
    geometry = sf::st_sfc(
      sf::st_polygon(list(matrix(
        c(0, 0, 1.5, 0, 1.5, 1, 0, 1, 0, 0),
        ncol = 2,
        byrow = TRUE
      ))),
      sf::st_polygon(list(matrix(
        c(0.5, 0, 2, 0, 2, 1, 0.5, 1, 0.5, 0),
        ncol = 2,
        byrow = TRUE
      )))
    ),
    crs = 4326
  )
}

test_that("hex assignment retains input order, character IDs, and inside records", {
  source_function("assert_unique_ids.R")
  source_function("assign_records_to_hex.R")

  hex_path <- tempfile(fileext = ".gpkg")
  on.exit(unlink(hex_path), add = TRUE)
  sf::st_write(make_test_hex_grid(), hex_path, quiet = TRUE)

  cleaned_records <- data.frame(
    record_name = c("inside_second", "outside", "inside_first"),
    longitude = c(1.5, 3, 0.5),
    latitude = c(0.5, 3, 0.5),
    species = c("sp_b", "sp_c", "sp_a")
  )

  assigned <- assign_records_to_hex(
    cleaned_records = cleaned_records,
    hex_path = hex_path,
    input_lonlat_crs = 4326
  )

  expect_identical(assigned$record_name, c("inside_second", "inside_first"))
  expect_identical(assigned$hex_id, c("h02", "h01"))
  expect_false("outside" %in% assigned$record_name)
  expect_false("record_id" %in% names(assigned))
})

test_that("hex assignment refuses overlapping cells and duplicated grid IDs", {
  source_function("assert_unique_ids.R")
  source_function("assign_records_to_hex.R")

  overlap_path <- tempfile(fileext = ".gpkg")
  duplicate_path <- tempfile(fileext = ".gpkg")
  on.exit(unlink(c(overlap_path, duplicate_path)), add = TRUE)
  sf::st_write(make_overlapping_test_hex_grid(), overlap_path, quiet = TRUE)
  sf::st_write(
    make_test_hex_grid(hex_ids = c("same", "same")),
    duplicate_path,
    quiet = TRUE
  )

  overlapping_record <- data.frame(
    longitude = 1,
    latitude = 0.5,
    species = "sp_a"
  )
  expect_error(
    assign_records_to_hex(overlapping_record, overlap_path, 4326),
    "more than one hexagon"
  )

  inside_record <- data.frame(
    longitude = 0.5,
    latitude = 0.5,
    species = "sp_a"
  )
  expect_error(
    assign_records_to_hex(inside_record, duplicate_path, 4326),
    "duplicated"
  )
})

test_that("hex assignment validates coordinate columns before spatial work", {
  source_function("assert_unique_ids.R")
  source_function("assign_records_to_hex.R")

  expect_error(
    assign_records_to_hex(
      data.frame(longitude = "1", latitude = 2),
      tempfile(fileext = ".gpkg"),
      4326
    ),
    "must be numeric"
  )
  expect_error(
    assign_records_to_hex(
      data.frame(longitude = NA_real_, latitude = 2),
      tempfile(fileext = ".gpkg"),
      4326
    ),
    "present and finite"
  )
})
