test_that("shapefile tracking includes required and optional components", {
  source_function("track_shapefile_files.R")

  directory <- tempfile("shapefile_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)

  paths <- file.path(directory, paste0("grid.", c("shp", "shx", "dbf", "prj")))
  expect_true(all(file.create(paths)))

  tracked <- track_shapefile_files(file.path(directory, "grid.shp"))

  expect_setequal(tracked, paths)
  expect_identical(main_shapefile_path(tracked), file.path(directory, "grid.shp"))
})

test_that("shapefile tracking rejects incomplete component sets", {
  source_function("track_shapefile_files.R")

  directory <- tempfile("shapefile_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)

  expect_true(file.create(file.path(directory, "grid.shp")))
  expect_true(file.create(file.path(directory, "grid.shx")))

  expect_error(
    track_shapefile_files(file.path(directory, "grid.shp")),
    "[.]dbf"
  )
})

test_that("tracked files are selected by path without relying on names", {
  source_function("track_shapefile_files.R")

  tracked <- unname(c(
    file.path("outputs", "table.csv"),
    file.path("outputs", "matrix.rds")
  ))

  expect_identical(
    tracked_file_path(tracked, file.path("outputs", "matrix.rds")),
    file.path("outputs", "matrix.rds")
  )
  expect_error(
    tracked_file_path(tracked, file.path("outputs", "missing.rds")),
    "exactly one"
  )
})

test_that("IBRA assignment preserves sites and rejects multiple matches", {
  source_function("assert_unique_ids.R")
  source_function("assign_ibra_subregions.R")

  ibra <- sf::st_sf(
    SUB_NAME_7 = c("west", "east"),
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

  ibra_path <- tempfile(fileext = ".gpkg")
  on.exit(unlink(ibra_path), add = TRUE)
  sf::st_write(ibra, ibra_path, quiet = TRUE)

  model_df <- data.frame(
    hex_id = c("h1", "h2"),
    x = c(0.5, 1.5),
    y = c(0.5, 0.5)
  )

  assigned <- assign_ibra_subregions(model_df, ibra_path, 4326)
  expect_identical(assigned$hex_id, model_df$hex_id)
  expect_identical(assigned$ibra_sub, c("west", "east"))

  overlapping_ibra <- ibra
  overlapping_ibra$geometry[[2]] <- sf::st_polygon(list(matrix(
    c(0.25, 0, 1.75, 0, 1.75, 1, 0.25, 1, 0.25, 0),
    ncol = 2,
    byrow = TRUE
  )))
  overlapping_path <- tempfile(fileext = ".gpkg")
  on.exit(unlink(overlapping_path), add = TRUE)
  sf::st_write(overlapping_ibra, overlapping_path, quiet = TRUE)

  expect_error(
    assign_ibra_subregions(model_df[1, ], overlapping_path, 4326),
    "exactly one"
  )
})

test_that("dbMEM selection validates alignment before computation", {
  source_function("select_gllvm_dbmem.R")

  species_matrix <- matrix(c(1, 0, 0, 1), nrow = 2)

  expect_error(
    select_gllvm_dbmem(
      species_matrix,
      coordinates = matrix(1:6, ncol = 2),
      seed = 1,
      permutations = 99,
      global_alpha = 0.05,
      selection_alpha = 0.01
    ),
    "aligned"
  )

  expect_error(
    select_gllvm_dbmem(
      species_matrix,
      coordinates = matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE),
      seed = 1,
      permutations = 99,
      global_alpha = 2,
      selection_alpha = 0.01
    ),
    "between zero and one"
  )
})
