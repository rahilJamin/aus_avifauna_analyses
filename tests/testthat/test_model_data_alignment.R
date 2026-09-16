test_that("LCBD preparation aligns character IDs and removes Noisy Miner", {
  source_function("assert_unique_ids.R")
  source_function("calculate_lcbd_beta_matrices.R")
  source_function("prepare_lcbd_model_data.R")

  species_matrix <- matrix(
    c(
      1, 1, 0,
      0, 0, 1,
      1, 1, 1
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      c("hex_01", "hex_02", "hex_03"),
      c("Manorina melanocephala", "sp_a", "sp_b")
    )
  )
  hex_metadata <- data.frame(
    hex_id = c("hex_01", "hex_02", "hex_03"),
    S_obs = c(2, 1, 3),
    sampling_units = c(6, 6, 6),
    n_checklists = c(10, 5, 9),
    miner_detection_rate = c(0.2, 0, 0.4)
  )
  landscape_data <- data.frame(
    hex_id = c("hex_01", "hex_02", "hex_03"),
    pland_urban_2500m = c(5, 8, 99),
    pland_agriculture_2500m = c(10, 20, 30),
    area_mn_woodland_2500m = c(100, 200, 300),
    cohesion_woodland_2500m = c(40, 50, 60),
    unused_metric_with_na = c(1, NA, 3)
  )
  grid <- sf::st_sf(
    hex_id = c("hex_01", "hex_02", "hex_03"),
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
      ))),
      sf::st_polygon(list(matrix(
        c(2, 0, 3, 0, 3, 1, 2, 1, 2, 0),
        ncol = 2,
        byrow = TRUE
      )))
    ),
    crs = 4326
  )

  hex_path <- tempfile(fileext = ".gpkg")
  on.exit(unlink(hex_path), add = TRUE)
  sf::st_write(grid, hex_path, quiet = TRUE)

  prepared <- prepare_lcbd_model_data(
    species_matrix = species_matrix,
    hex_metadata = hex_metadata,
    landscape_data = landscape_data,
    hex_path = hex_path,
    noisy_miner_species = "Manorina melanocephala",
    urban_cover_max = 10,
    model_coord_crs = 3857
  )

  expect_identical(prepared$model_df$hex_id, c("hex_01", "hex_02"))
  expect_identical(
    rownames(prepared$lcbd_species_matrix),
    prepared$model_df$hex_id
  )
  expect_false(
    "Manorina melanocephala" %in% colnames(prepared$lcbd_species_matrix)
  )
  expect_true(all(prepared$lcbd_species_matrix %in% c(0, 1)))
  expect_true(all(is.finite(prepared$model_df$replacement_ratio)))
  expect_true("log_effort" %in% names(prepared$model_df))
  expect_true("unused_metric_with_na" %in% names(prepared$model_df))
  expect_equal(
    prepared$attrition$value[
      prepared$attrition$step == "final_model_sites"
    ],
    nrow(prepared$model_df)
  )
  expect_identical(
    dimnames(prepared$beta_matrices$dist_mats$total),
    list(prepared$model_df$hex_id, prepared$model_df$hex_id)
  )
})

test_that("LCBD preparation rejects duplicate IDs before any spatial join", {
  source_function("assert_unique_ids.R")
  source_function("calculate_lcbd_beta_matrices.R")
  source_function("prepare_lcbd_model_data.R")

  species_matrix <- matrix(
    c(1, 0, 0, 1),
    nrow = 2,
    dimnames = list(c("h1", "h2"), c("sp_a", "sp_b"))
  )
  duplicated_metadata <- data.frame(
    hex_id = c("h1", "h1"),
    S_obs = 1,
    sampling_units = 6,
    n_checklists = 6,
    miner_detection_rate = 0
  )
  landscape <- data.frame(
    hex_id = c("h1", "h2"),
    pland_urban_2500m = 0,
    pland_agriculture_2500m = c(0, 1),
    area_mn_woodland_2500m = c(1, 2),
    cohesion_woodland_2500m = c(2, 3)
  )

  expect_error(
    prepare_lcbd_model_data(
      species_matrix,
      duplicated_metadata,
      landscape,
      tempfile(fileext = ".gpkg"),
      noisy_miner_species = "Manorina melanocephala",
      urban_cover_max = 10,
      model_coord_crs = 3857
    ),
    "duplicated"
  )
})
