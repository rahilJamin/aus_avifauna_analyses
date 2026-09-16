test_that("path guards recognise protected paths without prefix mistakes", {
  source_config()
  source_function("safe_output_writers.R")

  raw_root <- cfg$raw_data_dir

  expect_true(path_is_inside(raw_root, raw_root))
  expect_true(path_is_inside(
    file.path(raw_root, "species", "ala.csv"),
    raw_root
  ))
  expect_true(path_is_inside(
    file.path(raw_root, "species", "..", "spatial_grids", "grid.gpkg"),
    raw_root
  ))
  expect_false(path_is_inside(
    paste0(raw_root, "_backup", .Platform$file.sep, "copy.csv"),
    raw_root
  ))
  expect_false(path_is_inside(
    file.path(cfg$output_dir, "tables", "summary.csv"),
    raw_root
  ))

  expect_error(
    assert_not_raw_output_path(file.path(raw_root, "bad.csv")),
    "Refusing to write"
  )
  expect_identical(
    assert_not_raw_output_path(file.path(cfg$table_dir, "ok.csv")),
    file.path(cfg$table_dir, "ok.csv")
  )
})

test_that("path normalisation rejects invalid scalar inputs", {
  source_function("safe_output_writers.R")

  expect_error(normalise_path_for_check(character()), "one non-empty")
  expect_error(normalise_path_for_check(c("a", "b")), "one non-empty")
  expect_error(normalise_path_for_check(NA_character_), "one non-empty")
  expect_error(normalise_path_for_check(""), "one non-empty")
})
