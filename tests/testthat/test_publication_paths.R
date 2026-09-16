test_that("publication paths put each analysis figure in its own folder", {
  source_function("publication_output_paths.R")

  config <- list(
    table_dir = file.path("outputs", "tables"),
    figure_dir = file.path("outputs", "figures"),
    diagnostic_dir = file.path("outputs", "diagnostics")
  )
  paths <- publication_output_paths(config)

  expect_named(paths, c("lcbd", "regional", "gllvm"))
  expect_length(paths$lcbd, 4L)
  expect_length(paths$regional, 10L)
  expect_length(paths$gllvm, 10L)
  expect_equal(sum(startsWith(paths$lcbd,
                              file.path(config$figure_dir, "lcbd_plots"))), 1L)
  expect_equal(sum(startsWith(paths$regional,
                              file.path(config$figure_dir, "regional_plots"))), 5L)
  expect_equal(sum(startsWith(paths$gllvm,
                              file.path(config$figure_dir, "gllvm_plots"))), 5L)
  expect_true(file.path(config$diagnostic_dir,
                        "06_gllvm_basic_diagnostics.pdf") %in% paths$gllvm)
  expect_false(anyDuplicated(unlist(paths)) > 0L)
})
