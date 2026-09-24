test_that("CLUM and NVIS maps preserve category-table and water handling", {
  source_config()
  source_function("make_landscape_reclassification.R")
  rules <- cfg$landscape$classification
  clum <- data.frame(Value = 1001:1020, SIMPN = 1:20)
  mapped <- make_landscape_reclassification(clum, "SIMPN", rules$landuse_groups)
  expect_equal(mapped[, "from"], 1001:1020)
  expect_equal(unname(mapped[, "to"]),
               c(1, 1, 1, 2, 2, 3, 4, 4, 4, 4, 4, 4, 4, 5, 4, 5, 5, 5, NA, NA))
  lower_case <- clum
  names(lower_case) <- tolower(names(lower_case))
  expect_equal(make_landscape_reclassification(lower_case, "SIMPN", rules$landuse_groups), mapped)

  nvis <- data.frame(Value = c(1, 15, 29, 31, 99, 24, 28, 30))
  mapped <- make_landscape_reclassification(nvis, "Value", rules$vegetation_groups)
  expect_equal(unname(mapped[, "to"]), c(1, 1, 2, 3, 4, NA, NA, NA))
  # Value 30 is in the category table but unmapped: it stays NA. The fallback
  # class 5 is reserved for raster codes absent from the category table entirely.
  expect_error(make_landscape_reclassification(clum[c(1, 1), ], "SIMPN",
                                               rules$landuse_groups), "unique")
  expect_error(make_landscape_reclassification(clum["Value"], "SIMPN",
                                               rules$landuse_groups), "SIMPN")
})

test_that("raster tracking includes classification sidecars without other rasters", {
  source_function("track_raster_files.R")
  directory <- tempfile("landscape_sidecars_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  names <- c("map.tif", "map.tif.aux.xml", "map.tif.vat.dbf", "map.tfw",
             "map_other.tif", "map_other.tif.aux.xml")
  invisible(file.create(file.path(directory, names)))
  tracked <- track_raster_files(file.path(directory, "map.tif"))
  expect_setequal(basename(tracked), names[1:4])
})

test_that("landscape extraction follows current Chao2 sites and the shared grid", {
  # Capture commands without executing any target or opening research inputs.
  # Changes to site selection must flow into all three model branches.
  env <- new.env(parent = baseenv())
  source_config(envir = env)
  env$tar_target <- function(name, command, ...) {
    list(name = as.character(substitute(name)), command = substitute(command))
  }
  for (file in c("pipeline_data_processing.R", "pipeline_landscape.R",
                 "pipeline_lcbd.R", "pipeline_regional_convergence.R",
                 "pipeline_gllvm.R")) {
    source(file.path(find_repo_root(), "pipeline", file), local = env)
  }
  targets <- c(env$data_processing_targets, env$landscape_targets,
               env$model_data_targets, env$lcbd_targets,
               env$regional_convergence_targets, env$gllvm_targets)
  commands <- setNames(
    lapply(targets, function(target) target$command),
    vapply(targets, function(target) target$name, character(1))
  )
  expect_identical(anyDuplicated(names(commands)), 0L)
  upstream <- function(name, chain = character()) {
    if (name %in% chain) stop("Cyclic target dependency: ", name)
    dependencies <- intersect(all.vars(commands[[name]]), names(commands))
    unique(c(dependencies, unlist(lapply(dependencies, upstream, chain = c(chain, name)))))
  }
  expect_true(all(c("landscape_site_ids", "inext_matrix_bundle", "hex_grid_files",
                    "hex_grid", "raw_ala_file", "landuse_raster_files",
                    "vegetation_raster_files", "study_area_files") %in%
                    upstream("landscape_metrics_file")))
  expect_false("inext_matrix_bundle" %in% upstream("landscape_reclassified_files"))
  for (model in c("lcbd_model_bundle", "regional_convergence_results",
                  "gllvm_environment_model_file", "gllvm_fourth_corner_model_file")) {
    expect_true(all(c("landscape_metrics_file", "inext_matrix_bundle", "hex_grid") %in%
                      upstream(model)))
  }

  # The small derived ID target must preserve the complete retained set and ID
  # text, not intersect with a supplied historical list or coerce IDs to numbers.
  env$inext_matrix_bundle <- list(valid = data.frame(hex_id = c("hex_0000101", "001")))
  expect_identical(eval(commands$landscape_site_ids, env), c("hex_0000101", "001"))
  expect_false("landscape_site_ids_file" %in% names(commands))
  expect_null(env$cfg$landscape$site_ids_path)
  expect_identical(env$grid_settings, env$cfg$hex_grid)
})

make_landscape_test_jobs <- function() {
  # Reverse site order between sources to expose accidental positional joins.
  vegetation <- expand.grid(plot_id = c("b", "a", "c", "d"), class = 1:3,
                             metric = c("pland", "area_mn"), stringsAsFactors = FALSE)
  vegetation$value <- ifelse(vegetation$metric == "pland", 10 * vegetation$class, 4)
  vegetation <- vegetation[!(vegetation$plot_id == "b" & vegetation$class == 2), ]
  vegetation$value[vegetation$plot_id == "c"] <- NA_real_
  woodland <- expand.grid(plot_id = c("a", "b", "c", "d"), class = 1,
                           metric = c("ed", "cohesion", "area_mn"),
                           stringsAsFactors = FALSE)
  woodland$value <- c(1:4, 80:83, 100:103)
  woodland$value[woodland$plot_id == "d" & woodland$metric == "cohesion"] <- NA_real_
  landuse <- expand.grid(plot_id = c("d", "c", "b", "a"), class = 1:5,
                         metric = "pland", stringsAsFactors = FALSE)
  landuse$value <- 20
  landuse$value[landuse$plot_id == "c"] <- NA_real_
  list(vegetation = vegetation, woodland = woodland, landuse = landuse)
}

test_that("landscape assembly aligns sites and preserves zero and NA rules", {
  source_config()
  source_function("assemble_landscape_metrics.R")
  jobs <- make_landscape_test_jobs()
  result <- assemble_landscape_metrics(jobs$vegetation, jobs$woodland,
                                      jobs$landuse, cfg$landscape$extraction)
  selected <- result$selected_metrics
  expect_identical(selected$hex_id, c("b", "a"))
  expect_equal(selected$area_mn_woodland_2500m, c(101, 100))
  expect_equal(selected$pland_sclerophyll_2500m, c(0, 20))
  expect_equal(selected$PLAND_habitat_2500m, c(40, 60))
  expect_named(selected, c("hex_id", paste0(cfg$landscape$extraction$selected_metrics, "_2500m")))
  expect_true(is.na(result$all_metrics$pland_wet_dense_2500m[3]))
  expect_true(is.na(result$all_metrics$cohesion_woodland_2500m[4]))

  # When land-use coverage confirms that a site was sampled, an absent
  # vegetation class has zero PLAND; missing patch metrics remain missing.
  jobs$landuse$value[jobs$landuse$plot_id == "c" & jobs$landuse$class == 4] <- 100
  result <- assemble_landscape_metrics(jobs$vegetation, jobs$woodland,
                                      jobs$landuse, cfg$landscape$extraction)
  expect_equal(result$selected_metrics$PLAND_habitat_2500m[
    result$selected_metrics$hex_id == "c"], 0)

  jobs$woodland <- rbind(jobs$woodland, jobs$woodland[1, ])
  expect_error(assemble_landscape_metrics(jobs$vegetation, jobs$woodland,
                                         jobs$landuse, cfg$landscape$extraction), "Duplicate")
})

test_that("a vegetation class absent throughout a subset still has a PLAND column", {
  source_config()
  source_function("assemble_landscape_metrics.R")
  jobs <- make_landscape_test_jobs()
  jobs$vegetation <- jobs$vegetation[jobs$vegetation$class != 3, ]
  result <- assemble_landscape_metrics(jobs$vegetation, jobs$woodland,
                                      jobs$landuse, cfg$landscape$extraction)
  expect_equal(result$selected_metrics$pland_mallee_arid_2500m, c(0, 0))
  expect_true(is.na(result$all_metrics$pland_mallee_arid_2500m[3]))
})

test_that("landscape input and output paths are portable and separated", {
  source_config()
  source_function("safe_output_writers.R")
  paths <- cfg$landscape$paths
  generated <- unlist(paths[c("vegetation_classes", "vegetation_binary", "landuse_classes",
                              "selected_metrics")])
  expect_false(any(vapply(generated, path_is_inside, logical(1), parent = cfg$raw_data_dir)))
  expect_false(any(grepl("^([A-Za-z]:|/)", unlist(paths))))
  expect_true(all(vapply(paths[c("landuse", "vegetation")],
                         path_is_inside, logical(1), parent = cfg$raw_data_dir)))
  expect_true(path_is_inside(cfg$study_area_path, cfg$raw_data_dir))
})

test_that("current terra classification requires an explicit source-NA guard", {
  skip_if_not_installed("terra")
  source_config()
  source_function("make_landscape_reclassification.R")
  raster <- terra::rast(nrows = 1, ncols = 4, xmin = 0, xmax = 4, ymin = 0, ymax = 1)
  terra::values(raster) <- c(1, 30, 100, NA)
  mapping <- make_landscape_reclassification(data.frame(Value = c(1, 30)),
                                             "Value", cfg$landscape$classification$vegetation_groups)
  result <- terra::classify(raster, mapping, others = 5)
  # Current terra applies the fallback to source NA cells. The reconstruction
  # function must therefore restore the source raster's missing-cell mask.
  expect_equal(as.numeric(terra::values(result)), c(1, NA, 5, 5))
})

test_that("reclassification preserves categories and missing cells", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  source_config()
  source_function("safe_output_writers.R")
  source_function("make_landscape_reclassification.R")
  source_function("reclassify_landscape_rasters.R")
  directory <- tempfile("landscape_reclass_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  rules <- cfg$landscape$classification
  raster <- terra::rast(nrows = 4, ncols = 4, xmin = 500000, xmax = 500200,
                        ymin = 6000000, ymax = 6000200, crs = rules$crs)
  landuse <- raster
  terra::values(landuse) <- rep(101:104, 4)
  levels(landuse) <- data.frame(Value = 101:104, SIMPN = c(1, 7, 19, 20))
  vegetation <- raster
  terra::values(vegetation) <- rep(c(1, 30, 100, NA), 4)
  levels(vegetation) <- data.frame(Value = c(1, 30), MVG_NAME = c("forest", "unmapped"))
  landuse_path <- file.path(directory, "clum.tif")
  vegetation_path <- file.path(directory, "nvis.tif")
  area_path <- file.path(directory, "area.gpkg")
  terra::writeRaster(landuse, landuse_path)
  terra::writeRaster(vegetation, vegetation_path)
  polygon <- sf::st_as_sfc(sf::st_bbox(c(xmin = 500000, xmax = 500200,
                                        ymin = 6000000, ymax = 6000200), crs = sf::st_crs(rules$crs)))
  sf::st_write(sf::st_sf(id = 1, geometry = polygon), area_path, quiet = TRUE)
  output_paths <- list(
    vegetation_classes = file.path(directory, "new_veg.tif"),
    vegetation_binary = file.path(directory, "new_binary.tif"),
    landuse_classes = file.path(directory, "new_lu.tif")
  )
  reclassify_landscape_rasters(landuse_path, vegetation_path, area_path, rules, output_paths)
  read_values <- function(path) as.numeric(terra::values(terra::rast(path)))
  expect_equal(read_values(output_paths$vegetation_classes), rep(c(1, NA, 5, NA), 4))
  expect_equal(read_values(output_paths$vegetation_binary), rep(c(1, NA, 2, NA), 4))
  expect_equal(read_values(output_paths$landuse_classes), rep(c(1, 4, NA, NA), 4))
  output_paths$landuse_classes <- landuse_path
  expect_error(reclassify_landscape_rasters(landuse_path, vegetation_path, area_path,
                                            rules, output_paths), "overwrite an input")
})

test_that("landscape sampling uses projected centroids and gives 100 percent for a uniform raster", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  skip_if_not_installed("landscapemetrics")
  source_function("assert_unique_ids.R")
  source_function("prepare_landscape_sites.R")
  source_function("extract_landscape_metrics.R")
  square <- function(x) sf::st_polygon(list(matrix(
    c(x, 100, x + 100, 100, x + 100, 200, x, 200, x, 100), ncol = 2, byrow = TRUE
  )))
  grid <- sf::st_sf(hex_id = c("west", "east"),
                    geometry = sf::st_sfc(square(100), square(300), crs = 3577))
  grid_path <- tempfile(fileext = ".gpkg")
  raster_path <- tempfile(fileext = ".tif")
  on.exit(unlink(c(grid_path, raster_path)), add = TRUE)
  sf::st_write(grid, grid_path, quiet = TRUE)
  sites <- prepare_landscape_sites(grid_path, c("east", "west"), 3577)
  expect_identical(sites$hex_id, c("west", "east"))
  expect_equal(unname(sf::st_coordinates(sites)), matrix(c(150, 350, 150, 150), ncol = 2))
  expect_error(prepare_landscape_sites(grid_path, "missing", 3577), "absent")
  expect_error(prepare_landscape_sites(grid_path, c("west", "west"), 3577),
               "duplicated")
  expect_error(prepare_landscape_sites(grid_path, character(), 3577),
               "At least one")
  expect_error(prepare_landscape_sites(grid_path, "west", 4326), "projected")

  raster <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 500,
                        ymin = 0, ymax = 500, crs = "EPSG:3577")
  terra::values(raster) <- 1
  terra::writeRaster(raster, raster_path)
  result <- extract_landscape_metrics(raster_path, sites, 100, "lsm_c_pland", FALSE)
  expect_setequal(as.character(result$plot_id), sites$hex_id)
  expect_equal(result$value, c(100, 100))
})
