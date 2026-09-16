test_that("incidence units use the selected survey-unit definition", {
  source_function("add_incidence_unit.R")

  records <- data.frame(
    hex_id = c("h1", "h1"),
    date = as.Date(c("2020-01-15", "2020-02-01")),
    dataResourceName = c("resource_a", "resource_b"),
    samplingProtocol = c("protocol_a", "protocol_b")
  )

  month_year <- add_incidence_unit(records, "month_year")
  expect_identical(month_year$incidence_unit, c("2020-01", "2020-02"))

  checklist <- add_incidence_unit(records, "checklist_event")
  expect_identical(
    checklist$incidence_unit,
    c(
      "h1__2020-01-15__resource_a__protocol_a",
      "h1__2020-02-01__resource_b__protocol_b"
    )
  )

  expect_error(add_incidence_unit(records, "unknown"), "Unknown incidence_unit")
  expect_error(
    add_incidence_unit(
      transform(records, date = as.Date(NA_character_)),
      "month_year"
    ),
    "complete and valid"
  )
  expect_error(
    add_incidence_unit(records[c("date")], "checklist_event"),
    "requires"
  )
})

test_that("Noisy Miner detection uses the same incidence units as Chao2", {
  source_function("add_incidence_unit.R")
  source_function("calculate_miner_detection.R")

  hex_occ <- data.frame(
    hex_id = c("h1", "h1", "h1", "h1", "h2"),
    date = as.Date(c(
      "2020-01-01", "2020-01-01", "2020-01-20", "2020-02-01",
      "2020-03-01"
    )),
    dataResourceName = "resource_a",
    samplingProtocol = "protocol_a",
    species = c(
      "Manorina melanocephala", "Other species", "Other species",
      "Other species", "Other species"
    )
  )

  by_month <- calculate_miner_detection(
    hex_occ,
    noisy_miner_species = "Manorina melanocephala",
    incidence_unit = "month_year"
  )
  by_event <- calculate_miner_detection(
    hex_occ,
    noisy_miner_species = "Manorina melanocephala",
    incidence_unit = "checklist_event"
  )

  by_month <- by_month[order(by_month$hex_id), ]
  by_event <- by_event[order(by_event$hex_id), ]

  expect_identical(by_month$n_checklists, c(2L, 1L))
  expect_equal(by_month$miner_detection_rate, c(0.5, 0))
  expect_identical(by_event$n_checklists, c(3L, 1L))
  expect_equal(by_event$miner_detection_rate, c(1 / 3, 0))
  expect_true(all(by_event$n_miner_checklists <= by_event$n_checklists))
})

test_that("Chao2 filtering retains equality at every configured cutoff", {
  source_function("build_matrix_bundle.R")

  chao <- data.frame(
    hex_id = c("equal", "too_few_units", "too_few_species", "incomplete", "failed"),
    chao_status = c("ok", "ok", "ok", "ok", "failed"),
    sampling_units = c(6, 5, 6, 6, 100),
    S_obs = c(5, 5, 4, 5, 100),
    chao2_completeness = c(0.7, 0.9, 0.9, 0.69, 1)
  )

  retained <- retain_sampled_hexagons(
    chao,
    min_sampling_units = 6,
    min_species = 5,
    min_completeness = 0.7
  )

  expect_identical(retained$hex_id, "equal")
  expect_error(
    retain_sampled_hexagons(chao, 0, 5, 0.7),
    "positive whole"
  )
  expect_error(
    retain_sampled_hexagons(chao, 6, 5, 1.1),
    "between zero and one"
  )
})

test_that("one-hexagon Chao2 estimates preserve incidence invariants", {
  source_function("add_incidence_unit.R")
  source_function("estimate_chao2_for_hex.R")

  hex_data <- data.frame(
    hex_id = "h1",
    date = as.Date(c(
      "2020-01-01", "2020-02-01", "2020-03-01",
      "2020-01-01", "2020-01-01", "2020-02-01"
    )),
    species = c("sp_a", "sp_a", "sp_a", "sp_b", "sp_c", "sp_c")
  )

  estimate <- estimate_chao2_for_hex(hex_data, "month_year")

  expect_identical(estimate$chao_status, "ok")
  expect_equal(estimate$S_obs, 3)
  expect_equal(estimate$sampling_units, 3)
  expect_equal(estimate$Q1, 1)
  expect_equal(estimate$Q2, 1)
  expect_gte(estimate$chao2_estimate, estimate$S_obs)
  expect_true(estimate$chao2_completeness > 0)
  expect_lte(estimate$chao2_completeness, 1)
})

test_that("species matrices are binary, deduplicated, and follow retained-ID order", {
  source_function("make_species_matrix.R")

  hex_occ <- data.frame(
    hex_id = c("h1", "h1", "h1", "h2", "h3"),
    species = c("sp_a", "sp_a", "sp_b", "sp_b", "sp_c")
  )

  mat <- make_species_matrix(hex_occ, valid_hex_ids = c("h2", "h1"))

  expect_true(is.matrix(mat))
  expect_identical(rownames(mat), c("h2", "h1"))
  expect_true(all(mat %in% c(0, 1)))
  expect_equal(mat["h1", "sp_a"], 1L)
  expect_equal(mat["h2", "sp_a"], 0L)
  expect_false("sp_c" %in% colnames(mat))
  expect_error(
    make_species_matrix(hex_occ, c("h1", "h1")),
    "unique"
  )
  expect_error(
    make_species_matrix(hex_occ, c("h1", "missing")),
    "absent"
  )
})

test_that("matrix bundle keeps Chao2 metadata in species-matrix row order", {
  source_function("add_incidence_unit.R")
  source_function("calculate_miner_detection.R")
  source_function("estimate_chao2_for_hex.R")
  source_function("calculate_chao2_table.R")
  source_function("make_species_matrix.R")
  source_function("build_matrix_bundle.R")

  make_hex_records <- function(hex_id, species_names) {
    data.frame(
      hex_id = hex_id,
      date = rep(
        as.Date(c("2020-01-01", "2020-02-01", "2020-03-01")),
        times = c(2, 2, 1)
      ),
      dataResourceName = "resource",
      samplingProtocol = "protocol",
      species = c(
        species_names[1], species_names[2],
        species_names[1], species_names[2],
        species_names[2]
      )
    )
  }

  hex_occ <- rbind(
    make_hex_records("hex_b", c("sp_c", "sp_d")),
    make_hex_records("hex_a", c("Manorina melanocephala", "sp_a"))
  )

  bundle <- build_matrix_bundle(
    hex_occ = hex_occ,
    incidence_unit = "month_year",
    min_sampling_units = 1L,
    min_species = 1L,
    min_completeness = 0,
    noisy_miner_species = "Manorina melanocephala",
    show_progress = FALSE
  )

  expect_identical(
    as.character(bundle$hex_metadata$hex_id),
    rownames(bundle$species_matrix)
  )
  expect_identical(rownames(bundle$species_matrix), c("hex_a", "hex_b"))
  expect_true(all(bundle$species_matrix %in% c(0, 1)))
  expect_equal(
    bundle$hex_metadata$miner_detection_rate[
      bundle$hex_metadata$hex_id == "hex_a"
    ],
    2 / 3
  )
})

test_that("Podani-Jaccard components have the expected toy decomposition", {
  source_function("calculate_lcbd_beta_matrices.R")

  species_matrix <- matrix(
    c(
      1, 0,
      0, 1,
      1, 1
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("h1", "h2", "h3"), c("sp_a", "sp_b"))
  )

  beta <- calculate_lcbd_beta_matrices(species_matrix)
  total <- beta$dist_mats$total
  replacement <- beta$dist_mats$replacement
  richness_difference <- beta$dist_mats$richness_difference

  expected_total <- matrix(
    c(0, 1, 0.5, 1, 0, 0.5, 0.5, 0.5, 0),
    nrow = 3,
    byrow = TRUE,
    dimnames = dimnames(total)
  )
  expected_replacement <- matrix(
    c(0, 1, 0, 1, 0, 0, 0, 0, 0),
    nrow = 3,
    byrow = TRUE,
    dimnames = dimnames(total)
  )

  expect_equal(total, expected_total)
  expect_equal(replacement, expected_replacement)
  expect_equal(richness_difference, expected_total - expected_replacement)
  expect_equal(total, t(total))
  expect_equal(total, replacement + richness_difference)
  expect_equal(beta$site_mean_replacement, c(h1 = 0.5, h2 = 0.5, h3 = 0))
  expect_equal(
    beta$site_mean_richness_difference,
    c(h1 = 0.25, h2 = 0.25, h3 = 0.5)
  )
})

test_that("Podani-Jaccard calculation rejects malformed community matrices", {
  source_function("calculate_lcbd_beta_matrices.R")

  valid <- matrix(
    c(1, 0, 0, 1),
    nrow = 2,
    dimnames = list(c("h1", "h2"), c("sp_a", "sp_b"))
  )

  non_binary <- valid
  non_binary[1, 1] <- 2
  expect_error(calculate_lcbd_beta_matrices(non_binary), "presence-absence")

  empty_site <- valid
  empty_site[1, ] <- 0
  expect_error(calculate_lcbd_beta_matrices(empty_site), "at least one species")

  duplicated_sites <- valid
  rownames(duplicated_sites) <- c("h1", "h1")
  expect_error(calculate_lcbd_beta_matrices(duplicated_sites), "duplicated")
})
