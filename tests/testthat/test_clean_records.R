test_that("protocol classification handles structured, opportunistic, and unclear records", {
  source_function("classify_protocol.R")

  observed <- classify_protocol(
    data_resource = c(
      "BirdLife Australia, Birdata",
      "iNaturalist Australia",
      "Other resource",
      NA_character_
    ),
    sampling_protocol = c("2ha search", NA, "unclear survey", NA)
  )

  expect_identical(
    observed,
    c(
      "community_complete_or_structured",
      "opportunistic",
      "semi_structured_or_unclear",
      "semi_structured_or_unclear"
    )
  )
})

test_that("occurrence cleaning applies each filtering rule and records attrition", {
  source_function("classify_protocol.R")
  source_function("clean_occurrence_records.R")

  raw <- data.frame(
    record_id = c(
      "keep_threshold", "drop_absent", "drop_date", "drop_taxon",
      "drop_basis", "drop_protocol", "keep_missing_structured",
      "drop_missing_unclear", "drop_spatial", "duplicate_first",
      "duplicate_second", "tas_wrong", "tas_right", "mainland_wrong",
      "mainland_right", "drop_untrusted_missing", "keep_trusted_missing"
    ),
    occurrenceStatus = c(
      "PRESENT", "ABSENT", rep("PRESENT", 15)
    ),
    date = c(
      "2020-06-01", "2020-06-01", "2019-12-31", rep("2020-07-01", 14)
    ),
    species = c(
      "sp_keep", "sp_absent", "sp_old", NA, "sp_basis", "sp_protocol",
      "sp_structured", "sp_unclear", "sp_spatial", "sp_duplicate",
      "sp_duplicate", "Tas bird", "Tas bird", "Mainland bird",
      "Mainland bird", "sp_untrusted", "sp_trusted"
    ),
    longitude = c(1:9, 20, 20, 21:26),
    latitude = c(
      rep(-35, 11), -35, -42, -42, -35, -35, -35
    ),
    basisOfRecord = c(
      rep("HUMAN_OBSERVATION", 4), "PRESERVED_SPECIMEN",
      rep("HUMAN_OBSERVATION", 12)
    ),
    samplingProtocol = c(
      "2ha search", "survey", "survey", "survey", "survey",
      "excluded survey", "Area search", "casual survey", "survey",
      "survey", "survey", "survey", "survey", "survey", "survey",
      NA, NA
    ),
    dataResourceName = c(rep("Other", 15), "Other", "Trusted"),
    coordinateUncertaintyInMeters = c(
      2500, rep(100, 5), NA, NA, 100, 100, 100, 100, 100, 100, 100,
      100, 100
    ),
    stringsAsFactors = FALSE
  )

  raw_path <- tempfile(fileext = ".csv")
  on.exit(unlink(raw_path), add = TRUE)
  readr::write_csv(raw, raw_path)

  fake_coordinate_cleaner <- function(x, lon, lat, species, tests, value) {
    x$.sea <- x$record_id != "drop_spatial"
    x$.zer <- TRUE
    x$.equ <- TRUE
    x
  }

  rules <- list(
    trusted_resources = "Trusted",
    protocol_exclusions = "excluded survey",
    tasmanian_endemics = "Tas bird",
    mainland_obligates = "Mainland bird"
  )

  cleaned <- clean_occurrence_records(
    raw_path = raw_path,
    start_date = as.Date("2020-01-01"),
    end_date = as.Date("2020-12-31"),
    coord_precision = 2500,
    cleaning_rules = rules,
    protocol_filter = "all_retained",
    show_progress = FALSE,
    coordinate_cleaner = fake_coordinate_cleaner
  )

  expect_setequal(
    cleaned$records$record_id,
    c(
      "keep_threshold", "keep_missing_structured", "duplicate_first",
      "tas_right", "mainland_right", "keep_trusted_missing"
    )
  )
  expect_false("duplicate_second" %in% cleaned$records$record_id)
  expect_identical(
    cleaned$attrition$step,
    c(
      "raw_input",
      "present_and_date_window",
      "required_taxon_date_coordinates",
      "basis_of_record_filter",
      "protocol_screen_and_classification",
      "coordinate_cleaner_flags_added",
      "coordinate_uncertainty_and_spatial_validity",
      "deduplicated_lat_lon_date_species",
      "tasmanian_endemic_and_mainland_obligate_range_filter"
    )
  )
  expect_true(all(diff(cleaned$attrition$n_records) <= 0))
  expect_equal(tail(cleaned$attrition$n_records, 1), nrow(cleaned$records))
  expect_false(any(c(".sea", ".zer", ".equ") %in% names(cleaned$records)))

  structured <- clean_occurrence_records(
    raw_path = raw_path,
    start_date = as.Date("2020-01-01"),
    end_date = as.Date("2020-12-31"),
    coord_precision = 2500,
    cleaning_rules = rules,
    protocol_filter = "community_complete_or_structured",
    show_progress = FALSE,
    coordinate_cleaner = fake_coordinate_cleaner
  )
  expect_setequal(
    structured$records$record_id,
    c("keep_threshold", "keep_missing_structured")
  )
})

test_that("occurrence cleaning rejects invalid settings before coordinate cleaning", {
  source_function("classify_protocol.R")
  source_function("clean_occurrence_records.R")

  rules <- list(
    trusted_resources = character(),
    protocol_exclusions = "excluded",
    tasmanian_endemics = character(),
    mainland_obligates = character()
  )

  expect_error(
    clean_occurrence_records(
      raw_path = tempfile(fileext = ".csv"),
      start_date = as.Date("2020-01-01"),
      end_date = as.Date("2020-12-31"),
      coord_precision = 2500,
      cleaning_rules = rules,
      protocol_filter = "not_a_filter",
      show_progress = FALSE
    ),
    "Unknown protocol_filter"
  )
  expect_error(
    clean_occurrence_records(
      raw_path = tempfile(fileext = ".csv"),
      start_date = as.Date("2021-01-01"),
      end_date = as.Date("2020-12-31"),
      coord_precision = 2500,
      cleaning_rules = rules,
      show_progress = FALSE
    ),
    "ordered Date"
  )
  expect_error(
    clean_occurrence_records(
      raw_path = tempfile(fileext = ".csv"),
      start_date = as.Date("2020-01-01"),
      end_date = as.Date("2020-12-31"),
      coord_precision = -1,
      cleaning_rules = rules,
      show_progress = FALSE
    ),
    "non-negative"
  )
})
