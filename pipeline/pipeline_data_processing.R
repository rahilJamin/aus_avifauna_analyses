# =============================================================================
# pipeline_data_processing.R
# Targets for external inputs, occurrence cleaning, hex assignment, iNEXT/Chao2
# filtering, and aligned model-data preparation.
# =============================================================================

# Stage settings ----
# Keeping small configuration objects prevents an unrelated setting change from
# invalidating the whole pipeline. For example, changing a GLLVM seed should not
# make occurrence cleaning outdated.
data_input_paths <- cfg[c(
  "raw_ala_path",
  "hex_2_5_path",
  "selected_metrics_2_5_path",
  "trait_path",
  "ibra_path"
)]

data_output_paths <- cfg[c(
  "table_dir",
  "cleaned_records_path",
  "hex_occ_path",
  "chao_table_path",
  "hex_metadata_path",
  "species_matrix_path",
  "candidate_model_df_path",
  "model_df_path",
  "beta_matrices_path"
)]

cleaning_settings <- cfg[c(
  "analysis_id",
  "start_date",
  "end_date",
  "protocol_filter",
  "coord_precision_m"
)]

hex_settings <- cfg[c("input_lonlat_crs")]

matrix_settings <- cfg[c(
  "analysis_id",
  "incidence_unit",
  "min_sampling_units",
  "min_completeness",
  "min_species",
  "noisy_miner_species",
  "show_progress"
)]

model_data_settings <- cfg[c(
  "analysis_id",
  "noisy_miner_species",
  "urban_cover_max",
  "model_coord_crs"
)]

# Target definitions ----

data_processing_targets <- list(
# Required external inputs ----

  tar_target(
    raw_ala_file,
    data_input_paths$raw_ala_path,
    format = "file"
  ),
  tar_target(
    hex_grid_files,
    track_shapefile_files(data_input_paths$hex_2_5_path),
    format = "file"
  ),
  tar_target(
    landscape_metrics_file,
    data_input_paths$selected_metrics_2_5_path,
    format = "file"
  ),
  tar_target(
    species_traits_file,
    data_input_paths$trait_path,
    format = "file"
  ),
  tar_target(
    ibra_files,
    track_shapefile_files(data_input_paths$ibra_path),
    format = "file"
  ),


  # 01 Clean occurrence records ----

  tar_target(
    cleaned_records_bundle,
    clean_occurrence_records(
      raw_path = raw_ala_file,
      start_date = cleaning_settings$start_date,
      end_date = cleaning_settings$end_date,
      coord_precision = cleaning_settings$coord_precision_m,
      cleaning_rules = cleaning_rules,
      protocol_filter = cleaning_settings$protocol_filter
    )
  ),
  tar_target(
    cleaning_files,
    c(
      cleaned_records = write_csv_quiet(
        cleaned_records_bundle$records,
        data_output_paths$cleaned_records_path
      ),
      attrition = write_csv_quiet(
        dplyr::mutate(
          cleaned_records_bundle$attrition,
          analysis_id = cleaning_settings$analysis_id
        ),
        file.path(data_output_paths$table_dir, "01_attrition_clean_records.csv")
      )
    ),
    format = "file"
  ),


  # 02 Assign records to hexagons ----

  tar_target(
    hex_occurrences,
    assign_records_to_hex(
      cleaned_records = cleaned_records_bundle$records,
      hex_path = main_shapefile_path(hex_grid_files),
      input_lonlat_crs = hex_settings$input_lonlat_crs
    )
  ),
  tar_target(
    hex_assignment_files,
    {
      attrition <- tibble::tibble(
        step = c("cleaned_records", "assigned_to_2_5km_hex"),
        n_records = c(
          nrow(cleaned_records_bundle$records),
          nrow(hex_occurrences)
        ),
        n_hexagons = c(
          NA_integer_,
          dplyr::n_distinct(hex_occurrences$hex_id)
        ),
        n_species = c(
          dplyr::n_distinct(cleaned_records_bundle$records$species),
          dplyr::n_distinct(hex_occurrences$species)
        )
      )

      c(
        hex_occurrences = write_csv_quiet(
          hex_occurrences,
          data_output_paths$hex_occ_path
        ),
        attrition = write_csv_quiet(
          attrition,
          file.path(data_output_paths$table_dir, "02_attrition_hex_assignment.csv")
        )
      )
    },
    format = "file"
  ),


  # 03 Build the iNEXT species matrix ----

  tar_target(
    inext_matrix_bundle,
    build_matrix_bundle(
      hex_occ = hex_occurrences,
      incidence_unit = matrix_settings$incidence_unit,
      min_sampling_units = matrix_settings$min_sampling_units,
      min_species = matrix_settings$min_species,
      min_completeness = matrix_settings$min_completeness,
      noisy_miner_species = matrix_settings$noisy_miner_species,
      show_progress = matrix_settings$show_progress
    )
  ),
  tar_target(
    inext_matrix_files,
    c(
      chao_table = write_csv_quiet(
        inext_matrix_bundle$chao,
        data_output_paths$chao_table_path
      ),
      hex_metadata = write_csv_quiet(
        inext_matrix_bundle$hex_metadata,
        data_output_paths$hex_metadata_path
      ),
      species_matrix = save_rds_quiet(
        inext_matrix_bundle$species_matrix,
        data_output_paths$species_matrix_path
      ),
      attrition = write_csv_quiet(
        dplyr::mutate(
          inext_matrix_bundle$attrition,
          analysis_id = matrix_settings$analysis_id
        ),
        file.path(data_output_paths$table_dir, "03_attrition_inext_matrix.csv")
      )
    ),
    format = "file"
  ),


  # 03b Prepare aligned model data ----

  tar_target(
    landscape_data,
    read_rds_checked(
      landscape_metrics_file,
      "processed landscape metrics",
      required_class = "data.frame"
    )
  ),
  tar_target(
    model_data_bundle,
    prepare_lcbd_model_data(
      species_matrix = inext_matrix_bundle$species_matrix,
      hex_metadata = inext_matrix_bundle$hex_metadata,
      landscape_data = landscape_data,
      hex_path = main_shapefile_path(hex_grid_files),
      noisy_miner_species = model_data_settings$noisy_miner_species,
      urban_cover_max = model_data_settings$urban_cover_max,
      model_coord_crs = model_data_settings$model_coord_crs
    )
  ),
  tar_target(
    model_data_files,
    c(
      candidate_model_data = save_rds_quiet(
        model_data_bundle$candidate_model_df,
        data_output_paths$candidate_model_df_path
      ),
      model_data = save_rds_quiet(
        model_data_bundle$model_df,
        data_output_paths$model_df_path
      ),
      beta_matrices = save_rds_quiet(
        model_data_bundle$beta_matrices,
        data_output_paths$beta_matrices_path
      ),
      attrition = write_csv_quiet(
        dplyr::mutate(
          model_data_bundle$attrition,
          analysis_id = model_data_settings$analysis_id
        ),
        file.path(data_output_paths$table_dir, "03b_attrition_model_data.csv")
      )
    ),
    format = "file"
  )
)
