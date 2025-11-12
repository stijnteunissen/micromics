#' Unify and Format Metadata
#'
#' This function merges and formats metadata from various sources, including
#' QIIME metadata, experimental sample metadata, and qPCR or FCM data, to create a unified
#' metadata file for downstream analyses.
#'
#' @inheritParams create_folders
#'
#' @details
#' The function ensures that the metadata is unified and correctly formatted for further analysis by:
#' \itemize{
#'   \item Reading and processing the `metadata_extra.tsv` file, which must contain at least:
#'     \itemize{
#'       \item `SampleID`: A unique identifier for each sample.
#'       \item `sample_type`: Indicates whether the sample is a `sample`, `mock`, or `blank`.
#'       \item `DNA_Concentration`: The DNA concentration (in ng/µl).
#'     }
#'   \item Optionally processing and integrating qPCR or FCM data:
#'     \itemize{
#'       \item If qPCR data is available, calculating the mean for duplicates and merging it with `metadata_extra`.
#'       \item If FCM data is available, calculating the mean for duplicates and merging it with `metadata_extra`.
#'     }
#'   \item Reading the `metadata.tsv` file (QIIME metadata), ensuring it contains the `SampleID` column, and combining it with the processed `metadata_extra`.
#'   \item Writing the final combined metadata to a file. Note that the output file is named by concatenating the project name with `_metadata_formatted.tsv` and is saved in the project's `input_data` folder.
#' }
#'
#' All metadata files (QIIME metadata, experimental sample metadata, and qPCR/FCM data) must include
#' the `SampleID` column for proper merging. This column serves as the key to align data from multiple sources.
#'
#' @note
#' This function requires that the folder structure has been set up (using the `create_folders` function) before running.
#'
#' @return A data frame containing the unified metadata. The data frame is also saved as a file in the project's `input_data` folder.
#'
#' @examples
#' \dontrun{
#' # Process and unify metadata for a project
#' unified_metadata <- unify_metadata(projects)
#' }
#'
#' @export
unify_metadata <- function(projects) {

  log_message(paste("Step 2: Merge metadata: adding biomass (fcm or qpcr) data to the original metadata.", paste(projects, collapse = ", ")), log_file)

  # Define the project folder and destination folder paths
  project_name <- projects
  project_folder <- paste0(base_path, project_name)
  destination_folder <- paste0(project_folder, "/input_data/")

  # Read additional metadata
  metadata_extra_file <- list.files(destination_folder, pattern = "metadata_extra\\.(tsv|txt|csv)$", full.names = TRUE)[1]
  if (is.na(metadata_extra_file)) {
    message = paste("No metadata file (.tsv/.txt/.csv) found in ", destination_folder)
    stop(message)
    log_message(message, log_file)
  }

  ext = tolower(tools::file_ext(metadata_extra_file))
  metadata_extra <- if (ext == "csv") {
    read_csv(metadata_extra_file, col_names = TRUE, show_col_types = FALSE)
  } else {
    read_delim(metadata_extra_file, col_names = TRUE, show_col_types = FALSE)
  }

  # Initialize flags and data holders
  data_available <- FALSE
  qPCR_combined <- NULL
  FCM_combined <- NULL

  # Process multiple qPCR files
  qPCR_files <- list.files(destination_folder, pattern = "qPCR.*\\.(csv|tsv|txt)$", full.names = TRUE)
  if (length(qPCR_files) > 0) {
    qPCR_data <- qPCR_files %>%
      lapply(function(f) {
        ext <- tolower(file_ext(f))
        if (ext == "csv") {
          read_csv(f, col_names = TRUE, show_col_types = FALSE)
        } else {
          read_delim(f, col_names = TRUE, show_col_types = FALSE)
        }}) %>%
      bind_rows()

    # Check if the column 'sq_calc_mean' exists in the qPCR data
    if ("sq_calc_mean" %in% colnames(qPCR_data)) {
      # If it exists, use the existing sq_calc_mean values (taking the first occurrence per SampleID)
      qPCR_combined <- qPCR_data
    } else {
      # Otherwise, calculate the mean of SQ values and name the result 'sq_mean'
      qPCR_combined <- qPCR_data %>%
        group_by(SampleID) %>%
        summarise(sq_mean = mean(SQ, na.rm = TRUE), .groups = "drop")
    }

    # Merge the qPCR results with the extra metadata
    metadata_extra <- metadata_extra %>%
      left_join(qPCR_combined, by = "SampleID")

    data_available <- TRUE
  } else {
    log_message(paste("No qPCR data file found for project:", project_name), log_file)
  }

  # Process multiple FCM files
  FCM_files <- list.files(destination_folder, pattern = "fcm.*\\.(csv|tsv|txt)$", full.names = TRUE)

  if (length(FCM_files) > 0) {
    FCM_combined <- FCM_files %>%
      lapply(function(f) {
        ext <- tolower(file_ext(f))
        if (ext == "csv") {
          read_csv(f, col_names = TRUE, show_col_types = FALSE)
        } else {
          read_delim(f, col_names = TRUE, show_col_types = FALSE)
        }
      }) %>%
      bind_rows() %>%
      group_by(SampleID) %>%
      summarise(cells_per_ml = mean(cells_per_ml, na.rm = TRUE), .groups = "drop")

    metadata_extra <- metadata_extra %>%
      left_join(FCM_combined, by = "SampleID")

    data_available <- TRUE
  } else {
    log_message(paste("No FCM data file found for project:", project_name), log_file)
  }

  # If no data was available, log a warning and leave metadata unchanged
  if (!data_available) {
    warning_message <- paste("Warning: No qPCR or FCM data available for project:", project_name, "\n")
    log_message(warning_message, log_file)
  } else {
    # Write the updated metadata_extra back to a file
    # write_delim(metadata_extra, file = paste0(destination_folder, "/metadata_extra_updated.tsv"), delim = "\t", col_names = TRUE)
    message <- paste("Message: Metadata updated with available qPCR and/or FCM data for project:", project_name)
    log_message(message, log_file)
  }

  # Process main metadata files
  metadata_files <- list.files(destination_folder, pattern = "metadata\\.(tsv|txt|csv)$", full.names = TRUE)

  if (length(metadata_files) == 0) {
    error_message <- paste("Error: metadata file (.tsv/.txt/.csv) does not exist in", destination_folder)
    log_message(error_message, log_file)
    stop(error_message)
  }

  # Combine metadata files
  metadata_list <- lapply(metadata_files, function(f) {
    ext <- tolower(tools::file_ext(f))
    if (ext == "csv") {
      read_csv(f, col_names = TRUE, show_col_types = FALSE)
    } else {
      read_delim(f, col_names = TRUE, show_col_types = FALSE)
    }
  })

  metadata = bind_rows(metadata_list)


  # if (length(metadata_files) > 0) {
  #   metadata <- metadata_files %>%
  #     lapply(function(f) {
  #       ext <- tolower(file_ext(f))
  #       if (ext == "csv") {
  #         read_csv(f, col_names = TRUE, show_col_types = FALSE)
  #       } else {
  #         read_delim(f, col_names = TRUE, show_col_types = FALSE)
  #       }})
  # } else {
  #   error_message <- paste("Error: metadata file (.tsv/.txt/.csv) does not exist in", destination_folder)
  #   log_message(error_message, log_file)
  #   stop(error_message)
  # }

  # Rename '#SampleID' to 'SampleID' in metadata if it exists
  if ("#SampleID" %in% colnames(metadata)) {
    colnames(metadata)[colnames(metadata) == "#SampleID"] <- "SampleID"
  }

  # Ensure both files have a 'SampleID' column
  if (!("SampleID" %in% colnames(metadata))) {
    error_message <- paste0("Error: 'SampleID' column is missing in the metadata file.")
    log_message(error_message, log_file)
    stop(error_message)
  }
  if (!("SampleID" %in% colnames(metadata_extra))) {
    error_message <- paste0("Error: 'SampleID' column is missing in the metadata_extra file.")
    log_message(error_message, log_file)
    stop(error_message)
  }

  # Combine the metadata using an inner join on SampleID
  combined_metadata <- inner_join(metadata, metadata_extra, by = "SampleID")

  # if (any(grepl("\\s", colnames(combined_metadata)))) {
  #   error_message = paste0("Error: Spaces detected in metadata column names.")
  #   log_message(error_message, log_file)
  #   stop(error_message)
  # }
  #
  # bad_cells = which(apply(combined_metadata, 2, function(col)
  #   any(grepl("\\s", as.character(col)))))
  # if (length(bad_cells) > 0) {
  #   error_message = paste0("Error: Spaces detected in metadata values.")
  #   log_message(error_message, log_file)
  #   stop(error_message)
  # }

  for (col in names(combined_metadata)) {
    if (is.character(combined_metadata[[col]])) {
      has_comma = grepl(",", combined_metadata[[col]])
      if (any(has_comma)) {
        combined_metadata[[col]] = gsub(",", ".", combined_metadata[[col]])
      }
      if (all(grepl("^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$", combined_metadata[[col]]))) {
        combined_metadata[[col]] <- as.numeric(combined_metadata[[col]])
      }
    }
  }

  if (!"DNA_Concentration" %in% colnames(combined_metadata)) {
    error_message = paste0("Error: DNA_Concentration column not found in metadata.")
    log_message(error_message, log_file)
    stop(error_message)
  }

  zero = which(combined_metadata$DNA_Concentration == 0, arr.ind = TRUE)
  if (length(zero) > 0) {
    combined_metadata$DNA_Concentration[zero] = 0.000001
  }

  # Write the formatted metadata to a file
  write.table(combined_metadata, file = paste0(destination_folder, project_name, "_metadata_final.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

  return(combined_metadata)

  log_message("Metadata successfully merged.", log_file)
}
