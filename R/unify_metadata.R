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
  # Define the project folder and destination folder paths
  project_name = projects
  project_folder = paste0(base_path, project_name)
  destination_folder = paste0(project_folder, "/input_data/")

  # Add metadata
  metadata_extra_file = list.files(destination_folder, pattern = "metadata_extra\\.tsv$", full.names = TRUE)
  metadata_extra = read_delim(metadata_extra_file, col_names = TRUE, show_col_types = FALSE)

  # Initialize flags and data holders
  data_available = FALSE
  qPCR_combined = NULL
  FCM_combined = NULL

  # Process multiple qPCR files
  qPCR_files = list.files(destination_folder, pattern = "qPCR.*\\.csv$", full.names = TRUE)
  if (length(qPCR_files) > 0) {
    qPCR_combined = qPCR_files %>%
      lapply(read_delim, col_names = TRUE, show_col_types = FALSE) %>%
      bind_rows() %>%
      group_by(SampleID) %>%
      summarise(sq_mean = mean(SQ, na.rm = TRUE), .groups = "drop")

    metadata_extra = metadata_extra %>%
      left_join(qPCR_combined, by = "SampleID")

    data_available = TRUE
  } else {
    log_message(paste("No qPCR data found for project:", project_name), log_file)
  }

  # Process multiple FCM files
  FCM_files = list.files(destination_folder, pattern = "fcm.*\\.csv$", full.names = TRUE)

  if (length(FCM_files) > 0) {
    FCM_combined = FCM_files %>%
      lapply(read_delim, col_names = TRUE, show_col_types = FALSE) %>%
      bind_rows() %>%
      group_by(SampleID) %>%
      summarise(cells_per_ml = mean(cells_per_ml, na.rm = TRUE), .groups = "drop")

    metadata_extra = metadata_extra %>%
      left_join(FCM_combined, by = "SampleID")

    data_available = TRUE
  } else {
    log_message(paste("No FCM data found for project:", project_name), log_file)
  }

  # If no data was available, print a warning and leave metadata unchanged
  if (!data_available) {
    warning_message = paste("Warning: No qPCR or FCM data available for project:", project_name, "\n")
    log_message(warning_message, log_file)
  } else {
    # If data was available, write the updated metadata back to a file
    write_delim(metadata_extra, file = paste0(destination_folder, "/metadata_extra_updated.tsv"), delim = "\t", col_names = TRUE)
    message = paste("Message: Metadata updated with available qPCR and/or FCM data for project:", project_name)
    log_message(message, log_file)
  }

  # Process metadata
  metadata_files = list.files(destination_folder, pattern = "metadata\\.tsv$", full.names = TRUE)

  # Combine metadata files
  if (length(metadata_files) > 0) {
    metadata = metadata_files %>%
      lapply(read_delim, col_names = TRUE, show_col_types = FALSE) %>%
      bind_rows()
  } else {
    error_message = paste("Error: metadata.tsv file does not exist.")
    log_message(error_message, log_file)
    stop(error_message)
  }

  # Rename '#SampleID' to 'SampleID' in metadata if it exists
  if ("#SampleID" %in% colnames(metadata)) {
    colnames(metadata)[colnames(metadata) == "#SampleID"] <- "SampleID"
  }

  # Ensure both files have 'SampleID' column
  if (!("SampleID" %in% colnames(metadata))) {
    error_message = paste0("Error: 'SampleID' column is missing in the metadata file.")
    log_message(error_message, log_file)
    stop(error_message)
  }
  if (!("SampleID" %in% colnames(metadata_extra))) {
    error_message = paste0("Error: 'SampleID' column is missing in the metadata_extra file.")
    log_message(error_message, log_file)
    stop(error_message)
  }

  # Combine the metadata using a left join on SampleID
  combined_metadata = inner_join(metadata, metadata_extra, by = "SampleID")

  # Write the formatted metadata to a file
  write.table(combined_metadata, file = paste0(destination_folder, project_name, "_metadata_formatted.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

  return(combined_metadata)
}
