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
prepare_metadata <- function(project_id, base_path, norm_method = NULL, log_file) {

  log_message("Preparing metadata with combining biomass measurements (fcm or qpcr) with sample data.", status = "start", log_file)

  # Define data directories
  input_folder <- file.path(base_path, "input_data")
  qiime2_folder <- file.path(base_path, "qiime2_output")

  # Locate and validate metadata
  metadata_file <- list.files(qiime2_folder, pattern = "metadata.*\\.(tsv|txt|csv)$", full.names = TRUE)[1]
  if (is.na(metadata_file)) {
    error_message <- glue::glue("Error: Metadata file missing in {qiime2_folder}")
    log_message(error_message, status = "error", log_file)
    stop(error_message, call. = FALSE)
  }

  # Read metadata flexibly
  ext = tolower(tools::file_ext(metadata_file))
  metadata <- if (ext == "csv") {
    readr::read_csv(metadata_file, show_col_types = FALSE)
  } else {
    readr::read_delim(metadata_file, delim = "\t", show_col_types = FALSE)
  }

  # standardize QIIME2 sample ID headers
  if ("#SampleID" %in% colnames(metadata)) {
    metadata <- dplyr::rename(metadata, SampleID = `#SampleID`)
  }

  # Clean European decimal notation and auto-convert numeric columns safely
  metadata <- metadata %>%
    dplyr::mutate(dplyr::across(where(is.character), ~ {
      cleaned <- stringr::str_replace_all(.x, ",", ".")
      ifelse(stringr::str_detect(cleaned, "^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$"), as.numeric(cleaned), .x)
    }))

  # Metadata column validation --------------------------------------------

  # Required metadata columns
  required_col <- c("SampleID", "sample_or_control", "DNA_Concentration", "na_type")
  missing_col <- required_col[!required_col %in% colnames(metadata)]

  if (length(missing_col) > 0) {
    error_message <- glue::glue("Error: Metadata is missing required columns: {paste(missing_col, collapse = `, `)}")
    log_message(error_message, status = "error", log_file)
    stop(error_message, call. = FALSE)
  }

  # adjust zero values in DNA concentration safely
  metadata <- metadata %>%
    dplyr::mutate(DNA_Concentration = dplyr::if_else(DNA_Concentration == 0, 0.0000001, DNA_Concentration))

  # Adding Biomass data to metadata ---------------------------------------

  # Flow Cytometry (FCM)
  if (!is.null(norm_method) && tolower(norm_method) == "fcm") {
    if (!"cells_per_ml" %in% colnames(metadata)) {

      fcm_file <- list.files(qiime2_folder, pattern = "fcm.*\\.(csv|tsv|txt)$", full.names = TRUE)

      if (!is.na(fcm_file)) {
        # Read fcm flexibly
        ext <- tolower(tools::file_ext(fcm_file))
        fcm <- if (ext == "csv") {
          readr::read_csv(fcm_file, show_col_types = FALSE)
        } else {
          readr::read_delim(fcm_file, delim = "\t", show_col_types = FALSE)
        }
        fcm <- dplyr::select(fcm, SampleID, cells_per_ml)
        metadata <- dplyr::left_join(metadata, fcm, by = "SampleID")
        log_message("FCM biomass data successfully merged into metadata.", status = "info", log_file)
      } else {
        log_message("Warning: FCM normalization requested, but no fcm file found.", status = "warning", log_file)
      }
    }
  }

  # Quantitative PCR (qPCR)
  if (!is.null(norm_method) && tolower(norm_method) == "qpcr") {
    if (!"copies_per_ml" %in% colnames(metadata)) {

      qpcr_file <- list.files(qiime2_folder, pattern = "qpcr.*\\.(csv|tsv|txt)$|qPCR.*\\.(csv|tsv|txt)$", full.names = TRUE)

      if (!is.na(qpcr_file)) {
        # Read fcm flexibly
        ext <- tolower(tools::file_ext(qpcr_file))
        qpcr <- if (ext == "csv") {
          readr::read_csv(qpcr_file, show_col_types = FALSE)
        } else {
          readr::read_delim(qpcr_file, delim = "\t", show_col_types = FALSE)
        }
        qpcr <- dplyr::select(qpcr, SampleID, copies_per_ml)
        metadata <- dplyr::left_join(metadata, qpcr, by = "SampleID")
        log_message("qPCR biomass data successfully merged into metadata.", status = "info", log_file)
      } else {
        log_message("Warning: qPCR normalization requested, but no qPCR file found.", status = "warning", log_file)
      }
    }
  }

  # Save the updated metadata
  output_file <- file.path(input_folder, glue::glue("{project_id}_prepared_metadata.tsv"))
  readr::write_delim(metadata, file = output_file, delim = "\t")
  log_message("Metadata preparation completed successfully.", status = "success", log_file)
}
