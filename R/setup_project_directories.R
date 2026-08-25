#' Setup Project Directories
#'
#' This function creates a folder structure for the project, ensuring that all
#' necessary directories exist and that specific files required for downstream
#' analyses are present.
#'
#' @param project_id A character vector containing the name of the project.
#' @param base_path A character string indicating the base directory where the project folder is located. This path is used to locate the project folder and create the required subfolders.
#' @param log_file A character string specifying the path to the log file where warnings and errors will be recorded.
#'
#' @details
#' This function facilitates the setup of downstream analyses by:
#' \itemize{
#'   \item Creating a consistent directory structure for each project, including subfolders.
#'   \item Copying essential files from the `qiime2_output` folder into the `input_data` folder. These essential files include:
#'     \itemize{
#'       \item `table.qza`: The feature table from QIIME2.
#'       \item `rooted-tree.qza`: The phylogenetic tree used for diversity analysis.
#'       \item `classifier.qza`: The classifier used for taxonomy assignment.
#'       \item `metadata.tsv`: The QIIME2 sample metadata file required for analyses.
#'       \item `metadata_extra.tsv`: Any additional metadata provided for the samples.
#'       \item `pantaxa_stats_NCBI.tsv`: The reference database for copy number correction from [rrndb](https://rrndb.umms.med.umich.edu/downloads/).
#'       \item `prediction.RDS`: The predicted 16S copy numbers for each feature.
#'     }
#'   \item Checking for optional files that enhance analyses, such as:
#'     \itemize{
#'       \item `qPCR.csv`: Contains quantitative PCR data.
#'       \item `fcm.csv`: Contains flow cytometry data.
#'     }
#'   \item Logging warnings for missing optional files and errors for missing required files.
#' }
#'
#' The function ensures that downstream analyses—which rely on specific input
#' files (e.g., `table.qza`, `rooted-tree.qza`, etc.) have access to these files
#' in the correct directory structure. If any required files are missing from
#' the `qiime2_output` folder, the function stops execution and logs an error
#' message.
#'
#' @note
#' Each project folder must already exist and
#' must contain a subfolder named `qiime2_output`, which holds the outputs of
#' QIIME2 analysis and other necessary files. The function sets up the folder
#' structure for downstream analysis within this project folder.
#'
#' @return
#' None. This function is called for its side effects.
#'
#' @examples
#' \dontrun{
#' # Define the base path and log file location
#' project_id <- "project_name"
#' base_path <- "path/to/projects"
#' log_file <- "path/to/log_file.log"
#'
#' # Create folder structures for projects
#' create_folders(projects)
#' }
#'
#' @export
setup_project_directories <- function(project_id, base_path, log_file = log_file) {

  # Log the initialization process
  log_message(glue::glue("Initialzing project structure for: {project_id}"), status = "start", log_file)

  # Create the required directories if they don't exist
  if(!dir.exists(file.path(base_path, "input_data"))){dir.create(file.path(base_path, "input_data"))}
  if(!dir.exists(file.path(base_path, "csv_exports"))){dir.create(file.path(base_path, "csv_exports"))}
  if(!dir.exists(file.path(base_path, "raw_rds"))){dir.create(file.path(base_path, "raw_rds"))}
  if(!dir.exists(file.path(base_path, "clean_rds"))){dir.create(file.path(base_path, "clean_rds"))}
  if(!dir.exists(file.path(base_path, "messages"))){dir.create(file.path(base_path, "messages"))}
  if(!dir.exists(file.path(base_path, "figures"))){dir.create(file.path(base_path, "figures"))}

  # Define source and destination paths
  qiime2_folder <- file.path(base_path, "qiime2_output")
  input_folder <- file.path(base_path, "input_data")

  project_files <- list.files(qiime2_folder, full.names = TRUE)

  # Check for required files stop if missing
  required_files <- c("table.*\\.qza$", "classifier.*\\.qza", "metadata.*\\.(tsv|txt|csv)$")
  for (file_pattern in required_files) {
    if (!any(grepl(file_pattern, project_files))) {
      error_message <- glue::glue("Error: {file_pattern} does not exist in {qiime2_folder}\n")
      log_message(error_message, status = "error", log_file)
      stop(error_message, call. = FALSE)
    }
  }

  # Check for optional files only warning if missing
  optional_files <- c("rooted-tree.*\\.qza$", "qPCR.*\\.(tsv|txt|csv)$","fcm.*\\.(tsv|txt|csv)$", "prediction*\\.RDS$")
  for (file_pattern in optional_files) {
    if (!any(grepl(file_pattern, project_files))) {
      warning_message <- glue::glue("Warning: {file_pattern} does not exist in {qiime2_folder}\n")
      log_message(warning_message, status = "warning", log_file)
    }
  }

  # Copy the relevant files
  copy_patterns <- paste(
    "table.*\\.qza$", "rooted-tree.*\\.qza$", "classifier.*\\.qza$",
    "qPCR.*\\.(tsv|txt|csv)$", "fcm.*\\.(tsv|txt|csv)$", "prediction.*\\.RDS$",
    sep = "|")

  # Filter and copy the target files
  files_to_copy <- list.files(qiime2_folder, pattern = copy_patterns, full.names = TRUE)
  file.copy(files_to_copy, input_folder, overwrite = TRUE)

  # Log successful completion
  log_message("Folder structure successfully created.", status = "success", log_file)
}

