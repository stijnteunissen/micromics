#' Create a Phyloseq Object
#'
#' This function creates a `phyloseq` object using input data files such as
#' feature tables, taxonomic assignments, phylogenetic tree, and unified
#' metadata into a single `phyloseq` object for downstream analysis.
#'
#' @inheritParams create_folders
#'
#' @details
#' This function performs the following steps:
#' \itemize{
#'   \item Defines the paths to the required input files (feature table, rooted tree, taxonomy, and metadata).
#'   \item Searches for and retrieves these files from the `input_data` directory.
#'   \item Calls the `qza_to_phyloseq()` function from the `qiime2R` package to generate a `phyloseq` object based on the provided input files.
#'   \item Adds read count information to the sample metadata within the `phyloseq` object.
#' }
#'
#' The function assumes the following files are present in the `input_data` directory:
#' \itemize{
#'   \item `table.qza`: Feature table containing sample feature data.
#'   \item `rooted-tree.qza`: Phylogenetic tree.
#'   \item `classifier.qza`: Taxonomic classification file.
#'   \item `metadata_formatted.tsv`: Unified sample metadata.
#' }
#'
#' The resulting `phyloseq` object is essential for downstream analyses and integrates all input files into a single, structured object.
#' The created `phyloseq` object is saved as an RDS file named `<project_name>_phyloseq_uncleaned.rds` in the `output_data/rds_files/Before_cleaning_rds_files` directory.
#'
#' @return
#' A `phyloseq` object that integrates feature tables, taxonomy, phylogenetic trees, and metadata.
#' This object is also saved as an RDS file for further usage in downstream analyses.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' physeq_object <- creating_physeq_object(projects)
#' }
#'
#' @export
create_phyloseq = function(project_id, base_path, log_file) {

  log_message("Creating raw phyloseq object from QIIME2 artifacts and metadata", status = "start", log_file)

  # Define data and storage directories
  input_folder = file.path(base_path, "input_data")
  raw_rds_folder = file.path(base_path, "raw_rds")

  # Search for the required artifat and metadata files
  table_file <- list.files(input_folder, pattern = "table.*\\.qza$", full.names = TRUE, recursive = TRUE)
  rooted_tree_file <- list.files(input_folder, pattern = "rooted-tree.*\\.qza$", full.names = TRUE, recursive = TRUE)
  taxonomy_file <- list.files(input_folder, pattern = "classifier.*\\.qza", full.names = TRUE, recursive = TRUE)
  metadata_file <- list.files(input_folder, pattern = "metadata.*\\.tsv", full.names = TRUE)

  # Validate that critical components exist before attempting import
  if (length(table_file) == 0 || length(taxonomy_file) == 0 || length(metadata_file) == 0) {
    error_message <- "Error: Missing table, taxonomy, or metadata in input_data."
    log_message(error_message, status = "error", log_file)
    stop(error_message, call. = FALSE)
  }

  # Build the phyloseq object based on tree availability
  if (length(rooted_tree_file) == 0) {
    log_message("No phylogenetic tree found. Building phyloseq without tree.", status = "info", log_file)
    physeq <- qiime2R::qza_to_phyloseq(
      features = table_file,
      taxonomy = taxonomy_file,
      metadata = metadata_file
    )
  } else {
    log_message("Phylogenetic tree located. Building complete phyloseq object.", status = "info", log_file)
    physeq <- qiime2R::qza_to_phyloseq(
      features = table_file,
      tree = rooted_tree_file,
      taxonomy = taxonomy_file,
      metadata = metadata_file
    )
  }

  # Add total read count per sample to the sample_data (metadata)
  phyloseq::sample_data(physeq)$read_count <- phyloseq::sample_sums(physeq)

  # save the uncleaned phyloseq object as an RDS file
  output_file_path = file.path(raw_rds_folder, glue::glue("{project_id}_phyloseq_uncleaned.rds"))
  saveRDS(physeq, file = output_file_path)

  # Log completion and return the object to the R session
  log_message(glue::glue("Uncleaned phyloseq object stored successfully at: {output_file_path}"), status = "info", log_file)
  log_message("Phyloseq object generation complete.", status = "success", log_file)

  return(physeq)
}
