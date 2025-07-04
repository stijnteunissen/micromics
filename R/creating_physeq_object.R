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
creating_physeq_object = function(projects) {

  log_message(paste("Step 3: Creating phyloseq object (physeq): physeq is created using the table, rooted tree, classifier and metadata.", paste(projects, collapse = ", ")), log_file)

  project_name = projects
  project_folder = paste0(base_path, project_name)
  destination_folder = paste0(project_folder, "/input_data")
  output_folder_rds_files = paste0(project_folder, "/output_data/rds_files/Before_cleaning_rds_files/")

  # Search for the required files
  table_file <- list.files(destination_folder, pattern = "table.*\\.qza$", full.names = TRUE, recursive = TRUE)
  rooted_tree_file <- list.files(destination_folder, pattern = "rooted-tree.*\\.qza$", full.names = TRUE, recursive = TRUE)
  taxonomy_file <- list.files(destination_folder, pattern = "classifier.*\\.qza", full.names = TRUE, recursive = TRUE)
  metadata_file <- list.files(destination_folder, pattern = "metadata_final\\.tsv", full.names = TRUE)

  # Create the phyloseq object
  physeq <- qza_to_phyloseq(
    features = table_file,
    tree = rooted_tree_file,
    taxonomy = taxonomy_file,
    metadata = metadata_file
  )

  # Add read counts to metadata
  phyloseq::sample_data(physeq)$read_count = phyloseq::sample_sums(physeq)

  # save uncleaned psdata
  output_file_path = paste0(output_folder_rds_files, project_name, "_phyloseq_uncleaned.rds")
  saveRDS(physeq, file = output_file_path)
  log_message(paste("Uncleaned phyloseq object saved as .rds object in", output_file_path), log_file)

  return(physeq)

  log_message("Phyloseq successfully created.", log_file)

}
