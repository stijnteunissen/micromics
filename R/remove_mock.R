#' Remove Mock Features from a Phyloseq Object
#'
#' This function removes mock features from a `phyloseq` object.
#' These mock features, which can appear in other samples due to cross-contamination,
#' are removed to minimize their impact on the analysis samples.
#' In addition, the function filters the dataset to retain only samples without controls.
#' Users can choose whether to remove the mock features by setting the `mock` parameter.
#'
#' @inheritParams decontam
#'
#' @param mock_genera A vector of genera representing the taxa that make up the mock community.
#'                    These taxa are used to identify mock features to be removed from the `phyloseq` object.
#'
#' @param mock A logical value determining whether to filter out mock features.
#'   \itemize{
#'     \item `TRUE`: Remove mock features from the `phyloseq` object and retain only samples for downstream analysis.
#'     \item `FALSE`: Retain mock features. Use this option if no mock community is present in the dataset.
#'   }
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item If `mock = FALSE`, the function filters the dataset to retain only samples without controls, leaving the mock features intact.
#'         This option is suitable for datasets where no mock community is included.
#'   \item If `mock = TRUE`, the function:
#'     \itemize{
#'       \item Identifies mock features based on the provided `mock_genera`.
#'       \item Removes the mock features from the dataset.
#'       \item Retains only samples without controls.
#'     }
#' }
#'
#' @return
#' A filtered `phyloseq` object is returned and saved as an RDS file named
#' `project_name_phyloseq_asv_level_without_mock.rds` in the `output_data/rds_files/Before_cleaning_rds_files/` directory.
#'
#' @examples
#' \dontrun{
#' # Remove mock ASVs from the phyloseq object
#' physeq_no_mock <- remove_mock(physeq = physeq, mock_genera = c("Mock_Genus1", "Mock_Genus2"), mock = TRUE)
#'
#' # Retain mock ASVs but filter to only samples (or when no mock community is present)
#' physeq_no_filter <- remove_mock(physeq = physeq, mock_genera = c("Mock_Genus1", "Mock_Genus2"), mock = FALSE)
#' }
#'
#' @export
remove_mock = function(physeq, mock_genera, mock = TRUE, project_id, base_path, log_file) {

  log_message("Removing mock-specific ASVs", status = "start", log_file)

  # Define storage directory path
  raw_rds_folder <- file.path(base_path, "raw_rds")

  # Validate presence of required sample data
  if (!("sample_or_control" %in% colnames(phyloseq::sample_data(physeq)))) {
    error_message <- paste0("Error: 'sample_or_control' column is missing from the sample data.")
    log_message(error_message, status = "error", log_file)
    stop(error_message, call. = FALSE)
  }

  # Conditional proceessing based on mock presence
  if (mock == FALSE) {
    log_message("The 'mock' parameter is set to FALSE. Retaining true samples only, mock ASVs will not be filtered.", status = "info", log_file)

    # Filter dataset to isolate true biological samples and drop zero abundance ASVs
    physeq_filtered <- physeq %>%
      phyloseq::subset_samples(sample_or_control == "sample") %>%
      phyloseq::prune_taxa(taxa_sums(.) > 0, .)

  } else {
    log_message("Isolating mock community control profiles to identify cross-contaminating mock ASVs.", status = "info", log_file)

    # Isolate mock reference profiles and drop zero abundace ASVs
    physeq_mock <- physeq %>%
      phyloseq::subset_samples(sample_or_control == "mock") %>%
      phyloseq::prune_taxa(taxa_sums(.) > 0, .)

    # Extract specific ASVs belonging to the mock genera
    mock_ASVs <- physeq_mock %>%
      phyloseq::subset_taxa(Genus %in% mock_genera) %>%
      phyloseq::taxa_names()

    log_message(glue::glue("Identified {length(mock_asvs)} unique ASVs matching the specified mock community genera."), status = "info", log_file)

    # Prune mock control samples and completely strip out their corresponding ASVs from samples
    physeq_filtered <- physeq %>%
      phyloseq::subset_samples(sample_or_control == "sample") %>%
      phyloseq::prune_taxa(!taxa_names(.) %in% mock_ASVs, .)
  }

  # Save results
  output_rds_path <- file.path(raw_rds_folder, glue::glue("{project_id}_phyloseq_asv_without_mock.rds"))
  saveRDS(physeq_filtered, file = output_rds_path)

  # Log
  log_message(glue::glue("Filtered phyloseq object stored successfully at: {output_rds_path}"), status = "info", log_file)
  log_message("Mock sample and mock-specific ASV removal complete.", status = "success", log_file)

  return(physeq_filtered)
}
