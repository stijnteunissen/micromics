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
remove_mock = function(physeq = decontam_physeq,
                       mock_genera = mock_genera,
                       mock = TRUE) {

  psdata = physeq
  project_name = projects

  project_folder = paste0(base_path, project_name)
  output_folder_rds_files = paste0(project_folder, "/output_data/rds_files/Before_cleaning_rds_files/")

  if (mock == FALSE) {
    physeq_filtered = psdata
    message = paste0("message: mock parameter is FALSE. mock ASVs not filtered")
    log_message(message, log_file)
  } else if (mock == TRUE) {

    # if (!is.null(extra_genera)) {
    #   mock_genera = c(mock_genera, extra_genera)
    # } else
    #   mock_genera = mock_genera

    # extract mock OTU mock_genera
    mock_ASVs =
      psdata %>%
      subset_taxa(Genus %in% mock_genera) %>%
      taxa_names()

    # Filter phyloseq to remove mock samples
    physeq_no_mock =
      psdata %>%
      subset_samples(., sample_or_control == "sample") %>%
      prune_taxa(taxa_sums(.) > 0, .)

    physeq_filtered =
      physeq_no_mock %>%
      prune_taxa(!taxa_names(.) %in% mock_ASVs, .)
  }

  output_file_path = paste0(output_folder_rds_files, project_name, "_phyloseq_asv_level_without_mock.rds")
  saveRDS(physeq_filtered, file = output_file_path)
  log_message(paste("Phyloseq object without mock saved as .rds object in", output_file_path), log_file)

  return(physeq_filtered)
}
