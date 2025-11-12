#' Convert Phyloseq Data to Tibble and Export by Taxonomic Level
#'
#' This function transforms a `phyloseq` object into a tibble for each specified taxonomic level
#' using the \code{psmelt} function from the \pkg{phyloseq} package.
#'
#' @inheritParams group_tax
#'
#' @details
#' The primary task of this function is to convert a `phyloseq` object into a tibble at each specified
#' taxonomic level using the \code{psmelt} function.
#'
#' @return
#' The function saves multiple `phyloseq`-derived tibbles as RDS files and returns a list of tibbles:
#'   \itemize{
#'     \item **If \code{norm_method} is \code{NULL}:**
#'           A tibble of copy number–corrected counts is saved as
#'           \code{<project_name>_psmelt_<tax>_level_copy_number_corrected_counts.rds} for each taxonomic level.
#'
#'     \item **If \code{norm_method = "fcm"}:**
#'           Two tibbles are saved for each taxonomic level:
#'           \itemize{
#'             \item A tibble of copy number–corrected counts.
#'             \item A tibble of FCM-normalized, rarefied counts saved as
#'                   \code{<project_name>_psmelt_<tax>_level_fcm_normalised_cell_concentration_rarefied.rds}.
#'           }
#'
#'     \item **If \code{norm_method = "qpcr"}:**
#'           Two tibbles are saved for each taxonomic level:
#'           \itemize{
#'             \item A tibble of copy number–corrected counts.
#'             \item A tibble of qPCR-normalized, rarefied counts saved as
#'                   \code{<project_name>_psmelt_<tax>_level_qpcr_normalised_cell_concentration_rarefied.rds}.
#'           }
#'   }
#'
#' @examples
#' \dontrun{
#'   # Export data without normalization (copy number–corrected counts only)
#'   result <- psdata_to_tibble(physeq = rarefied_genus_physeq)
#'
#'   # Export data with flow cytometry normalization
#'   result <- psdata_to_tibble(physeq = rarefied_genus_physeq, norm_method = "fcm")
#'
#'   # Export data with qPCR normalization
#'   result <- psdata_to_tibble(physeq = rarefied_genus_physeq, norm_method = "qpcr")
#' }
#'
#' @export
psdata_to_tibble = function(physeq = rarefied_genus_physeq,
                            norm_method = NULL,
                            taxrank = c("Phylum", "Class", "Order", "Family", "Genus")) {

  log_message(paste("Step 12: Convert psdata to tibble.", paste(projects, collapse = ", ")), log_file)

  project_name = projects
  project_folder = paste0(base_path, project_name)
  output_folder_rds_files = paste0(project_folder, "/output_data/rds_files/After_cleaning_rds_files/")
  if (!dir.exists(output_folder_rds_files)) {dir.create(output_folder_rds_files, recursive = TRUE)}

  results = list()

  for (tax in taxrank) {
    log_message(paste("Processing taxonomic level:", tax), log_file)
    tax_folder = file.path(output_folder_rds_files, paste0(tax))
    if(!dir.exists(tax_folder)){dir.create(tax_folder)}

    if (is.null(norm_method)) {
      psdata = physeq[[paste0("psdata_copy_number_corrected_", tax)]]
      psdata_tibble = psdata %>% psmelt() %>% as_tibble()

      output_file_path = paste0(tax_folder, "/", project_name, "_psmelt_", tax, "_level_copy_number_corrected_counts.rds")
      saveRDS(psdata_tibble, file = output_file_path)
      log_message(paste("psmelt data copy number corrected counts", tax, "level saved as .rds object in", output_file_path), log_file)

      results[[paste0("psmelt_copy_number_corrected_", tax)]] <- psdata_tibble
    } else if (norm_method == "fcm") {
      psdata = physeq[[paste0("psdata_copy_number_corrected_", tax)]]
      psdata_tibble = psdata %>% psmelt() %>% as_tibble()

      output_file_path = paste0(tax_folder, "/", project_name, "_psmelt_", tax, "_level_copy_number_corrected_counts.rds")
      saveRDS(psdata_tibble, file = output_file_path)
      log_message(paste("psmelt data copy number corrected counts", tax, "level saved as .rds object in", output_file_path), log_file)

      psdata_fcm = physeq[[paste0("psdata_fcm_norm_rarefied_", tax)]]
      psdata_tibble_fcm = psdata_fcm %>% psmelt() %>% as_tibble()

      output_file_path = paste0(tax_folder, "/", project_name, "_psmelt_", tax, "_level_fcm_normalised_cell_concentration_rarefied.rds")
      saveRDS(psdata_tibble_fcm, file = output_file_path)
      log_message(paste("psmelt data fcm normalised cell concentration (cells per ml/gram sample)", tax, "level rarefied saved as .rds object in", output_file_path), log_file)

      results[[paste0("psmelt_copy_number_corrected_", tax)]] <- psdata_tibble
      results[[paste0("psmelt_fcm_norm_rarefied_", tax)]] <- psdata_tibble_fcm
    } else if (norm_method == "qpcr") {
      psdata = physeq[[paste0("psdata_copy_number_corrected_", tax)]]
      psdata_tibble = psdata %>% psmelt() %>% as_tibble()

      output_file_path = paste0(tax_folder, "/", project_name, "_psmelt_", tax, "_level_copy_number_corrected_counts.rds")
      saveRDS(psdata_tibble, file = output_file_path)
      log_message(paste("psmelt data copy number corrected counts", tax, "level saved as .rds object in", output_file_path), log_file)

      psdata_qpcr = physeq[[paste0("psdata_qpcr_norm_rarefied_", tax)]]
      psdata_tibble_qpcr = psdata_qpcr %>% psmelt() %>% as_tibble()

      output_file_path = paste0(tax_folder, "/", project_name, "_psmelt_", tax, "_level_qpcr_normalised_celL_concentration_rarefied.rds")
      saveRDS(psdata_tibble_qpcr, file = output_file_path)
      log_message(paste("psmelt qpcr normalised cell concentration (cells per ml/gram sample)", tax, "level rarefied saved as .rds object in", output_file_path), log_file)

      results[[paste0("psmelt_copy_number_corrected_", tax)]] <- psdata_tibble
      results[[paste0("psmelt_qpcr_norm_rarefied_", tax)]] <- psdata_tibble_qpcr
    }
  }
  return(results)

  log_message("Successfully converted to a tibble.", log_file)
}
