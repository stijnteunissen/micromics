#' Aggregate Phyloseq Data by Taxonomic Level
#'
#' This function aggregates ASV data at specified taxonomic levels (e.g., Phylum, Class, Order, Family, or Genus)
#' using the \code{tax_glom} function from the \pkg{phyloseq} package.
#'
#' @inheritParams rarefying
#'
#' @param taxrank A character vector indicating the taxonomic levels at which to group the data.
#'
#' @details
#' The function applies the \code{tax_glom} function to group ASVs at each specified taxonomic level. It creates
#' a dedicated folder for each taxonomic level under the output directory and saves the aggregated data as RDS files.
#'
#' @return
#' The function saves multiple `phyloseq` objects as RDS files.
#' The aggregated objects are saved in the output directory `output_data/rds_files/After_cleaning_rds_files/`.
#'
#' @examples
#' \dontrun{
#' # Aggregate data using flow cytometry normalization
#' result <- group_tax(physeq = rarefied_asv_physeq, norm_method = "fcm")
#'
#' # Aggregate data using qPCR normalization
#' result <- group_tax(physeq = rarefied_asv_physeq, norm_method = "qpcr")
#' }
#'
#' @export
group_tax = function(physeq = rarefied_asv_physeq,
                     norm_method = NULL,
                     taxrank = c("Phylum", "Class", "Order", "Family", "Genus"),
                     copy_correction = TRUE) {

  log_message(paste("Step 11: Tax glom: ASVs are merged at taxonomic ranks.", paste(projects, collapse = ", ")), log_file)

  project_name = projects
  project_folder = paste0(base_path, project_name)
  output_folder_rds_files = paste0(project_folder, "/output_data/rds_files/After_cleaning_rds_files/")
  if(!dir.exists(output_folder_rds_files)){dir.create(output_folder_rds_files)}

  if (copy_correction == TRUE) {
    cc = ""
  } else if (copy_correction == FALSE) {
    cc = "without_"
  }

  results = list()

  for (tax in taxrank) {
    log_message(paste("Processing taxonomic level:", tax), log_file)
    tax_folder = file.path(output_folder_rds_files, paste0(tax))
    if(!dir.exists(tax_folder)){dir.create(tax_folder)}

    if (is.null(norm_method)) {
      psdata = physeq[["psdata_asv_copy_number_corrected"]]
      psdata_tax = tax_glom(psdata, taxrank = tax)

      output_file_path = paste0(tax_folder, "/", project_name, "_phyloseq_", tax, "_level_", cc, "copy_number_corrected_counts.rds")
      saveRDS(psdata_tax, file = output_file_path)
      log_message(paste("phyloseq data copy number corrected counts", tax, "level saved as .rds object in", output_file_path), log_file)

      results[[paste0("psdata_copy_number_corrected_", tax)]] <- psdata_tax
    } else if (!is.null(norm_method) && norm_method == "fcm") {
      psdata = physeq[["psdata_asv_copy_number_corrected"]]
      psdata_tax = tax_glom(psdata, taxrank = tax)

      output_file_path = paste0(tax_folder, "/", project_name, "_phyloseq_", tax, "_level_copy_number_corrected_counts.rds")
      saveRDS(psdata_tax, file = output_file_path)
      log_message(paste("phyloseq data copy number corrected counts", tax, "level saved as .rds object in", output_file_path), log_file)

      psdata_fcm = physeq[["psdata_asv_fcm_norm_rarefied"]]
      psdata_tax_fcm = tax_glom(psdata_fcm, taxrank = tax)

      output_file_path = paste0(tax_folder, "/", project_name, "_phyloseq_", tax, "_level_fcm_normalised_cell_concentration_rarefied.rds")
      saveRDS(psdata_tax_fcm, file = output_file_path)
      log_message(paste("Phyloseq data fcm normalised cell concentration (cells per ml/gram sample)", tax, " level rarefied saved as .rds object in", output_file_path), log_file)

      results[[paste0("psdata_copy_number_corrected_", tax)]] <- psdata_tax
      results[[paste0("psdata_fcm_norm_rarefied_", tax)]] <- psdata_tax_fcm
    } else if (!is.null(norm_method) && norm_method == "qpcr") {
      psdata = physeq[["psdata_asv_copy_number_corrected"]]
      psdata_tax = tax_glom(psdata, taxrank = tax)

      output_file_path = paste0(tax_folder, "/", project_name, "_phyloseq_", tax, "_level_copy_number_corrected_counts.rds")
      saveRDS(psdata_tax, file = output_file_path)
      log_message(paste("phyloseq data copy_number corrected counts", tax, " level saved as .rds object in", output_file_path), log_file)

      psdata_qpcr = physeq[["psdata_asv_qpcr_norm_rarefied"]]
      psdata_tax_qpcr = tax_glom(psdata_qpcr, taxrank = tax)

      output_file_path = paste0(tax_folder, "/", project_name, "_phyloseq_", tax, "_level_qpcr_normalised_cell_concentration_rarefied.rds")
      saveRDS(psdata_tax_qpcr, file = output_file_path)
      log_message(paste("phyloseq qpcr normalised cell concentration (cells per ml/gram sample)", tax, " level rarefied saved as .rds object in", output_file_path), log_file)

      results[[paste0("psdata_copy_number_corrected_", tax)]] <- psdata_tax
      results[[paste0("psdata_qpcr_norm_rarefied_", tax)]] <- psdata_tax_qpcr
    }
  }

  return(results)

  log_message("Tax glom successfully aplied", log_file)
}
