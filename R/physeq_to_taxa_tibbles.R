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
physeq_to_taxa_tibbles = function(physeq, norm_method = NULL, copy_correction = TRUE, taxrank = c("Phylum", "Class", "Order", "Family", "Genus"), project_id, base_path, log_file) {

  log_message("Starting taxonomic agglomeration and conversion to tibble", status = "start", log_file)

  # Define output directory path
  clean_rds_folder <- file.path(base_path, "clean_rds")

  results = list()

  for (tax in taxrank) {
    log_message(glue::glue("Processing taxonomic level: {tax}"), status = "info", log_file)

    # Create a director for the specific taxonomic rank
    tax_folder = file.path(clean_rds_folder, tax)
    if (!dir.exists(tax_folder)) {
      dir.create(tax_folder, recursive = TRUE)
    }

    physeq_rmp = physeq[["physeq_rmp_rarefied"]]
    physeq_rmp_glom = phyloseq::tax_glom(physeq_rmp, taxrank = tax)

    psdata_rmp <- physeq_rmp_glom %>%
      phyloseq::psmelt() %>%
      tibble::as_tibble()

    if (copy_correction == TRUE) {
      saveRDS(file = file.path(tax_folder, glue::glue("{project_id}_psdata_{tax}_copy_corrected_counts_rarefied.rds")), object = psdata_rmp)
    } else {
      saveRDS(file = file.path(tax_folder, glue::glue("{project_id}_psdata_{tax}_without_copy_corrected_counts_rarefied.rds")), object = psdata_rmp)
    }

    results[[glue::glue("psdata_rmp_{tax}")]] <- psdata_rmp

    if (!is.null(norm_method)) {
      physeq_qmp = physeq[["physeq_qmp_rarefied"]]
      physeq_qmp_glom = phyloseq::tax_glom(physeq_qmp, taxrank = tax)

      psdata_qmp <- physeq_qmp_glom %>%
        phyloseq::psmelt() %>%
        tibble::as_tibble()

      saveRDS(file = file.path(tax_folder, glue::glue("{project_id}_psdata_{tax}_biomass_normalised_counts_rarefied.rds")), object = psdata_qmp)

      results[[glue::glue("psdata_qmp_{tax}")]] <- psdata_qmp

    }
  }
  return(results)
}
