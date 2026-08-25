#' Normalize Phyloseq Data
#'
#' This function applies normalization to a `phyloseq` object, converting
#' data to absolute values based on 16S rRNA copy numbers and sample biomass,
#' using either flow cytometry (FCM) data or qPCR data. The function can apply
#' copy number correction prior to biomass normalization.
#'
#' @inheritParams remove_mock
#'
#' @param norm_method A character string specifying the normalization method. Options are:
#'   \itemize{
#'     \item `"fcm"`: Normalize based on flow cytometry data, converting abundances to cell concentrations (cells/mL or per gram sample).
#'     \item `"qpcr"`: Normalize based on qPCR data, converting abundances to cell equivalents (cells/mL or per gram sample).
#'     \item `NULL`: Apply only copy number correction without further normalization.
#'   }
#'
#' @param copy_correction A logical value indicating whether the data should be corrected for
#'   the predicted 16S rRNA copy numbers prior to biomass normalization. Options are:
#'   \itemize{
#'     \item `TRUE`: Both relative and absolute abundances are corrected using the predicted copy numbers.
#'     \item `FALSE`: Abundances are not corrected by copy number. Note that qPCR normalization requires
#'           copy number correction to provide absolute data.
#'   }
#'
#' @details
#' The function follows these steps based on the chosen parameters:
#'
#' 1. **Copy Number Correction (if `copy_correction = TRUE`):**
#'    - Correct ASV abundances by dividing each count by its predicted 16S rRNA copy number.
#'      The prediction is based on the method described in
#'      ["Accounting for 16S rRNA copy number prediction uncertainty and its implications in bacterial diversity analyses"](https://dx.doi.org/10.1038/s43705-023-00266-0).
#'    - This correction adjusts for variability in 16S rRNA gene copy numbers across taxa, enabling
#'      the calculation of cell equivalents.
#'
#' 2. **FCM Normalization (`norm_method = "fcm"`):**
#'    - When `copy_correction = TRUE`: The copy number–corrected abundances are multiplied by the FCM data,
#'      where FCM data (cells per mL or per gram) is included in the metadata or a file with "fcm" in the name,
#'      with the column `cells_per_ml`.
#'    - When `copy_correction = FALSE`: The raw abundances are multiplied by the FCM data without prior copy number correction.
#'
#' 3. **qPCR Normalization (`norm_method = "qpcr"`):**
#'    - The qPCR data, provided in 16S copies per mL or per gram sample (included in the metadata or
#'      a file with "qpcr" in the name and column `sq_calc_mean`), is used together with copy number
#'      predictions to calculate absolute abundances.
#'
#' ### DNA vs. RNA Normalization
#' The interpretation of normalized data depends on the nucleic acid type:
#'   - **DNA:** Normalized abundances usually represent **cells per mL (or per gram)**, assuming one genome copy per cell.
#'   - **RNA:** Normalized abundances often represent **copies per cell equivalent per mL (or per gram)**;
#'      RNA reflects transcriptional activity and may vary considerably with cell condition.
#'
#' @references
#' Gao, Y., & Wu, M. (2023). Accounting for 16S rRNA copy number prediction uncertainty
#' and its implications in bacterial diversity analyses. *ISME Communications, 3*(1), 59.
#' doi:[10.1038/s43705-023-00266-0](https://dx.doi.org/10.1038/s43705-023-00266-0)
#'
#' @return
#' The function saves multiple `phyloseq` objects as RDS files:
#'   \itemize{
#'     \item `<project_name>_phyloseq_asv_level_without_copy_number_corrected_counts.rds`: if `copy_correction = FALSE`;
#'           a phyloseq object with uncorrected counts.
#'     \item `<project_name>_phyloseq_asv_level_copy_number_corrected_counts.rds`: if `copy_correction = TRUE`;
#'           a phyloseq object with counts corrected by predicted 16S copy numbers.
#'     \item `<project_name>_phyloseq_asv_level_fcm_normalised_cell_concentration.rds`: if `norm_method = "fcm"` and `copy_correction = TRUE`;
#'           a phyloseq object with abundances normalized based on FCM data and copy number correction.
#'     \item `<project_name>_phyloseq_asv_level_fcm_normalised_cell_concentration_without_copy_number_corrected_count.rds`: if `norm_method = "fcm"` and `copy_correction = FALSE`;
#'           a phyloseq object with FCM normalization applied without copy number correction.
#'     \item `<project_name>_phyloseq_asv_level_qpcr_normalised_cell_concentration.rds`: if `norm_method = "qpcr"`;
#'           a phyloseq object with abundances normalized to cell equivalents using qPCR data.
#'   }
#'
#' The relative phyloseq object (without biomass normalization) is saved in the
#' `output_data/rds_files/After_cleaning_rds_files/ASV` directory, and all biomass-normalized
#' phyloseq objects are saved in the `output_data/rds_files/Before_cleaning_rds_files` directory.
#'
#' @examples
#' \dontrun{
#' # Apply only copy number correction
#' result <- normalize_data(physeq = physeq, norm_method = NULL)
#'
#' # Normalize using flow cytometry (FCM) data
#' result <- normalize_data(physeq = physeq, norm_method = "fcm", copy_correction = TRUE)
#'
#' # Normalize using qPCR data
#' result <- normalize_data(physeq = physeq, norm_method = "qpcr", copy_correction = TURE)
#' }
#'
#' @export
normalise_data = function(physeq, norm_method = NULL, copy_correction = TRUE, project_id, base_path, log_file) {

  log_message("Starting copy number correction and biomass normalisation", status = "start", log_file)

  # Define output directory paths
  input_folder <- file.path(base_path, "input_data")
  raw_rds_folder <- file.path(base_path, "raw_rds")
  figure_folder <- file.path(base_path, "figures")

  # Initialize placeholders for the final outputs
  physeq_copy_corrected <- NULL
  physeq_biomass_normalised <- NULL

  # Execute rRNA gene copy number correction via RasperGade16S profiles
  if (copy_correction == TRUE) {
    log_message("Executing 16S rRNA gene copy number correction via RasperGade profiles.", status = "info", log_file)

    # Locate and read RasperGade16S RDS file
    raspergade_file <- list.files(input_folder, pattern = "prediction.*\\.RDS", full.names = TRUE)
    if (length(raspergade_file) == 0) {
      error_message <- "Error: RaserGade output file wiht 16s copies predictions file missing from input_data."
      log_message(error_message, status = "error", log_file)
      stop(error_message, call. = FALSE)
    }

    raspergade_rds <- readRDS(raspergade_file)

    # Clean raspergade output to a standaridised dataframe
    raspergade_df <- raspergade_rds$discrete %>%
      dplyr::rename(
        OTU = label,
        copy_number = x,
        probability = probs
      ) %>%
      dplyr::select(OTU, copy_number, probability) %>%
      dplyr::mutate(OTU = as.character(OTU))

    # Extract otu table and convert to a tibble
    otu_df <- data.frame(otu_table(physeq))
    otu_df$OTU <- rownames(otu_df)
    otu_tibble <- as_tibble(otu_df) %>%
      dplyr::select(OTU, everything())

    # Merge abundance data with copy number values
    joined_otu_tibble <- otu_tibble %>%
      dplyr::inner_join(raspergade_df, by = "OTU")

    # Apply copy number correction
    corrected_otu_tibble <- joined_otu_tibble %>%
      dplyr::mutate(across(
        .cols = c(everything(), -OTU, -copy_number, -probability),
        .fns = ~ ceiling(. / copy_number)
      )) %>%
      dplyr::select(-copy_number, -probability)

    # Convert data back to a valid phyloseq object
    corrected_matrix <- as.data.frame(corrected_otu_tibble[, -1])
    rownames(corrected_matrix) <- corrected_otu_tibble$OTU

    otu_corrected <- phyloseq::otu_table(corrected_matrix, taxa_are_rows = TRUE)
    phyloseq::taxa_names(physeq_copy_corrected) <- corrected_otu_tibble$OTU

    physeq_copy_corrected <- physeq
    phyloseq::otu_table(physeq_copy_corrected) <- otu_corrected
  }

  # FCM Normalisation
  if (!is.null(norm_method) && norm_method == "fcm") {

    log_message("Flow Cytometry (FCM) biomass normalisation.", status = "info", log_file)

    # Use copy number corrected data if copy_correction is TRUE, otherwise use the original phyloseq
    fcm_input <- if (copy_correction) physeq_copy_corrected else physeq

    # Convert otu table to a dataframe
    otu_df_fcm <- data.frame(phyloseq::otu_table(fcm_input))
    otu_df_fcm$OTU <- rownames(otu_df_fcm)

    # Pivot the OTU matrix into a long format
    fcm_long <- otu_df_fcm %>%
      tidyr::pivot_longer(
        cols = -OTU,
        names_to = "SampleID",
        values_to = "Abundance"
        )

    # Extract sample metadata
    sample_df_fcm <- data.frame(phyloseq::sample_data(fcm_input))
    sample_df_fcm$SampleID <- rownames(sample_df_fcm)

    # Combine abundance data with sample metadata and calculate scaling factors per individual sample
    joined_tibble_fcm <- fcm_long %>%
      dplyr::inner_join(sample_df_fcm, by = "SampleID") %>%
      dplyr::group_by(SampleID) %>%
      dplyr::mutate(
        max_cells_per_ml = max(cells_per_ml, na.rm = TRUE),
        scale_factor = ifelse(max_cells_per_ml > 1e7, 10 ^ ceiling(log10(max_cells_per_ml / 1e7)), 1),
        scaled_cells_per_ml = cells_per_ml / scale_factor)

    # Calculate relative abundance and normalise wiht biomass data
    joined_tibble_fcm_norm <- joined_tibble_fcm %>%
      dplyr::group_by(SampleID) %>%
      dplyr::mutate(relative_abund = Abundance / sum(Abundance)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(norm_abund = ceiling(relative_abund * scaled_cells_per_ml)) %>%
      dplyr::select(OTU, norm_abund, SampleID)

    # Reshape the normalised data back into a wide matrix format
    fcm_norm_wide <- joined_tibble_fcm_norm %>%
      tidyr::pivot_wider(
        names_from = SampleID,
        values_from = norm_abund
        )

    # Reconstruct the phyloseq OTU table from the wide formatted data
    fcm_norm_matrix <- as.data.frame(fcm_norm_wide[, -1])
    rownames(fcm_norm_matrix) <- fcm_norm_wide$OTU

    otu_fcm_norm <- phyloseq::otu_table(fcm_norm_matrix, taxa_are_rows = TRUE)
    phyloseq::taxa_names(otu_fcm_norm) <- fcm_norm_wide$OTU

    # Generate normalised phyloseq object
    physeq_biomass_normalised <- fcm_input
    phyloseq::otu_table(physeq_biomass_normalised) <- otu_fcm_norm

    # Add the calculated scaling factors to the sample metadata
    updated_sample_df <- as.data.frame(phyloseq::sample_data(physeq_biomass_normalised))
    original_row_names <- rownames(updated_sample_df)
    updated_sample_df$SampleID <- original_row_names

    # Extract unique scaling parameters per sample identifier
    scaling_factors_df <- joined_tibble_fcm %>%
      dplyr::select(SampleID, scaled_cells_per_ml, scale_factor) %>%
      dplyr::distinct()

    # Merge scaling metrics back into the metadata
    final_sample_df <- data.frame(updated_sample_df) %>%
      dplyr::left_join(scaling_factors_df, by = "SampleID") %>%
      dplyr::select(-SampleID)

    rownames(final_sample_df) <- original_row_names
    phyloseq::sample_data(physeq_biomass_normalised) <- phyloseq::sample_data(final_sample_df)

    # qPCR Normalisation
  } else if (!is.null(norm_method) && norm_method == "qpcr" && copy_correction == TRUE) {

    log_message("qPCR biomass normalisation.", status = "info", log_file)

    # Convert otu table to a dataframe
    otu_df_qpcr <- data.frame(otu_table(physeq))
    otu_df_qpcr$OTU <- rownames(otu_df_qpcr)

    # Pivot the OTU matrix into a long format
    qpcr_long <- otu_df_qpcr %>%
      tidyr::pivot_longer(
        cols = -OTU,
        names_to = "SampleID",
        values_to = "Abundance")

    # Extract sample metadata
    sample_df_qpcr <- data.frame(sample_data(physeq))
    sample_df_qpcr$SampleID <- rownames(sample_df_qpcr)

    # Combine sample metadata with copy number
    joined_tibble_qpcr <- qpcr_long %>%
      dplyr::inner_join(sample_df_qpcr, by = "SampleID") %>%
      dplyr::inner_join(raspergade_df, by = "OTU")

    # Calculate scaling factor based on the maximum copies per ml
    joined_tibble_qpcr_scaled <- joined_tibble_qpcr %>%
      dplyr::group_by(SampleID) %>%
      dplyr::mutate(
        max_copies_per_ml = max(copies_per_ml, na.rm = TRUE),
        scale_factor = ifelse(max_copies_per_ml > 1e7, 10 ^ ceiling(log10(max_copies_per_ml / 1e7)), 1),
        scaled_copies_per_ml = copies_per_ml / scale_factor
      )

    # Calculate normalised abundance
    joined_tibble_qpcr_norm <- joined_tibble_qpcr_scaled %>%
      dplyr::group_by(SampleID) %>%
      dplyr::mutate(
        relative_abundance = Abundance / sum(Abundance),
        absolute_abundance_qpcr = relative_abundance * scaled_copies_per_ml,
        norm_abund = ceiling(absolute_abundance_qpcr / copy_number)
      ) %>%
      dplyr::select(OTU, norm_abund, SampleID)

    # Reshape the normalised data back into a wide matrix format
    qpcr_norm_wide <- joined_tibble_qpcr_norm %>%
      tidyr::pivot_wider(
        names_from = SampleID,
        values_from = norm_abund
      )

    # Reconstruct the phyloseq OTU table from the wide fromatted data
    qpcr_norm_matrix <- as.data.frame(qpcr_norm_wide[, -1])
    rownames(qpcr_norm_matrix) <- qpcr_norm_wide$OTU

    otu_qpcr_norm <- phyloseq::otu_table(qpcr_norm_matrix, taxa_are_rows = TRUE)
    phyloseq::taxa_names(otu_qpcr_norm) <- qpcr_norm_wide$OTU

    # Generate the normalised phyloseq object
    physeq_biomass_normalised <- qpcr_input
    phyloseq::otu_table(physeq_biomass_normalised) <- otu_qpcr_norm

    # Add the calculated scaling factors back to the final sample metadata
    updated_sample_df <- as.data.frame(phyloseq::sample_data(physeq_biomass_normalised))
    original_row_names <- rownames(updated_sample_df)
    updated_sample_df$SampleID <- original_row_names

    # Extract unique scaling parameters per sample identifier
    scaling_factors_df <- joined_tibble_qpcr_scaled %>%
      dplyr::select(SampleID, scaled_copies_per_ml, scale_factor) %>%
      dplyr::distinct()

    # Merge scaling metrics back into the metadata and strictly preserve row names
    final_sample_df <- data.frame(updated_sample_df) %>%
      dplyr::left_join(scaling_factors_df, by = "SampleID") %>%
      dplyr::select(-SampleID)

    rownames(final_sample_df) <- original_row_names
    phyloseq::sample_data(physeq_biomass_normalised) <- phyloseq::sample_data(final_sample_df)

    # Handle configuration errors safely
  } else if (!is.null(norm_method) && norm_method == "qpcr" && copy_correction == FALSE) {
    error_message <- "Error: norm_method is set to 'qpcr' but copy_correction is FALSE. The qPCR normalisation method requires copy number correction to be enabled."
    log_message(error_message, status = "error", log_file)
    stop(error_message, call. = FALSE)
  }

  # Generate validation plots if copy number correction is enabled
  if (copy_correction == TRUE) {

    log_message("Raspergade copy prediction compared to rrnDB database and plot generation.", status = "info", log_file)

    # Path to the database
    zip_file_path <- file.path(input_folder, "rrnDB-5.9_pantaxa_stats_NCBI.tsv.zip")

    # Define the official secure URL for the rrnDB pan-taxa archive
    rrndb_url <- "https://rrndb.umms.med.umich.edu/downloads/rrnDB-5.9_pantaxa_stats_NCBI.tsv.zip"

    # Automatically download the database archive if it does not exist locally
    if (!file.exists(zip_file_path)) {
      log_message("rrnDB archive not found locally. downloading zip file.", status = "info", log_file)

      # Execute secure transfer using binary mode to protect zip structure
      download.file(
        url = rrndb_url,
        destfile = zip_file_path,
        mode = "wb",
        quiet = TRUE
      )
    }

    # Verify that the required database archive exists and is valid
    if (!file.exists(zip_file_path) || file.info(zip_file_path)$size == 0) {
      error_message <- "Error: Downloaded rrnDB ZIP file is missing or corrupted."
      log_message(error_message, status = "error", log_file)
      base::stop(error_message, call. = FALSE)
    }

    # Extract the downloaded database files to the input directory
    unzip(zip_file_path, exdir = input_folder)
    rrndb_tsv_file <- list.files(input_folder, pattern = "pantaxa_stats_NCBI\\.tsv$", full.names = TRUE)
    rrndb_database_tsv <- readr::read_tsv(rrndb_tsv_file, show_col_types = FALSE)

    # Filter for genus level and match naming
    rrndb_database <- rrndb_database_tsv %>%
      dplyr::filter(rank == "genus") %>%
      dplyr::select(Genus = name, dplyr::everything())

    # Extract taxonomy information directly from the original phyloseq object
    tax_tibble <- as.data.frame(phyloseq::tax_table(physeq)) %>%
      tibble::rownames_to_column("OTU") %>%
      dplyr::select(OTU, Genus) %>%
      dplyr::as_tibble()

    # Join the taxonomic data, RasperGade predictions, and reference database stats
    validation_table <- tax_tibble %>%
      dplyr::left_join(raspergade_df, by = "OTU") %>%
      dplyr::left_join(rrndb_database, by = "Genus")

    # Filter and categorize the comparison profiles based on probability bounds
    validation_plot_data <- validation_table %>%
      dplyr::select(OTU, Genus, mean, copy_number, probability) %>%
      dplyr::filter(!base::is.na(mean) & !base::is.na(copy_number)) %>%
      dplyr::mutate(probability_rate = dplyr::case_when(
        probability > 0.9 ~ "High (> 0.9)",
        probability >= 0.5 & probability <= 0.9 ~ "Medium (>= 0.5 & <= 0.9)",
        probability < 0.5 ~ "Low (< 0.5)"
      ))

    # Perform statistical verification via Pearson correlation coefficient
    cor_test <- stats::cor.test(validation_plot_data$mean, validation_plot_data$copy_number, method = "pearson")
    r_val <- cor_test$estimate
    label_text <- base::sprintf("r = %.2f", r_val)

    # Build the diagnostic comparison visualization
    copy_number_comparison <- ggplot(validation_plot_data, aes(x = mean, y = copy_number, colour = probability_rate)) +
      geom_point(size = 1) +
      geom_smooth(method = "lm", se = FALSE, colour = "black", aes(group = 1)) +
      geom_abline(intercept = 0, slope = 1, colour = "red") +
      scale_colour_manual(values = c("High (> 0.9)" = "darkgreen", "Medium (>= 0.5 & <= 0.9)" = "orange", "Low (< 0.5)" = "red4")) +
      xlim(0, 15) +
      ylim(0, 15) +
      theme_bw() +
      labs(
        title = "ASV copy number comparison",
        x = "Mean copy number reference rrnDB",
        y = "Predicted copy number per genus",
        colour = "Probability rate"
      ) +
      theme(
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = c(0.05, 0.95),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "white", colour = "grey80")
      ) +
      annotate("text", x = 2.25, y = 11, label = label_text, hjust = 1, vjust = 1, size = 5)

    # Export visualisation
    ggplot2::ggsave(filename = file.path(figure_folder, glue::glue("{project_id}_copy_number_comparison.png")), plot = copy_number_comparison, width = 8, height = 8)
    ggplot2::ggsave(filename = file.path(figure_folder, glue::glue("{project_id}_copy_number_comparison.pdf")), plot = copy_number_comparison, width = 8, height = 8)
  }

  if (copy_correction == FALSE) {
    physeq_copy_corrected = physeq
  }

  if (!is.null(norm_method)) {
    if (copy_correction == TRUE) {
      saveRDS(file = file.path(raw_rds_folder, glue::glue("{project_id}_phyloseq_copy_corrected_counts.rds")), object = physeq_copy_corrected)
      saveRDS(file = file.path(raw_rds_folder, glue::glue("{project_id}_phyloseq_biomass_normalised_counts.rds")), object = physeq_biomass_normalised)

      return(list(physeq_copy_corrected = physeq_copy_corrected, physeq_biomass_normalised = physeq_biomass_normalised))
    } else {
      saveRDS(file = file.path(raw_rds_folder, glue::glue("{project_id}_phyloseq_without_copy_corrected_counts.rds")), object = physeq_copy_corrected)
      saveRDS(file = file.path(raw_rds_folder, glue::glue("{project_id}_phyloseq_biomass_normalised_counts.rds")), object = physeq_biomass_normalised)

      return(list(physeq_copy_corrected = physeq_copy_corrected, physeq_biomass_normalised = physeq_biomass_normalised))
    }
  } else if (copy_correction == TRUE) {
    saveRDS(file = file.path(raw_rds_folder, glue::glue("{project_id}_phyloseq_copy_corrected_counts.rds")), object = physeq_copy_corrected)

    return(list(physeq_copy_corrected = physeq_copy_corrected))

  } else {
    saveRDS(file = file.path(raw_rds_folder, glue::glue("{project_id}_phyloseq_without_copy_corrected_counts.rds")), object = physeq_copy_corrected)

    return(list(physeq_copy_corrected = physeq_copy_corrected))
  }
}
