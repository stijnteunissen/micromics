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
normalise_data = function(physeq = without_mock_physeq,
                          norm_method = NULL,
                          copy_correction = TRUE) {

  log_message(paste("Step 9: Copy number correction: Correct the copy number for relative data and biomass normalisation for absolute data.", paste(projects, collapse = ", ")), log_file)

  psdata = physeq
  project_name = projects

  project_folder = paste0(base_path, project_name)
  destination_folder = paste0(project_folder, "/input_data/")
  figure_folder_pdf = paste0(project_folder, "/figures/PDF_figures/")
  if(!dir.exists(figure_folder_pdf)) { dir.create(figure_folder_pdf) }
  figure_folder_png = paste0(project_folder, "/figures/PNG_figures/")
  if(!dir.exists(figure_folder_png)) { dir.create(figure_folder_png) }
  output_folder_rds_files_before = paste0(project_folder, "/output_data/rds_files/Before_cleaning_rds_files/")
  output_folder_rds_files_after = paste0(project_folder, "/output_data/rds_files/After_cleaning_rds_files/")
  output_asv_rds_files = paste0(output_folder_rds_files_after, "ASV/")
  if(!dir.exists(output_asv_rds_files)){dir.create(output_asv_rds_files)}

  # raspergate
  df_psdata = data.frame(otu_table(psdata))
  df_psdata$OTU = rownames(df_psdata)
  pstibble =
    as_tibble(df_psdata) %>%
    select(OTU, everything())

  rasperGade16S_file = list.files(destination_folder, pattern = "prediction\\.RDS$", full.names = TRUE)
  rasperGade16S_rds = readRDS(rasperGade16S_file)

  raspergade_df = rasperGade16S_rds$discrete %>%
    dplyr::rename(OTU = label,
           copy_number = x,
           probability = probs) %>%
    select(OTU, copy_number, probability) %>%
    as_tibble()

  joined_pstibble =
    pstibble %>%
    inner_join(., raspergade_df, by = "OTU")

  corrected_joined_pstibble =
    joined_pstibble %>%
    rowwise() %>%
    mutate(across(c(everything(), -OTU, -copy_number), ~ . / copy_number)) %>%
    mutate(across(c(everything(), -OTU, -copy_number), ~ ceiling(.))) %>%
    select(-copy_number)

  otu_corrected = otu_table(data.frame(corrected_joined_pstibble[, -1]), taxa_are_rows = TRUE)
  taxa_names(otu_corrected) = corrected_joined_pstibble$OTU

  psdata_copy_number_corrected = psdata
  otu_table(psdata_copy_number_corrected) = otu_corrected

  # FCM Normalization
  if (!is.null(norm_method) && norm_method == "fcm") {

    # Use copy number corrected data if copy_correction is TRUE, otherwise use the original psdata
    fcm_input <- if (copy_correction) psdata_copy_number_corrected else psdata

    df_psdata_fcm <- data.frame(otu_table(fcm_input))
    df_psdata_fcm$OTU <- rownames(df_psdata_fcm)

    psdata_fcm_long <- df_psdata_fcm %>%
      pivot_longer(cols = -OTU,
                   names_to = "SampleID",
                   values_to = "Abundance")

    df_sample_data_fcm <- data.frame(sample_data(fcm_input))
    df_sample_data_fcm$SampleID <- rownames(df_sample_data_fcm)

    joined_pstibble_fcm <- psdata_fcm_long %>%
      inner_join(df_sample_data_fcm, by = "SampleID") %>%
      group_by(SampleID) %>%
      mutate(max_cells_per_ml = max(cells_per_ml, na.rm = TRUE),
             scale_factor = ifelse(max_cells_per_ml > 1e7, 10 ^ ceiling(log10(max_cells_per_ml / 1e7)), 1),
             cells_per_ml = cells_per_ml / scale_factor)

    joined_pstibble_fcm_norm <- joined_pstibble_fcm %>%
      group_by(SampleID) %>%
      mutate(relative_abund = Abundance / sum(Abundance)) %>%
      ungroup() %>%
      mutate(norm_abund = ceiling(relative_abund * cells_per_ml)) %>%
      select(OTU, norm_abund, SampleID)

    fcm_norm_wide <- joined_pstibble_fcm_norm %>%
      pivot_wider(names_from = SampleID,
                  values_from = norm_abund)

    otu_fcm_norm <- otu_table(data.frame(fcm_norm_wide[, -1]), taxa_are_rows = TRUE)
    taxa_names(otu_fcm_norm) <- fcm_norm_wide$OTU

    psdata_fcm_norm <- fcm_input
    otu_table(psdata_fcm_norm) <- otu_fcm_norm

    # add scale factor to sample data
    df_tmp = data.frame(sample_data(psdata_fcm_norm))
    df_tmp = df_tmp %>% mutate(SampleID = rownames(df_tmp))

    # get unique scale factor values per sample
    df_sf = joined_pstibble_fcm %>%
      select(SampleID, scale_factor) %>%
      distinct()

    # if `sq_calc_mean` already exist in sample data, rename it to `orgin_sq_clac_mean`
    if ("cells_per_ml" %in% colnames(df_tmp)) {
      colnames(df_tmp)[colnames(df_tmp) == "cells_per_ml"] <- "orgin_cells_per_ml"
    }

    # merge into sample data
    df_tmp2 = df_tmp %>%
      left_join(df_sf, by = "SampleID")

    # restore rownames
    rownames(df_tmp2) = rownames(sample_data(psdata_fcm_norm))
    sample_data(psdata_fcm_norm) = sample_data(df_tmp2)

    # qpcr normalisatie
  } else if (!is.null(norm_method) && norm_method == "qpcr" && copy_correction == TRUE) {
    df_psdata_qpcr = data.frame(otu_table(psdata))
    df_psdata_qpcr$OTU = rownames(df_psdata_qpcr)

    psdata_qpcr_long =
      df_psdata_qpcr %>%
      pivot_longer(cols = -OTU,
                   names_to = "SampleID",
                   values_to = "Abundance")

    df_sample_data_qpcr = data.frame(sample_data(psdata))
    df_sample_data_qpcr$SampleID = rownames(df_sample_data_qpcr)

    joined_pstibble_qpcr =
      psdata_qpcr_long %>%
      inner_join(., df_sample_data_qpcr, by = "SampleID")

    # Merge `copy_number` from raspergade16s data
    joined_pstibble_qpcr_copy_number =
      joined_pstibble_qpcr %>%
      inner_join(., raspergade_df %>% select(OTU, copy_number), by = "OTU")

    # cell equivlants abundance in plaats van relative
    # abdance is cell equivlants
    # aangeven wat de eenheiden zijn vor relatief als het relatief als is of absoluut en de juiste eenheid kan per gram of ml

    results_list = list()

    if ("dna" %in% joined_pstibble_qpcr_copy_number$na_type) {
      joined_pstibble_qpcr_norm_dna =
        joined_pstibble_qpcr_copy_number %>%
        filter(na_type == "dna") %>%
        group_by(SampleID)

      if (!("sq_calc_mean" %in% colnames(joined_pstibble_qpcr_norm_dna))) {
        joined_pstibble_qpcr_norm_dna =
          joined_pstibble_qpcr_norm_dna %>%
          mutate(a = (sq_mean / insert_volume) * dilution_factor,
                 b = a * elution_volume_ul,
                 sq_calc_mean = b / sample_volume_or_mass) %>%
          ungroup()
        results_list$dna = joined_pstibble_qpcr_norm_dna
      } else {
        results_list$dna = joined_pstibble_qpcr_norm_dna %>% ungroup()
        log_message(paste("sq_calc_mean already exists, skipping its calculation for dna samples."), log_file)
      }
    }

    if ("rna" %in% joined_pstibble_qpcr_copy_number$na_type) {
      joined_pstibble_qpcr_norm_rna =
        joined_pstibble_qpcr_copy_number %>%
        filter(na_type == "rna") %>%
        group_by(SampleID)

      if (!("sq_calc_mean" %in% colnames(joined_pstibble_qpcr_norm_rna))) {
      joined_pstibble_qpcr_norm_rna =
        joined_pstibble_qpcr_norm_rna %>%
        mutate(a = (sq_mean / insert_volume) * dilution_factor,
               b = a * Final_volume_of_cDNA,
               c = b * (RNA_volume_used_Dnase_treatment / Dnase_Volume_used_cDNA_synthesis),
               d = c * (elution_volume_ul / RNA_volume_used_Dnase_treatment),
               sq_calc_mean = d / sample_volume_or_mass) %>%
        ungroup()
      results_list$rna = joined_pstibble_qpcr_norm_rna
      } else {
        results_list$rna = joined_pstibble_qpcr_norm_rna %>% ungroup()
        log_message("sq_calc_mean already exists, skipping its calculation for RNA samples.", log_file)
      }
    }

    if (!is.null(results_list$dna) && !is.null(results_list$rna)) {
      joined_pstibble_combined2 =
        bind_rows(results_list$dna, results_list$rna) %>%
        mutate(max_sq_calc_mean = max(sq_calc_mean, na.rm = TRUE),
               scale_factor = ifelse(max_sq_calc_mean > 1e7, 10 ^ ceiling(log10(max_sq_calc_mean / 1e7)), 1),
               sq_calc_mean = sq_calc_mean / scale_factor)

      joined_pstibble_combined =
        joined_pstibble_combined2 %>%
        group_by(SampleID) %>%
        mutate(relative_abundance = Abundance / sum(Abundance),
               absolute_abundance_qpcr = relative_abundance * sq_calc_mean,
               norm_abund = ceiling(absolute_abundance_qpcr / copy_number)) %>%
        ungroup() %>%
        select(OTU, norm_abund, SampleID)

    } else if (!is.null(results_list$dna)) {
      joined_pstibble_combined2 =
        results_list$dna %>%
        mutate(max_sq_calc_mean = max(sq_calc_mean, na.rm = TRUE),
               scale_factor = ifelse(max_sq_calc_mean > 1e7, 10 ^ ceiling(log10(max_sq_calc_mean / 1e7)), 1),
               sq_calc_mean = sq_calc_mean / scale_factor)

      joined_pstibble_combined =
        joined_pstibble_combined2 %>%
        group_by(SampleID) %>%
        mutate(relative_abundance = Abundance / sum(Abundance),
               absolute_abundance_qpcr = relative_abundance * sq_calc_mean,
               norm_abund = ceiling(absolute_abundance_qpcr / copy_number)) %>%
        ungroup() %>%
        select(OTU, norm_abund, SampleID)

    } else if (!is.null(results_list$rna)) {
      joined_pstibble_combined2 =
        results_list$rna %>%
        mutate(max_sq_calc_mean = max(sq_calc_mean, na.rm = TRUE),
               scale_factor = ifelse(max_sq_calc_mean > 1e7, 10 ^ ceiling(log10(max_sq_calc_mean / 1e7)), 1),
               sq_calc_mean = sq_calc_mean / scale_factor)

      joined_pstibble_combined =
        joined_pstibble_combined2 %>%
        group_by(SampleID) %>%
        mutate(relative_abundance = Abundance / sum(Abundance),
               absolute_abundance_qpcr = relative_abundance * sq_calc_mean,
               norm_abund = ceiling(absolute_abundance_qpcr / copy_number)) %>%
        ungroup() %>%
        select(OTU, norm_abund, SampleID)
    }

    qpcr_norm_wide =
      joined_pstibble_combined %>%
      pivot_wider(names_from = SampleID,
                  values_from = norm_abund)

    otu_qpcr_norm = otu_table(data.frame(qpcr_norm_wide[, -1]), taxa_are_rows = TRUE)
    taxa_names(otu_qpcr_norm) = qpcr_norm_wide$OTU
    psdata_qpcr_norm = psdata
    otu_table(psdata_qpcr_norm) = otu_qpcr_norm

    # Add sq_calc_mean to sample_data
    df_tmp = data.frame(sample_data(psdata_qpcr_norm))
    df_tmp = df_tmp %>% mutate(SampleID = rownames(df_tmp))

    # Get unique values of SampleID and sq_calc_mean
    df_sq = joined_pstibble_combined2 %>%
      select(SampleID, sq_calc_mean, scale_factor) %>%
      distinct()

    # if `sq_calc_mean` already exist in sample data, rename it to `orgin_sq_clac_mean`
    if ("sq_calc_mean" %in% colnames(df_tmp)) {
      colnames(df_tmp)[colnames(df_tmp) == "sq_calc_mean"] <- "orgin_sq_calc_mean"
    }

    # Merge with sample_data
    df_tmp2 = df_tmp %>%
      left_join(df_sq, by = "SampleID")

    # Rownames restore
    rownames(df_tmp2) = rownames(sample_data(psdata_qpcr_norm))
    sample_data(psdata_qpcr_norm) = sample_data(df_tmp2)

  } else if (!is.null(norm_method) && norm_method == "qpcr" && copy_correction == FALSE) {
    log_message(paste("Error: norm_method is set to 'qpcr' but copy_correction is FALSE. The qPCR normalization method requires copy number correction to be enabled."), log_file)
    stop("Error: norm_method is set to 'qpcr' but copy_correction is FALSE. The qPCR normalization method requires copy number correction to be enabled.")
  }

  if (copy_correction == TRUE) {
  # Downloading the rndb-database
  zip_file_path <- file.path(destination_folder, "rrnDB-5.9_pantaxa_stats_NCBI.tsv.zip")
  download.file("https://rrndb.umms.med.umich.edu/downloads/rrnDB-5.9_pantaxa_stats_NCBI.tsv.zip",
                destfile = zip_file_path, mode = "wb")

  # Check if the file exists
  if (!file.exists(zip_file_path) || file.info(zip_file_path)$size == 0) {
    stop("Error: Downloaded ZIP file is missing or corrupted.")
  }

  # unzip the downloaded database
  unzip(zip_file_path, exdir = destination_folder)
  rrndb_database_tsv_file = list.files(destination_folder, pattern = "pantaxa_stats_NCBI\\.tsv$", full.names = TRUE)
  rrndb_database_tsv = read_tsv(rrndb_database_tsv_file)
  rrndb_database = rrndb_database_tsv %>% filter(rank == "genus") %>% select(Genus = "name", everything())

  # preparing the phyloseq data
  # psmelt = psdata %>% psmelt() %>% as_tibble()
  # psmelt = psmelt %>% select(OTU, Genus)

  tax = as.data.frame(tax_table(psdata)) %>%
    tibble::rownames_to_column("OTU") %>%
    select(OTU, Genus) %>% as_tibble()

  proj_tabel = left_join(tax, raspergade_df, by = "OTU")
  final_tabel = left_join(proj_tabel, rrndb_database, by = "Genus")

  plot_data = final_tabel %>%
    select(OTU, Genus, mean, copy_number, probability) %>%
    filter(!is.na(mean) & !is.na(copy_number)) %>%
    mutate(probability_rate = case_when(
      probability > 0.9 ~ "High (> 0.9)",
      probability >= 0.5 & probability <= 0.9 ~ "Medium (>= 0.5 & <= 0.9)",
      probability < 0.5 ~ "Low (< 0.5)"))

  cor_test = cor.test(plot_data$mean, plot_data$copy_number, method = "pearson")
  r_val = cor_test$estimate
  p_val = cor_test$p.value
  label_text = sprintf("r = %.2f", r_val)

  copy_number_comparison =
    ggplot(plot_data, aes(x = mean, y = copy_number, color = probability_rate)) +
    geom_point(aes(color = probability_rate), size = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "black", aes(group = 1)) +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    scale_color_manual(values = c("High (> 0.9)" = "darkgreen", "Medium (>= 0.5 & <= 0.9)" = "orange", "Low (< 0.5)" = "red4")) +
    #facet_wrap(~ probability_rate) +
    xlim(0, 15) +
    ylim(0, 15) +
    theme_bw() +
    labs(
      title = "ASV copy number comparison",
      x = "mean copy number reference rrnDB",
      y = "Predicted copy number per genus",
      color = "Probability rate") +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 12),
          legend.position = c(0.05, 0.95),
          legend.justification = c(0, 1),
          legend.background = element_rect(fill = "white",  colour = "grey80")) +
    annotate("text", x = 2.25, y = 11, label = label_text,
             hjust = 1, vjust = 1, size = 5)

  # Save as PDF
  figure_file_path_pdf <- file.path(figure_folder_png, paste0(project_name, "_copy_number_comparison.pdf"))
  ggsave(filename = figure_file_path_pdf, plot = copy_number_comparison, width = 6, height = 6)

  log_message(paste("Copy number comparison saved as .pdf object in", figure_file_path_pdf), log_file)

  # Save as PNG
  figure_file_path_png <- file.path(figure_folder_pdf, paste0(project_name, "_copy_number_comparison.png"))
  ggsave(filename = figure_file_path_png, plot = copy_number_comparison, width = 6, height = 6, dpi = 600)

  log_message(paste("Copy number comparison saved as .png object in", figure_folder_pdf), log_file)
  }

  if (copy_correction == FALSE) {
  psdata_copy_number_corrected = psdata
  output_file_path = paste0(output_asv_rds_files, project_name, "_phyloseq_asv_level_without_copy_number_corrected_counts.rds")
  saveRDS(psdata, file = output_file_path)
  log_message(paste("phyloseq data without copy number corrected counts asv level saved as .rds object in", output_file_path), log_file)

  } else if (copy_correction == TRUE) {
  output_file_path = paste0(output_asv_rds_files, project_name, "_phyloseq_asv_level_copy_number_corrected_counts.rds")
  saveRDS(psdata_copy_number_corrected, file = output_file_path)
  log_message(paste("phyloseq data copy number corrected counts asv level saved as .rds object in", output_file_path), log_file)

  }

  if (is.null(norm_method)) {
    log_message("No normalization method specified for absolute data. Only copy number correction applied.", log_file)
  }

  if (norm_method == "fcm") {
    if (copy_correction == TRUE) {
    output_file_path = paste0(output_folder_rds_files_before, project_name, "_phyloseq_asv_level_fcm_normalised_cell_concentration.rds")
    saveRDS(psdata_fcm_norm, file = output_file_path)
    log_message(paste("Phyloseq data fcm normalised cell concentration (cells per ml/gram sample) asv level saved as .rds object in", output_file_path), log_file)

    } else if (copy_correction == FALSE) {
    output_file_path = paste0(output_folder_rds_files_before, project_name, "_phyloseq_asv_level_fcm_normalised_cell_concentration_withou_copy_number_corrected_count.rds")
    saveRDS(psdata_fcm_norm, file = output_file_path)
    log_message(paste("Phyloseq data fcm normalised cell concentration (cells per ml/gram sample) asv level saved as .rds object in", output_file_path), log_file)
    }

    return(list(psdata_asv_copy_number_corrected = psdata_copy_number_corrected, psdata_asv_fcm_norm = psdata_fcm_norm))

  } else if (norm_method == "qpcr") {
    output_file_path = paste0(output_folder_rds_files_before, project_name, "_phyloseq_asv_level_qpcr_normalised_cell_concentration.rds")
    saveRDS(psdata_qpcr_norm, file = output_file_path)
    log_message(paste("phyloseq qpcr normalised cell concentration (cells per ml/gram sample) asv level saved as .rds object in", output_file_path), log_file)

    return(list(psdata_asv_copy_number_corrected = psdata_copy_number_corrected, psdata_asv_qpcr_norm = psdata_qpcr_norm))

  } else {
    log_message("Error: Invalid normalization method specified. Use 'fcm' or 'qpcr'.")
    stop("Error: Invalid normalization method specified. Use 'fcm' or 'qpcr'.")
  }

  if (copy_correction == TRUE) {
    log_message("Copy number correction successfully applied.", log_file)
  } else if (copy_correction == FALSE) {
    log_message("Copy number correction has not been applied.", log_file)
  }

}


