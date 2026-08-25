#' Decontaminate a Phyloseq Object Using Specified Methods
#'
#' This function removes contamination from a phyloseq object using the
#' [`decontam`](https://benjjneb.github.io/decontam/vignettes/decontam_intro.html) package.
#' It supports the frequency, prevalence, or both methods for contaminant identification.
#' Note that the prevalence method relies on the presence of blank samples;
#' if blank samples are not available, only the frequency method can be performed.
#'
#' @inheritParams resolve_tree
#' @param decon_method A character string specifying the contamination removal method.
#'   Possible values are:
#'   \itemize{
#'     \item `frequency`: Identifies contaminants by examining the distribution of sequence feature frequencies
#'           as a function of the input DNA concentration.
#'     \item `prevalence`: Identifies contaminants by comparing the prevalence (presence/absence across samples)
#'           of sequence features in true samples versus negative controls (blanks).
#'     \item `both`: Applies both frequency and prevalence methods sequentially.
#'           Taxa flagged by either method are considered contaminants.
#'           This option requires that blank samples are available.
#'   }
#' @param blank A logical value indicating whether blank samples were included in the dataset.
#'   \itemize{
#'     \item `TRUE`: Blank samples are present, allowing the use of the `both` method.
#'     \item `FALSE`: No blank samples are available; in this case, only the `frequency` method can be applied.
#'   }
#'
#' @return A phyloseq object with contaminants removed. The decontaminated object is saved as an RDS file named
#' `<project_name>_phyloseq_asv_level_decontam.rds` in the `output_data/rds_files/Before_cleaning_rds_files` directory.
#'
#' @details
#' The function uses the `decontam` package to remove contaminants based on the specified method.
#' If `both` is chosen, the function applies the frequency method first and then the prevalence method,
#' flagging any taxa identified by either method as contaminants.
#' In addition, diagnostic plots (showing read counts and contaminant prevalence) are generated and saved as PDF files.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' decontam(physeq = phyloseq_data, decon_method = "both", blank = TRUE)
#' }
#'
#' @export
decontam =  function(physeq, decon_method = c("frequency", "prevalence", "both"), blank = TRUE, project_id, base_path, log_file) {

  log_message("Starting microbiome decontamination.", status = "start", log_file)

  # Match argument input choices
  decon_method <- match.arg(decon_method)

  # Define coutput directory paths
  figure_folder <- file.path(base_path, "figures")
  raw_rds_folder <- file.path(base_path, "raw_rds")

  # Initialize empty containers for tracking contaminant
  contam_taxa_freq <- character(0)
  contam_taxa_prev <- character(0)

  # Validate presence of required sample data
  if (!("sample_or_control" %in% colnames(phyloseq::sample_data(physeq)))) {
    error_message <- paste0("Error: 'sample_or_control' column is missing from the sample data.")
    log_message(error_message, status = "error", log_file)
    stop(error_message, call. = FALSE)
  }

  # Stop if combination (both) filtering is requested but blank is missing
  if (decon_method == "both" && (!blank || !"blank" %in% phyloseq::sample_data(physeq)$sample_or_control)) {
    error_message <- paste0("Error: 'decon_method' set to 'both' but 'blank' sample is either not present or blank parameter is FALSE.")
    log_message(error_message, status = "error", log_file)
    stop(error_message, call. = FALSE)
  }

  # Adjust 'decon_method' based on blank parameter
  if (blank == FALSE) {
    decon_method <- "frequency"
    log_message("The 'blank' parameter is FALSE. Forcing decontamination method to 'frequency'.", status = "warning", log_file)
  }

  # Execute Frequency based decontamination -------------------------------
  if (decon_method %in% c("frequency", "both")) {
    contam_df_freq <- decontam::isContaminant(physeq, method = "frequency", conc = "DNA_Concentration")

    if (any(contam_df_freq$contaminant)) {
      found_count <- sum(contam_df_freq$contaminant)
      log_message(glue::glue("Frequency method identified {found_count} contaminant ASVs."), status = "info", log_file)
      contam_taxa_freq <- rownames(contam_df_freq)[contam_df_freq$contaminant == TRUE]
    } else {
      log_message("No contaminants detected via DNA concentration frequency patterns.", status = "info", log_file)
    }
  }

  # Execute Prevalence based decontamination ------------------------------
  if (decon_method %in% c("prevalence", "both")) {
    phyloseq::sample_data(physeq)$is_neg <- phyloseq::sample_data(physeq)$sample_or_control == "blank"
    contam_df_prev <- decontam::isContaminant(physeq, method = "prevalence", neg = "is_neg", threshold = 0.5)

    if (any(contam_df_prev$contaminant)) {
      found_count <- sum(contam_df_prev$contaminant)
      log_message(glue::glue("Prevalance method identified {found_count} contaminant ASVs."), status = "info", log_file)
      contam_taxa_prev <- rownames(contam_df_prev)[contam_df_prev$contaminant == TRUE]
    } else {
      log_message("No contaminants detected via negative control prevalence checks.", status = "info", log_file)
    }
  }

  # Combine filtering arrays ----------------------------------------------
  if (decon_method == "both") {
    all_OTUs = union(rownames(contam_df_freq), rownames(contam_df_prev))

    contam_df_both <- tibble::tibble(
      OTU = all_OTUs,
      frequency_contaminant = contam_df_freq$contaminant[match(all_OTUs, rownames(contam_df_freq))],
      prevalence_contaminant = contam_df_prev$contaminant[match(all_OTUs, rownames(contam_df_prev))]
    ) %>%
      dplyr::mutate(both_contaminant = frequency_contaminant | prevalence_contaminant)

    contam_df_both_unique <- contam_df_both %>%
      dplyr::distinct(OTU, .keep_all = TRUE) %>%
      dplyr::mutate(contaminant = both_contaminant)

    contam_taxa <- unique(c(contam_taxa_freq, contam_taxa_prev))

  } else if (decon_method == "frequency") {
    contam_taxa <- contam_taxa_freq
  } else if (decon_method == "prevalence") {
    contam_taxa <- contam_taxa_prev
  }

  # Prune identified contaminants out of the dataset
  physeq_no_contam <- phyloseq::prune_taxa(!phyloseq::taxa_names(physeq) %in% contam_taxa, physeq)

  # Save results
  output_rds_path <- file.path(raw_rds_folder, glue::glue("{project_id}_phyloseq_asv_decontam.rds"))
  saveRDS(physeq_no_contam, file = output_rds_path)
  log_message(glue::glue("Decontaminated object successfully stored at: {output_rds_path}"), status = "info", log_file)

  # Generate quality control plot -----------------------------------------
  metadata <- as.data.frame(as(phyloseq::sample_data(physeq), "data.frame"))
  library_df <- metadata %>%
    dplyr::mutate(
      read_count = phyloseq::sample_sums(physeq),
      index = dplyr::row_number(read_count)
    )

  plot_library_size <- ggplot2::ggplot(library_df, aes(x = index, y = read_count, color = sample_or_control)) +
    geom_point() +
    labs(
      x = "Sample Index",
      y = "Read Count",
      title = "Read Count Plot"
    ) +
    theme_minimal()

  ggplot2::ggsave(filename = file.path(figure_folder, glue::glue("{project_id}_decontam_library_size.png")), plot = plot_library_size, width = 8, height = 8)
  ggplot2::ggsave(filename = file.path(figure_folder, glue::glue("{project_id}_decontam_library_size.pdf")), plot = plot_library_size, width = 8, height = 8)

  # Generate Prevalence Plots if blanks are present
  if (blank == TRUE) {
    presence_absence = phyloseq::transform_sample_counts(physeq, function(Abundance) 1 * (Abundance > 0))
    presence_absence_neg = phyloseq::prune_samples(phyloseq::sample_data(presence_absence)$sample_or_control == "blank", presence_absence)
    presence_absence_pos = phyloseq::prune_samples(phyloseq::sample_data(presence_absence)$sample_or_control != "blank", presence_absence)

    # Establish plotting target veactors dynamically
    contaminant_logical <- switch(decon_method,
                                  "frequency" = contam_df_freq$contaminant,
                                  "prevalence" = contam_df_prev$contaminant,
                                  "both" = contam_df_both_unique$contaminant)

    presence_absence_df <- data.frame(
      presence_absence_pos = phyloseq::taxa_sums(presence_absence_pos),
      presence_absence_neg = phyloseq::taxa_sums(presence_absence_neg),
      contaminant = contaminant_logical
    )

    plot_prevalence = ggplot2::ggplot(presence_absence_df, aes(x = presence_absence_neg, y = presence_absence_pos, color = contaminant)) +
      geom_point() +
      labs(
        x = "Prevalence (Negative Controls)",
        y = "Prevalence (True Samples)",
        title = glue::glue("{stringr::str_to_title(decon_method)} Isolation Plot: {project_id}")
      ) +
      theme_minimal()

    # Save decontam method plot
    ggplot2::ggsave(filename = file.path(figure_folder, glue::glue("{project_id}_decontam_method_{decon_method}.png")), plot = plot_prevalence, width = 8, height = 8)
    ggplot2::ggsave(filename = file.path(figure_folder, glue::glue("{project_id}_decontam_method_{decon_method}.pdf")), plot = plot_prevalence, width = 8, height = 8)

    log_message("Decontamination diagnostic visualizations generated and exported.", status = "info", log_file)
  }

  log_message("Microbiome decontamination protocol completed successfully.", status = "success", log_file)
  return(physeq_no_contam)
}
