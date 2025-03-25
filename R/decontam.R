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
decontam =  function(physeq = resolved_tree_physeq,
                     decon_method = c("frequency", "prevalence", "both"),
                     blank = TRUE) {

  psdata = physeq
  project_name = projects

  project_folder = paste0(base_path, project_name)
  figure_folder = paste0(project_folder, "/figures/")
  output_folder_csv_files = paste0(project_folder, "/output_data/csv_files/")
  output_folder_rds_files = paste0(project_folder, "/output_data/rds_files/Before_cleaning_rds_files/")

  contam_taxa_freq = character(0)
  contam_taxa_prev = character(0)

  if (!("sample_or_control" %in% colnames(sample_data(psdata)))) {
    error_message = paste0("error: 'sample_or_control' column is missing from the sample data.")
    log_message(error_message, log_file)
    stop(error_message)
  }

  # error if decon_method is "both" but blank is FALSE or missing
  if (decon_method == "both" && (!blank || !"blank" %in% sample_data(psdata)$sample_or_control)) {
    error_message = paste0("error: 'decon_method' set to 'both' but 'blank' samples are either not present or blank parameter is FALSE.")
    log_message(error_message, log_file)
    stop(error_message)
  }

  # # error if blank is missing
  # if (!"blank" %in% sample_data(psdata)$sample_or_control) {
  #   error_message = paste0("error: blank sample is not present.")
  #   log_message(error_message, log_file)
  #   stop(error_message)
  # }

  # adjust decon_method based on blank parameter
  if (blank == FALSE) {
    decon_method = "frequency"
    message = paste0("message: 'blank' parameter is FALSE using method 'frequency'.")
    log_message(message, log_file)
  }

  # Decontaminate with frequency method
  if (decon_method == "frequency" || decon_method == "both") {
    contamdf.freq <- isContaminant(psdata, method = "frequency", conc = "DNA_Concentration")
    print(paste(project_name, "- Method frequency:"))
    if (any(contamdf.freq$contaminant)) {
      print(table(contamdf.freq$contaminant))
      contam_taxa_freq = rownames(contamdf.freq)[contamdf.freq$contaminant == TRUE]
    } else {
      message = paste0("No contaminants found using frequency method.")
      log_message(message, log_file)
    }
  }

  # Decontaminate with prevalence method
  if (decon_method == "prevalence" || decon_method == "both") {
    sample_data(psdata)$is.neg = sample_data(psdata)$sample_or_control == "blank"
    contamdf.prev = isContaminant(psdata, method = "prevalence", neg = "is.neg", threshold = 0.5)
    print(paste(project_name, "- Method prevalence:"))
    if (any(contamdf.prev$contaminant)) {
      print(table(contamdf.prev$contaminant))
      contam_taxa_prev = rownames(contamdf.prev)[contamdf.prev$contaminant == TRUE]
    } else {
      message = paste0("No contaminants found using prevalence method.")
      log_message(message, log_file)
    }
  }

  # Combine both methods
  if (decon_method == "both") {
    all_OTUs = union(rownames(contamdf.freq), rownames(contamdf.prev))
    contamdf.both = tibble(
      OTU = all_OTUs,
      frequency_contaminant = contamdf.freq$contaminant[match(all_OTUs, rownames(contamdf.freq))],
      prevalence_contaminant = contamdf.prev$contaminant[match(all_OTUs, rownames(contamdf.prev))])
    contamdf.both = contamdf.both %>%
      mutate(both_contaminant = frequency_contaminant | prevalence_contaminant)
    contamdf.both_unique = contamdf.both %>%
      distinct(OTU, .keep_all = TRUE) %>%
      mutate(contaminant = both_contaminant)
    contam_taxa = unique(c(contam_taxa_freq, contam_taxa_prev))
  } else if (decon_method == "frequency") {
    contam_taxa = contam_taxa_freq
  } else if (decon_method == "prevalence") {
    contam_taxa = contam_taxa_prev
  }

  physeq_no_contam = prune_taxa(!taxa_names(psdata) %in% contam_taxa, psdata)

  output_file_path = paste0(output_folder_rds_files, project_name, "_phyloseq_asv_level_decontam.rds")
  saveRDS(physeq_no_contam, file = output_file_path)
  log_message(paste("Decontam phyloseq object saved as .rds object in", output_file_path), log_file)

  df = as.data.frame(sample_data(psdata))
  df$read_count = sample_sums(psdata)
  df = df[order(df$read_count),]
  df$Index = seq(nrow(df))

  p1 = ggplot(df, aes(x = Index, y = read_count, color = sample_or_control)) +
    geom_point() +
    xlab("Sample Index") +
    ylab("Read Count") +
    ggtitle(paste("Read Count Plot for", project_name))

  figure_file_path = paste0(figure_folder, project_name, "_decontam_library_size.pdf")
  ggsave(filename = figure_file_path, plot = p1)
  log_message(paste("Read Count plot saved as .pdf object in", figure_file_path), log_file)

  print(p1)

  if (blank == TRUE) {

    presence_absence = transform_sample_counts(psdata, function(Abundance) 1 * (Abundance > 0))
    presence_absence_neg = prune_samples(sample_data(presence_absence)$sample_or_control == "blank", presence_absence)
    presence_absence_pos = prune_samples(sample_data(presence_absence)$sample_or_control != "blank", presence_absence)

    if (decon_method == "frequency") {
      presence_absence_df = data.frame(
        presence_absence_pos = taxa_sums(presence_absence_pos),
        presence_absence_neg = taxa_sums(presence_absence_neg),
        contaminant = contamdf.freq$contaminant
      )
    } else if (decon_method == "prevalence") {
      presence_absence_df = data.frame(
        presence_absence_pos = taxa_sums(presence_absence_pos),
        presence_absence_neg = taxa_sums(presence_absence_neg),
        contaminant = contamdf.prev$contaminant
      )
    } else if (decon_method == "both") {
      presence_absence_df = data.frame(
        presence_absence_pos = taxa_sums(presence_absence_pos),
        presence_absence_neg = taxa_sums(presence_absence_neg),
        contaminant = contamdf.both_unique$contaminant
      )
    }

    p2 = ggplot(presence_absence_df, aes(x = presence_absence_neg, y = presence_absence_pos, color = contaminant)) +
      geom_point() +
      xlab("Prevalence (Negative Controls)") +
      ylab("Prevalence (True Samples)") +
      ggtitle(paste(decon_method, "Plot for", project_name))

    figure_file_path = paste0(figure_folder, project_name, "_decontam_method_", decon_method, ".pdf")

    ggsave(filename = figure_file_path, plot = p2)
    log_message(paste("Plot saved as .pdf object in", figure_file_path), log_file)
    print(p2)
  }

  return(physeq_no_contam)
}
