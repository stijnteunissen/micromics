#' Generate a Heatmap of Relative Abundance
#'
#' This function creates a heatmap of relative abundance data from a phyloseq
#' object (or similar data frame) at the genus level. It calculates the relative
#' abundance of each taxon per sample and groups taxa with low relative
#' abundance (below a defined threshold) into an "Other" category. The heatmap
#' is then facetted based on additional sample metadata if available and saved
#' as a PDF.
#'
#' @param physeq A phyloseq object containing normalized genus-level data. The default is \code{rarefied_genus_psmelt}.
#' @param ntaxa An integer specifying the maximum number of taxa to display individually. Taxa below the threshold are grouped into "Other". If \code{NULL}, \code{ntaxa} is set to 23.
#' @param norm_method A character string specifying the normalization method. If \code{NULL}, the function uses the provided \code{physeq} directly. If set to \code{"fcm"} or \code{"qpcr"}, the function extracts the corresponding \code{psmelt_copy_number_corrected_} data based on the taxonomic rank.
#' @param taxrank A character string indicating the taxonomic rank to use for grouping taxa. The default is \code{"Genus"}.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Sets up project folder paths for figures and output data.
#'   \item Extracts and processes the input data to compute the relative abundance (in percentage) of each taxon per sample.
#'   \item Groups taxa with a mean relative abundance below a defined cutoff into an "Other" category.
#'   \item Optionally orders the data by \code{Sample_Date} if that factor is present in the metadata.
#'   \item Creates a base heatmap using \code{ggplot2}, with samples on the x-axis and taxa on the y-axis. The fill color reflects the relative abundance, and text labels are added for values exceeding a threshold.
#'   \item If more than one \code{na_type} is present (e.g., both DNA and RNA), separate heatmaps are generated for each and then combined.
#'   \item Saves the final heatmap as a PDF file in the project's figures folder.
#' }
#'
#' @return A \code{ggplot} object representing the heatmap of relative abundance.
#'
#' @examples
#' \dontrun{
#'   # Generate a heatmap using default parameters
#'   heatmap_plot <- heatmap(physeq = rarefied_genus_psmelt)
#'
#'   # Generate a heatmap with a specified number of taxa and a normalization method
#'   heatmap_plot <- heatmap(
#'     physeq = rarefied_genus_psmelt,
#'     ntaxa = 20,
#'     norm_method = "fcm",
#'     taxrank = c("Phylum", "Class", "Order", "Family", "Genus")
#'   )
#' }
#'
#' @export
heatmap = function(physeq = rarefied_genus_psmelt,
                   ntaxa = NULL,
                   norm_method = NULL,
                   taxrank = c("Phylum", "Class", "Order", "Family", "Genus"),
                   date_factor = NULL) {

  log_message(paste("Step 14: Creating heatmap.", paste(projects, collapse = ", ")), log_file)

  base_heatmap = function(plot_data, x_value, abund_value, legend_name, x_label = "Sample", tax_column) {
    ggplot(plot_data, aes(x = Sample,
                          y = !!sym(tax_column))) +
      geom_tile(aes(fill = !!sym(abund_value)), color = NA) +
      scale_fill_gradient(low = "white", high = "darkred", name = legend_name) +
      labs(x = x_label,
           y = tax_column) +
      theme_classic() +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_markdown(size = 10),
            axis.ticks.x = element_blank(),
            strip.placement = "outside",
            strip.background = element_blank(),
            #strip.text = element_text(face = "bold", angle = 90, vjust = 0.5, hjust = 0),
            strip.text = element_text(face = "bold"),
            ggh4x.facet.nestline = element_line(colour = "black")) +
      scale_x_discrete(expand = c(0, 0)) +
      geom_text(aes(label = ifelse(mean_rel_abund > 3, paste0(round(mean_rel_abund, 0)), ifelse(mean_rel_abund == 0, ".", "")),
                    color = ifelse(mean_rel_abund > 50, "#D3D3D3", "black")),
                size = 3)
  }

  facet_add = function(present_factors, include_na_type = FALSE) {
    if (include_na_type) {
      factors_to_add = c("na_type", present_factors)
    } else {
      factors_to_add = present_factors
    }

    if (!is.null(factors_to_add) && length(factors_to_add) > 0) {
      return(facet_nested(cols = vars(!!!syms(factors_to_add)), scales = "free", space = "free", nest_line = element_line(linetype = 1)))
    } else {
      return(NULL)
    }
  }

  project_name = projects

  project_folder = paste0(base_path, project_name)
  figure_folder_pdf = paste0(project_folder, "/figures/PDF_figures/")
  if(!dir.exists(figure_folder_pdf)) { dir.create(figure_folder_pdf, recursive = TRUE) }
  figure_folder_png = paste0(project_folder, "/figures/PNG_figures/")
  if(!dir.exists(figure_folder_png)) { dir.create(figure_folder_png, recursive = TRUE) }

  if (is.null(ntaxa)) {
    ntaxa = 23
  }

  for (tax in taxrank) {
    log_message(paste("Processing taxonomic level:", tax), log_file)

    # Choose the appropriate dataset based on norm_method
    if (is.null(norm_method)) {
      copy_number_corrected_data = physeq[[paste0("psmelt_copy_number_corrected_", tax)]]
    } else if (norm_method == "fcm" || norm_method == "qpcr") {
      copy_number_corrected_data = physeq[[paste0("psmelt_copy_number_corrected_", tax)]]
    }

    # Create output folders (barplot folder and tax folder)
    heatmap_folder_png = paste0(figure_folder_png, "Heatmap/")
    if(!dir.exists(heatmap_folder_png)) { dir.create(heatmap_folder_png, recursive = TRUE) }
    tax_folder_png = paste0(heatmap_folder_png, tax, "/")
    if(!dir.exists(tax_folder_png)) { dir.create(tax_folder_png, recursive = TRUE) }

    heatmap_folder_pdf = paste0(figure_folder_pdf, "Heatmap/")
    if(!dir.exists(heatmap_folder_pdf)) { dir.create(heatmap_folder_pdf, recursive = TRUE) }
    tax_folder_pdf = paste0(heatmap_folder_pdf, tax, "/")
    if(!dir.exists(tax_folder_pdf)) { dir.create(tax_folder_pdf, recursive = TRUE) }

  variable_columns = intersect(present_variable_factors, colnames(copy_number_corrected_data))
  factor_columns = unique(c(variable_columns))
  present_factors = if (length(factor_columns) > 0) factor_columns else NULL

  # relatieve abudnatie onder 1% wordt geroepeerd onder "Others"
  genus_abund_rel =
    copy_number_corrected_data %>%
    group_by(Sample, !!!syms(tax), na_type, !!!syms(present_factors)) %>%
    summarise(abund = sum(Abundance), .groups = "drop") %>%
    group_by(Sample) %>%
    mutate(mean_rel_abund = abund/sum(abund) * 100) %>%
    ungroup() %>%
    group_by(Sample, !!sym(tax)) %>%
    mutate(!!sym(tax) := str_replace(!!sym(tax), "(.*)_unclassified", "Unclassified *\\1*")) %>%
    mutate(!!sym(tax) := case_when(
      str_detect(!!sym(tax), "Genus of") ~ str_replace(!!sym(tax), "Genus of (\\S+)", "Genus of *\\1*"),
      str_detect(!!sym(tax), "(\\S+)\\s+(\\S+)") ~ str_replace(!!sym(tax), "(\\S+)\\s+(\\S+)", "*\\1* (*\\2*)"),
      TRUE ~ str_replace(!!sym(tax), "^(\\S*)$", "*\\1*")
    ))

  if (!is.null(date_factor) && date_factor %in% present_factors) {
    genus_abund_rel <- genus_abund_rel %>%
      mutate(!!sym(date_factor) := as.Date(!!sym(date_factor), format = "%d/%m/%Y")) %>%
      arrange(!!sym(date_factor))
  }

  legend_cutoff_rel =
    genus_abund_rel %>%
    group_by(!!sym(tax)) %>%
    summarise(mean_rel_abund = mean(mean_rel_abund, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_rel_abund)) %>%
    slice_head(n = ntaxa)

  legend_cutoff_names_rel = legend_cutoff_rel %>% pull(!!sym(tax))
  legend_cutoff_value_rel = legend_cutoff_rel %>% pull(mean_rel_abund) %>% min()

  genus_pool_rel =
    genus_abund_rel %>%
    group_by(!!sym(tax)) %>%
    summarise(pool = max(mean_rel_abund) < legend_cutoff_value_rel,
              mean = mean(mean_rel_abund), .groups = "drop")

  plot_data_rel =
    inner_join(genus_abund_rel, genus_pool_rel, by = tax) %>%
    mutate(!!sym(tax) := if_else(!!sym(tax) %in% legend_cutoff_names_rel,
                                         !!sym(tax),
                                         glue("Other max.<{round(legend_cutoff_value_rel, 1)}%"))) %>%
    group_by(Sample, !!!syms(tax), na_type, !!!syms(present_factors)) %>%
    summarise(mean_rel_abund = sum(mean_rel_abund),
              median = median(mean), .groups = "drop") %>%
    mutate(!!sym(tax) := factor(!!sym(tax)),
           !!sym(tax) := fct_reorder(!!sym(tax), median, .desc = TRUE))

  na_types = unique(plot_data_rel$na_type)

  if (length(na_types) == 1) {
    heatmap_relative =
      base_heatmap(plot_data_rel, "Sample", "mean_rel_abund", legend_name = "Relative\nAbundance (%)", x_label = "Sample", tax) +
      scale_color_identity() +
      facet_add(present_factors)

    n_samples <- length(unique(plot_data_rel$Sample))
    fig.width <- max(12, n_samples * 0.40)

    figure_file_path = paste0(tax_folder_png, project_name, "_heatmap_relative_", tax, "_level.png")
    ggsave(filename = figure_file_path, plot = heatmap_relative, width = fig.width, height = 10, limitsize = FALSE, dpi = 600)
    log_message(paste("Relative heatmap saved as .png object in", figure_file_path), log_file)

    figure_file_path = paste0(tax_folder_pdf, project_name, "_heatmap_relative_", tax, "_level.pdf")
    ggsave(filename = figure_file_path, plot = heatmap_relative, width = fig.width, height = 10, limitsize = FALSE)
    log_message(paste("Relative heatmap saved as .pdf object in", figure_file_path), log_file)

  } else if (length(na_types) == 2) {
    plot_data_rel_dna =
      plot_data_rel %>%
      filter(na_type == "dna")

    heatmap_relative_dna =
      base_heatmap(plot_data_rel_dna, "Sample", "mean_rel_abund", legend_name = "Relative\nAbundance (%)", x_label = "Sample", tax) +
      scale_color_identity() +
      facet_add(present_factors, include_na_type = TRUE) +
      theme(legend.position = "right")

    plot_data_rel_rna =
      plot_data_rel %>%
      filter(na_type == "rna")

    heatmap_relative_rna =
      base_heatmap(plot_data_rel_rna, "Sample", "mean_rel_abund", legend_name = "Relative\nAbundance (%)", x_label = "Sample", tax) +
      scale_color_identity() +
      facet_add(present_factors, include_na_type = TRUE) +
      theme(legend.position = "right")

    num_samples_dna = plot_data_rel_dna %>% distinct(Sample) %>% nrow()
    num_samples_rna = plot_data_rel_rna %>% distinct(Sample) %>% nrow()

    dna_width = num_samples_dna
    rna_width = num_samples_rna + 1
    total_width = dna_width + rna_width

    rel_widths = c(dna_width / total_width, rna_width / total_width)

    legend = get_legend(heatmap_relative_dna + theme(legend.position = "right"))

    combined_heatmap_rel = plot_grid(heatmap_relative_dna, heatmap_relative_rna,
                                     ncol = 2, labels = c("A", "B"), rel_widths = rel_widths)

    heatmap_relative = plot_grid(combined_heatmap_rel, legend, ncol = 2, rel_widths = c(3, 0.5))

    n_samples_dna <- length(unique(plot_data_rel_dna$Sample))
    fig.width_dna <- max(12, n_samples_dna * 0.40)

    n_samples_rna <- length(unique(plot_data_rel_rna$Sample))
    fig.width_rna <- max(12, n_samples_rna * 0.40)

    figure_file_path = paste0(tax_folder_png, project_name, "_heatmap_relative_dna_", tax, "_level.png")
    ggsave(filename = figure_file_path, plot = heatmap_relative_dna, width = fig.width_dna, height = 10, limitsize = FALSE, dpi = 600)
    log_message(paste("Relative heatmap saved as .png object in", figure_file_path), log_file)

    figure_file_path = paste0(tax_folder_pdf, project_name, "_heatmap_relative_dna_", tax, "_level.pdf")
    ggsave(filename = figure_file_path, plot = heatmap_relative_dna, width = fig.width_dna, height = 10, limitsize = FALSE)
    log_message(paste("Relative heatmap saved as .pdf object in", figure_file_path), log_file)

    figure_file_path = paste0(tax_folder_png, project_name, "_heatmap_relative_rna_", tax, "_level.png")
    ggsave(filename = figure_file_path, plot = heatmap_relative_rna, width = fig.width_rna, height = 10, limitsize = FALSE, dpi = 600)
    log_message(paste("Relative heatmap saved as .png object in", figure_file_path), log_file)

    figure_file_path = paste0(tax_folder_pdf, project_name, "_heatmap_relative_rna_", tax, "_level.pdf")
    ggsave(filename = figure_file_path, plot = heatmap_relative_rna, width = fig.width_rna, height = 10, limitsize = FALSE)
    log_message(paste("Relative heatmap saved as .pdf object in", figure_file_path), log_file)
    }
  }
  log_message("Heatmap successfully plotted.", log_file)
}

