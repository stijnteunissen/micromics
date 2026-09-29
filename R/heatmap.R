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
heatmap = function(psdata, taxrank, ntaxa = 23, facet_vars, plot_width, project_id, base_path, log_file) {

  # Define internal heatmap function
  base_heatmap = function(plot_data, x_value, abund_value, legend_name, x_label = "Sample", tax_column, facet_vars) {

    # Initialize ggplot heatmap with dynamci taxonomic column
    p <- ggplot(plot_data, aes(x = Sample,
                               y = !!sym(tax_column))) +
      geom_tile(aes(fill = !!sym(abund_value)), color = NA) +
      scale_fill_gradient(low = "white", high = "darkred", name = legend_name) +
      labs(x = x_label,
           y = tax_column) +
      theme_classic() +
      theme(axis.text.y = element_markdown(size = 10),
            strip.placement = "outside",
            strip.background = element_blank(),
            #strip.text = element_text(face = "bold", angle = 90, vjust = 0.5, hjust = 0),
            strip.text = element_text(face = "bold"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
            ggh4x.facet.nestline = element_line(colour = "black")) +
      scale_x_discrete(expand = c(0, 0)) +
      geom_text(aes(label = ifelse(!!sym(abund_value) > 3, paste0(round(!!sym(abund_value), 0)), ifelse(!!sym(abund_value) == 0, ".", "")),
                    color = ifelse(!!sym(abund_value) > 50, "#D3D3D3", "black")),
                size = 3)

    # Add nested facets condittionally if vairables are provided
    if (!is.null(facet_vars)) {
      p <- p +
        facet_nested(
          cols = vars(!!!syms(facet_vars)),
          scales = "free_x",
          space = "free_x",
          nest_line = element_line(linetype = 1, color = "black", linewidth = 0.6),
          strip = strip_nested(size = "variable")
        )
    }
    return(p)
  }

  # Set up directory paths
  figures_folder <- file.path(base_path, "03_figures")
  heatmap_folder <- file.path(figures_folder, "heatmap")
  if(!dir.exists(heatmap_folder)) { dir.create(heatmap_folder, recursive = TRUE) }

  # Process data sequentially for each taxonomic rank
  for (tax in taxrank) {
    psdata_rmp <- psdata[[paste0("psdata_rmp_", tax)]]

    # Summarise abundances and calculate relative percentages
    abund_rmp <- psdata_rmp %>%
      dplyr::group_by(Sample, !!sym(tax), na_type, across(-c(Abundance, OTU, read_count))) %>%
      dplyr::summarise(abund = sum(Abundance), .groups = "drop") %>%
      dplyr::group_by(Sample) %>%
      dplyr::mutate(rel_abund = abund/sum(abund) * 100) %>%
      dplyr::ungroup()

    # Determine cut-off threshold on top abundance taxa
    legend_cutoff_rmp <- abund_rmp %>%
      dplyr::group_by(!!sym(tax)) %>%
      dplyr::summarise(mean_rel_abund = mean(rel_abund, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(desc(mean_rel_abund)) %>%
      dplyr::slice_head(n = ntaxa)

    # Define label for low abundance taxa
    legend_cutoff_value_rmp <- legend_cutoff_rmp %>% pull(mean_rel_abund) %>% min()
    other_label <- glue("Other max.<{round(legend_cutoff_value_rmp, 2)}%")

    # Reorder and lump lower abundace taxa into other cateory
    plot_data_rmp <- abund_rmp %>%
      dplyr::group_by(!!sym(tax)) %>%
      dplyr::mutate(total_genus_abund = mean(rel_abund, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        !!sym(tax) := fct_reorder(!!sym(tax), total_genus_abund, .desc = TRUE),
        !!sym(tax) := fct_lump_n(!!sym(tax), n = ntaxa, w = total_genus_abund, other_level = other_label)) %>%
      dplyr::group_by(Sample, !!sym(tax), na_type, across(-c(total_genus_abund, rel_abund))) %>%
      dplyr::summarise(rel_abund = sum(rel_abund), .groups = "drop")

    # Clean an italicize names specifically for genus level
    if (tax == "Genus") {
      plot_data_rmp_cleaned <- plot_data_rmp %>%
        dplyr::group_by(Sample, Genus) %>%
        dplyr::mutate(
          Genus = if_else(str_detect(Genus, other_label), as.character(Genus),
                          case_when(
                            str_detect(Genus, "Genus of Candidatus (\\S+) bacterium (\\S+)") ~ str_replace(Genus, "Genus of Candidatus (\\S+) bacterium (\\S+)", "Genus of *Candidatus* *\\1* (\\2)"),
                            str_detect(Genus, "Genus of Candidatus (\\S+)") ~ str_replace(Genus, "Genus of Candidatus (\\S+)", "Genus of *Candidatus* *\\1*"),
                            str_detect(Genus, "(.*)_unclassified") ~ str_replace(Genus, "(.*)_unclassified", "Unclassified *\\1*"),
                            str_detect(Genus, "Genus of") ~ str_replace(Genus, "Genus of (\\S+)", "Genus of *\\1*"),
                            str_detect(Genus, "(\\S+)\\s+(\\S+)") ~ str_replace(Genus, "(\\S+)\\s+(\\S+)", "*\\1* (*\\2*)"),
                            TRUE ~ str_replace(Genus, "^(\\S*)$", "*\\1*")))) %>%
        dplyr::ungroup()
    } else {
      plot_data_rmp_cleaned <- plot_data_rmp
    }

    # Extract final character order of features excluding other category
    tax_order <- plot_data_rmp_cleaned %>%
      dplyr::group_by(!!sym(tax)) %>%
      dplyr::summarise(avg = mean(rel_abund), .groups = "drop") %>%
      dplyr::filter(!!sym(tax) != other_label) %>%
      dplyr::arrange(desc(avg)) %>%
      dplyr::pull(!!sym(tax)) %>%
      as.character()

    # Re-level factor to lock in visual display order
    final_levels <- c(tax_order, other_label)
    plot_data_rmp_cleaned <- plot_data_rmp_cleaned %>%
      dplyr::mutate(!!sym(tax) := factor(!!sym(tax), levels = final_levels))

    # Call base function to assemble heatmap plot object
    heatmap_rmp <- base_heatmap(plot_data_rmp_cleaned, "Sample", "rel_abund", legend_name = "Relative\nAbundance (%)", x_label = "Sample", tax, facet_vars) +
      scale_color_identity() +
      theme(legend.position = "right")

    # figures saved as png
    ggsave(filename = file.path(heatmap_folder, glue::glue("heatmap_rmp_{tax}.png")), plot = heatmap_rmp, width = plot_width, height = 10, dpi = 600)
    ggsave(filename = file.path(heatmap_folder, glue::glue("heatmap_rmp_{tax}.pdf")), plot = heatmap_rmp, width = plot_width, height = 10)
  }
}
