#' Create Barplots for Relative and Absolute Abundance
#'
#' This function generates barplots for microbial data at the genus level. It supports both
#' relative and absolute abundance data and can include facets based on available metadata factors.
#' The resulting plots can be saved as PDF files, and the underlying data can be exported as
#' CSV and RDS files.
#'
#'
#' @param ntaxa An integer specifying the maximum number of taxa to display in the barplot. Default is 23.
#' @param norm_method A string indicating the normalization method used for absolute abundance
#' data. Options are `"fcm"` (flow cytometry) or `"qpcr"` (quantitative PCR).
#' (relative abundance only).
#' @param sample_matrix An optional matrix specifying the sample structure or metadata.
#' @param group_by_factor with this option you can separtate de barplot for factors
#'
#' @details
#' - Relative abundance plots show proportions of taxa in each sample, with taxa having a mean
#' relative abundance below 1% grouped as "Other".
#' - Absolute abundance plots use normalized cell equivalents (`norm_method = "fcm"` or `"qpcr"`)
#' to display the number of cells per mL for each taxon.
#' - Facets are added based on metadata factors present in the phyloseq object.
#' - Taxa labels are styled to include genus and species names, if available.
#'
#' @return The function generates and saves barplots as PDF files in the project’s `figures/` folder.
#' It also saves the processed data as CSV and RDS files in the corresponding `output_data/` folders.
#' Additionally, the function outputs the plot object for further customization if needed.
#'
#' @examples
#' # Example usage
#' barplot(
#'   physeq = rarefied_genus_psmelt,
#'   ntaxa = 20,
#'   colorset = my_colors,
#'   norm_method = "fcm",
#'   sample_matrix = sample_metadata
#' )
#'
#' @export
barplot <- function(psdata, norm_method, taxrank, ntaxa = 23, facet_vars, plot_width, project_id, base_path, log_file) {

  # Internal function for creating basic barplots with a dynamic taxonomic column
  base_barplot <- function(plot_data, x_value, y_value, colorset, tax_column,
                           x_label = "Sample", y_label, facet_vars) {
    p <- ggplot(plot_data,
                aes(x = !!sym(x_value),
                    y = !!sym(y_value),
                    fill = !!sym(tax_column))) +
      geom_bar(stat = "identity") +
      scale_fill_manual(name = tax_column, values = colorset) +
      guides(fill = guide_legend(nrow = 8)) +
      theme_classic(base_size = 14) +
      labs(x = x_label,
           y = y_label,
           fill = tax_column) +
      theme(
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(),
        legend.key.size = unit(5, "pt"),
        legend.position = "bottom",
        strip.background = element_rect(colour = "white"),
        ggh4x.facet.nestline = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
        panel.spacing.x = unit(0.3, "lines"),
        strip.text = element_text(face = "bold", size = 8),
      )

    # Add nested facets conditionally if variables are provided
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
  export_folder <- file.path(base_path, "02_exports")

  barplot_folder <- file.path(figures_folder, "barplot")
  if(!dir.exists(barplot_folder)) { dir.create(barplot_folder, recursive = TRUE) }

  barplot_table_folder <- file.path(export_folder, "barplot")
  if(!dir.exists(barplot_table_folder)) { dir.create(barplot_table_folder, recursive = TRUE) }

  # Process data sequentially for each taxonomic rank
  for (tax in taxrank) {
    psdata_rmp <- psdata[[paste0("psdata_rmp_", tax)]]
    psdata_qmp <- psdata[[paste0("psdata_qmp_", tax)]]

    # Summarise abundances and calculate relative percentage
    abund_rmp <- psdata_rmp %>%
      dplyr::group_by(Sample, !!sym(tax), na_type, across(-c(Abundance, OTU, read_count))) %>%
      dplyr::summarise(abund = sum(Abundance), .groups = "drop") %>%
      dplyr::group_by(Sample) %>%
      dplyr::mutate(rel_abund = abund/sum(abund) * 100) %>%
      dplyr::ungroup()

    # Determine cutoff threshold based on top taxa
    legend_cutoff_rmp <- abund_rmp %>%
      dplyr::group_by(!!sym(tax)) %>%
      dplyr::summarise(mean_rel_abund = mean(rel_abund, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(desc(mean_rel_abund)) %>%
      dplyr::slice_head(n = ntaxa)

    # Define label for low abundance lumped category
    legend_cutoff_value_rmp <- legend_cutoff_rmp %>% pull(mean_rel_abund) %>% min()
    other_label <- glue("Other max.<{round(legend_cutoff_value_rmp, 2)}%")

    # Reorder and lump lower abundace taxa into other categroy
    plot_data_rmp <- abund_rmp %>%
      dplyr::group_by(!!sym(tax)) %>%
      dplyr::mutate(total_genus_abund = mean(rel_abund, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        !!sym(tax) := fct_reorder(!!sym(tax), total_genus_abund, .desc = TRUE),
        !!sym(tax) := fct_lump_n(!!sym(tax), n = ntaxa, w = total_genus_abund, other_level = other_label)) %>%
      dplyr::group_by(Sample, !!sym(tax), na_type, across(-c(total_genus_abund, rel_abund))) %>%
      dplyr::summarise(rel_abund = sum(rel_abund), .groups = "drop")

    all_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
    present_ranks <- intersect(all_ranks, colnames(plot_data_rmp))

    export_table_rmp <- abund_rmp %>%
      dplyr::select(all_of(present_ranks), Sample, rel_abund) %>%
      tidyr::pivot_wider(names_from = Sample, values_from = rel_abund, values_fill = 0) %>%
      dplyr::group_by(across(all_of(present_ranks))) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(mean = mean(c_across(where(is.numeric)), na.rm = TRUE)) %>%
      dplyr::arrange(desc(mean)) %>%
      dplyr::select(-mean) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, digits = 3)))

    readr::write_csv(export_table_rmp, file.path(barplot_table_folder, glue::glue("barplot_rmp_data_{tax}.csv")))

    # Clean and italicize names specifically for Genus level
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

    # Generate color palette for features
    top_names <- final_levels[final_levels != other_label]
    num_top_colors <- length(top_names)
    paired_extender <- colorRampPalette(brewer.pal(12, "Paired"))
    extended_paired_colors <- paired_extender(num_top_colors)
    my_colors <- setNames(extended_paired_colors, top_names)
    my_colors[other_label] <- "#D3D3D3"

    # Call base function to assemble heatmap plot object
    barplot_rmp <- base_barplot(plot_data_rmp_cleaned, "Sample", "rel_abund", my_colors, tax, x_label = "Sample", y_label = "Relative Abundance (%)", facet_vars) +
      ggtitle("Relative Abundance") +
      scale_y_continuous(expand = c(0, 0))

    # figures saved as png and pdf
    ggsave(filename = file.path(barplot_folder, glue::glue("barplot_rmp_{tax}.png")), plot = barplot_rmp, width = plot_width, height = 10, dpi = 600)
    ggsave(filename = file.path(barplot_folder, glue::glue("barplot_rmp_{tax}.pdf")), plot = barplot_rmp, width = plot_width, height = 10)

    # qmp
    if (!is.null(norm_method)) {
      abund_qmp <- psdata_qmp %>%
        dplyr::group_by(Sample, !!sym(tax), na_type, across(-c(Abundance, OTU, read_count))) %>%
        dplyr::summarise(norm_abund = sum(Abundance), .groups = "drop") %>%
        dplyr::ungroup()

      legend_cutoff_qmp <- abund_qmp %>%
        dplyr::group_by(!!sym(tax)) %>%
        dplyr::summarise(mean_norm_abund = mean(norm_abund, na.rm = TRUE), .groups = "drop") %>%
        dplyr::arrange(desc(mean_norm_abund)) %>%
        dplyr::slice_head(n = ntaxa)

      legend_cutoff_value_qmp <- legend_cutoff_qmp %>% pull(mean_norm_abund) %>% min()
      scale_legend_cutoff_value_qmp <- 10^floor(log10(legend_cutoff_value_qmp))
      scaled_legend_cutoff_value_qmp <- round(legend_cutoff_value_qmp / scale_legend_cutoff_value_qmp, 1)
      other_label <- glue("Other max.<{scaled_legend_cutoff_value_qmp}x10^{log10(scale_legend_cutoff_value_qmp)} cell equivalents")

      plot_data_qmp <- abund_qmp %>%
        dplyr::group_by(!!sym(tax)) %>%
        dplyr::mutate(total_genus_abund = mean(norm_abund, na.rm = TRUE)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          !!sym(tax) := fct_reorder(!!sym(tax), total_genus_abund, .desc = TRUE),
          !!sym(tax) := fct_lump_n(!!sym(tax), n = ntaxa, w = total_genus_abund, other_level = other_label)) %>%
        dplyr::group_by(Sample, !!sym(tax), na_type, across(-c(total_genus_abund, norm_abund))) %>%
        dplyr::summarise(norm_abund = sum(norm_abund), .groups = "drop")

      all_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
      present_ranks <- intersect(all_ranks, colnames(plot_data_rmp))

      export_table_qmp <- abund_qmp %>%
        dplyr::select(all_of(present_ranks), Sample, norm_abund) %>%
        tidyr::pivot_wider(names_from = Sample, values_from = norm_abund, values_fill = 0) %>%
        dplyr::group_by(across(all_of(present_ranks))) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(mean = mean(c_across(where(is.numeric)), na.rm = TRUE)) %>%
        dplyr::arrange(desc(mean)) %>%
        dplyr::select(-mean) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, digits = 3)))

      readr::write_csv(export_table_qmp, file.path(barplot_table_folder, glue::glue("barplot_qmp_data_{tax}.csv")))

      if (tax == "Genus") {
        plot_data_qmp_cleaned <- plot_data_qmp %>%
          dplyr::group_by(Sample, Genus) %>%
          dplyr::mutate(
            Genus = if_else(str_detect(Genus, "^Other"), as.character(Genus),
                            case_when(
                              str_detect(Genus, "Genus of Candidatus (\\S+) bacterium (\\S+)") ~ str_replace(Genus, "Genus of Candidatus (\\S+) bacterium (\\S+)", "Genus of *Candidatus* *\\1* (\\2)"),
                              str_detect(Genus, "Genus of Candidatus (\\S+)") ~ str_replace(Genus, "Genus of Candidatus (\\S+)", "Genus of *Candidatus* *\\1*"),
                              str_detect(Genus, "(.*)_unclassified") ~ str_replace(Genus, "(.*)_unclassified", "Unclassified *\\1*"),
                              str_detect(Genus, "Genus of") ~ str_replace(Genus, "Genus of (\\S+)", "Genus of *\\1*"),
                              str_detect(Genus, "(\\S+)\\s+(\\S+)") ~ str_replace(Genus, "(\\S+)\\s+(\\S+)", "*\\1* (*\\2*)"),
                              TRUE ~ str_replace(Genus, "^(\\S*)$", "*\\1*")))) %>%
          dplyr::ungroup()
      } else {
        plot_data_qmp_cleaned <- plot_data_qmp
      }

      genus_order <- plot_data_qmp_cleaned %>%
        dplyr::group_by(!!sym(tax)) %>%
        dplyr::summarise(avg = mean(norm_abund), .groups = "drop") %>%
        dplyr::filter(!!sym(tax) != other_label) %>%
        dplyr::arrange(desc(avg)) %>%
        dplyr::pull(!!sym(tax)) %>%
        as.character()

      final_levels <- c(genus_order, other_label)
      plot_data_qmp_cleaned <- plot_data_qmp_cleaned %>%
        dplyr::mutate(!!sym(tax) := factor(!!sym(tax), levels = final_levels))

      top_names <- final_levels[final_levels != other_label]
      num_top_colors <- length(top_names)

      paired_extender <- colorRampPalette(brewer.pal(12, "Paired"))
      extended_paired_colors <- paired_extender(num_top_colors)

      my_colors <- setNames(extended_paired_colors, top_names)
      my_colors[other_label] <- "#D3D3D3"

      barplot_qmp <- base_barplot(plot_data_qmp_cleaned, "Sample", "norm_abund", my_colors, tax, x_label = "Sample", y_label = "Cell equivalents (Cells/ml) sample", facet_vars) +
        ggtitle("Absolute Abundance") +
        scale_y_continuous(labels = function(x) {
          ifelse(x == 0, "0", sapply(x, function(num) {
            base <- floor(log10(abs(num)))
            mantissa <- num / 10^base
            ifelse(base == 0, as.character(mantissa),
                   as.expression(bquote(.(round(mantissa, 1)) ~ "×" ~ 10^.(base))))
          }))
        }, expand = c(0, 0), limits = c(0, NA))

      # figures saved as png and pdf
      ggsave(filename = file.path(barplot_folder, glue::glue("barplot_qmp_{tax}.png")), plot = barplot_qmp, width = plot_width, height = 10, dpi = 600)
      ggsave(filename = file.path(barplot_folder, glue::glue("barplot_qmp_{tax}.pdf")), plot = barplot_qmp, width = plot_width, height = 10)
    }
  }
}
