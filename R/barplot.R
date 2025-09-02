#' Create Barplots for Relative and Absolute Abundance
#'
#' This function generates barplots for microbial data at the genus level. It supports both
#' relative and absolute abundance data and can include facets based on available metadata factors.
#' The resulting plots can be saved as PDF files, and the underlying data can be exported as
#' CSV and RDS files.
#'
#' @inheritParams group_tax
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
barplot = function(physeq = rarefied_genus_psmelt,
                   ntaxa = NULL,
                   norm_method = NULL,
                   sample_matrix = NULL,
                   group_by_factor = NULL,
                   taxrank = c("Phylum", "Class", "Order", "Family", "Genus"),
                   date_factor = NULL) {

  log_message(paste("Step 12: Creating barplot.", paste(projects, collapse = ", ")), log_file)

  # Internal function for creating basic barplots with a dynamic taxonomic column
  base_barplot = function(plot_data, x_value, y_value, colorset, tax_column,
                          x_label = "Sample", y_label = "Cell equivalents (Cells/ml) sample") {
    p = ggplot(plot_data, aes(x = !!sym(x_value),
                              y = !!sym(y_value),
                              fill = !!sym(tax_column))) +
      geom_bar(stat = "identity") +
      scale_fill_manual(name = tax_column, values = colorset) +
      theme_classic(base_size = 14) +
      labs(x = x_label, y = y_label, fill = tax_column) +
      theme(axis.ticks.x = element_blank(),
            legend.text = element_markdown(),
            legend.key.size = unit(5, "pt"),
            legend.position = "bottom",
            strip.background = element_rect(colour = "white"),
            strip.text = element_text(face = "bold"),
            ggh4x.facet.nestline = element_line(colour = "black"))

    if (ntaxa > 23) {
      p = p + guides(fill = guide_legend(nrow = 14))
    } else {
      p = p + guides(fill = guide_legend(nrow = 8))
    }

    if (!is.null(present_factors)) {
      p = p + theme(axis.text.x = element_blank())
    } else {
      p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0))
    }
    return(p)
  }

  # Function to add facets if there are present factors
  facet_add = function(present_factors) {
    if (!is.null(present_factors) && length(present_factors) > 0) {
      return(
        facet_nested(
          cols = vars(!!!syms(present_factors)), scales = "free_x", space = "free",
          nest_line = element_line(linetype = 1)
        )
      )
    } else {
      return(NULL)
    }
  }

  # Global directories and project variables (assumed to be available)
  project_name = projects
  project_folder = paste0(base_path, project_name)
  figure_folder_pdf = paste0(project_folder, "/figures/PDF_figures/")
  if(!dir.exists(figure_folder_pdf)) { dir.create(figure_folder_pdf) }
  figure_folder_png = paste0(project_folder, "/figures/PNG_figures/")
  if(!dir.exists(figure_folder_png)) { dir.create(figure_folder_png) }
  destination_folder = paste0(project_folder, "/input_data/")
  output_folder_csv_files = paste0(project_folder, "/output_data/csv_files/")
  output_folder_rds_files = paste0(project_folder, "/output_data/rds_files/After_cleaning_rds_files/")

  # Define the color palette
  dark2_colors <- brewer.pal(8, "Dark2")
  paired_colors <- brewer.pal(12, "Paired")
  set_colors <- brewer.pal(8, "Set2")
  set1_colors <- brewer.pal(8, "Set1")
  spectral_colors <- brewer.pal(11, "Spectral")
  additional_palette <- unique(c(
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5",
    "#393b79", "#9c755f", "#e7298a", "#66c2a5", "#fc8d62",
    "#8da0cb", "#e78ac3", "#a6d854", "#ffed6f", "#ffff33", "#fdbf6f", "#ff7f00", "#6a3d9a", "#b15928",
    "#e41a1c", "#377eb8", "#4daf4a", "#ff6a4d",
    "#c6dbef", "#fdae61", "#fee08b", "#91bfdb", "#d73027", "#4575b4", "#313695", "#ffcc00", "#a1d99b",
    "#ff99cc", "#32cd32", "#ff6347", "#20b2aa", "#c71585", "#3cb371", "#6495ed", "#9b59b6",
    "#2ecc71", "#e74c3c", "#3498db", "#f39c12", "#8e44ad", "#16a085", "#f1c40f", "#d35400", "#27ae60",
    "#2980b9", "#e67e22"
  ))
  colorset <- unique(c(dark2_colors, paired_colors, set_colors, set1_colors, spectral_colors, additional_palette))

  if (is.null(ntaxa)) {
    ntaxa = 23
  }

  # Iterate over the desired taxonomic levels
  for (tax in taxrank) {
    log_message(paste("Processing taxonomic level:", tax), log_file)

    # Choose the appropriate dataset based on norm_method
    if (is.null(norm_method)) {
      copy_number_corrected_data = physeq[[paste0("psmelt_copy_number_corrected_", tax)]]
    } else if (norm_method == "fcm" || norm_method == "qpcr") {
      copy_number_corrected_data = physeq[[paste0("psmelt_copy_number_corrected_", tax)]]
    }

    # Create output folders (barplot folder and tax folder)
    barplot_folder_png = paste0(figure_folder_png, "Barplot/")
    if(!dir.exists(barplot_folder_png)) { dir.create(barplot_folder_png) }
    tax_folder_png = paste0(barplot_folder_png, tax, "/")
    if(!dir.exists(tax_folder_png)) { dir.create(tax_folder_png) }

    barplot_folder_pdf = paste0(figure_folder_pdf, "Barplot/")
    if(!dir.exists(barplot_folder_pdf)) { dir.create(barplot_folder_pdf) }
    tax_folder_pdf = paste0(barplot_folder_pdf, tax, "/")
    if(!dir.exists(tax_folder_pdf)) { dir.create(tax_folder_pdf) }

    variable_columns = intersect(present_variable_factors, colnames(copy_number_corrected_data))
    factor_columns = unique(c(variable_columns))
    present_factors = if (length(factor_columns) > 0) factor_columns else NULL

    # Define the complete taxonomic order and extract all higher levels than the current tax
    taxonomy_order <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
    current_tax <- tax   # current taxonomic column (e.g., "Family" or "Tax_label")
    current_index <- match(current_tax, taxonomy_order)
    higher_levels <- taxonomy_order[1:(current_index)]

    # Calculate relative abundance and modify taxonomic labels
    genus_abund_rel =
      copy_number_corrected_data %>%
      group_by(Sample, !!!syms(higher_levels), na_type, !!!syms(present_factors)) %>%
      summarise(abund = sum(Abundance), .groups = "drop") %>%
      group_by(Sample) %>%
      mutate(mean_rel_abund = abund/sum(abund) * 100) %>%
      ungroup() %>%
      group_by(Sample, !!sym(current_tax)) %>%
      mutate(!!sym(current_tax) := str_replace(!!sym(current_tax), "(.*)_unclassified", "Unclassified *\\1*")) %>%
      mutate(!!sym(current_tax) := case_when(
        str_detect(!!sym(current_tax), "Genus of") ~ str_replace(!!sym(current_tax), "Genus of (\\S+)", "Genus of *\\1*"),
        str_detect(!!sym(current_tax), "(\\S+)\\s+(\\S+)") ~ str_replace(!!sym(current_tax), "(\\S+)\\s+(\\S+)", "*\\1* (*\\2*)"),
        TRUE ~ str_replace(!!sym(current_tax), "^(\\S*)$", "*\\1*")
      ))

    # if date factor is presetn it is put in correct format and chronological order
    if (!is.null(date_factor) && date_factor %in% present_factors) {
      genus_abund_rel <- genus_abund_rel %>%
        mutate(!!sym(date_factor) := as.Date(!!sym(date_factor), format = "%d/%m/%Y")) %>%
        arrange(!!sym(date_factor))
    }

    # create a wide table with tax rank and relative abundance
    genus_abund_wide_rel =
      genus_abund_rel %>%
      select(!!!syms(higher_levels), Sample, mean_rel_abund) %>%
      pivot_wider(names_from = Sample, values_from = mean_rel_abund, values_fill = 0) %>%
      group_by_at(vars(!!!syms(higher_levels))) %>%
      rowwise() %>%
      mutate(mean = mean(c_across(where(is.numeric)))) %>%
      arrange(desc(mean)) %>%
      select(-mean) %>%
      ungroup() %>%
      mutate(across(where(is.numeric), ~round(.x, digits = 2))) %>%
      select(!!!syms(higher_levels), everything())

    # create a folder for each taxonomic rank table
    barplot_csv_folder = paste0(output_folder_csv_files, "Barplot/")
    if(!dir.exists(barplot_csv_folder)) { dir.create(barplot_csv_folder) }
    tax_csv_folder = paste0(barplot_csv_folder, tax, "/")
    if(!dir.exists(tax_csv_folder)) { dir.create(tax_csv_folder) }

    # save the wide table
    output_file_path = paste0(tax_csv_folder, project_name, "_barplot_relative_data_", tax, "_level.csv")
    write_csv(genus_abund_wide_rel, file = output_file_path, col_names = TRUE)
    log_message(paste("Relative data saved as .csv object in", output_file_path), log_file)

    # calculate the top ntaxa shown in the legend
    legend_cutoff_rel =
      genus_abund_rel %>%
      group_by(!!sym(current_tax)) %>%
      summarise(mean_rel_abund = mean(mean_rel_abund, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(mean_rel_abund)) %>%
      slice_head(n = ntaxa)

    legend_cutoff_names_rel = legend_cutoff_rel %>% pull(!!sym(current_tax))
    legend_cutoff_value_rel = legend_cutoff_rel %>% pull(mean_rel_abund) %>% min()

    genus_pool_rel =
      genus_abund_rel %>%
      group_by(!!sym(current_tax)) %>%
      summarise(pool = max(mean_rel_abund) < legend_cutoff_value_rel,
                mean = mean(mean_rel_abund), .groups = "drop")

    colorset <- unique(c(dark2_colors, paired_colors, set_colors, set1_colors, spectral_colors, additional_palette))

    # assign a color to all top ntaxa
    inner_join(genus_abund_rel, genus_pool_rel, by = current_tax) %>%
      mutate(!!sym(current_tax) := if_else(!!sym(current_tax) %in% legend_cutoff_names_rel,
                                           !!sym(current_tax),
                                           glue("Other max.<{round(legend_cutoff_value_rel, 2)}%"))) %>%
      group_by(!!sym(current_tax)) %>%
      summarise(mean_rel_abund = sum(mean_rel_abund),
                median = median(mean), .groups = "drop") %>%
      mutate(!!sym(current_tax) := factor(!!sym(current_tax)),
             !!sym(current_tax) := fct_reorder(!!sym(current_tax), median, .desc = TRUE)) %>%
      arrange(desc(median)) %>%
      pull(!!sym(current_tax)) %>%
      unique() -> names(colorset)

    colorset[glue("Other max.<{round(legend_cutoff_value_rel, 2)}%")] <- "#D3D3D3"

    # creating the plot data for relative barplots
    plot_data_rel =
      inner_join(genus_abund_rel, genus_pool_rel, by = current_tax) %>%
      mutate(!!sym(current_tax) := if_else(!!sym(current_tax) %in% legend_cutoff_names_rel,
                                           !!sym(current_tax),
                                           glue("Other max.<{round(legend_cutoff_value_rel, 2)}%"))) %>%
      group_by(Sample, !!!syms(higher_levels), na_type, !!!syms(present_factors)) %>%
      summarise(mean_rel_abund = sum(mean_rel_abund),
                median = median(mean), .groups = "drop") %>%
      mutate(!!sym(current_tax) := factor(!!sym(current_tax)),
             !!sym(current_tax) := fct_reorder(!!sym(current_tax), median, .desc = TRUE))

    tax_rds_folder = paste0(output_folder_rds_files, tax, "/")
    if(!dir.exists(tax_rds_folder)) { dir.create(tax_rds_folder) }

    output_file_path = paste0(tax_rds_folder, project_name, "_pstibble_relative_data.rds")
    saveRDS(plot_data_rel, file = output_file_path)
    log_message(paste("Relative pstibble data saved as .rds object in", output_file_path), log_file)

    na_types = unique(plot_data_rel$na_type)

    # creating relative barplots for number of na_types
    if (length(na_types) == 1) {
      barplot_relative =
        base_barplot(plot_data_rel, "Sample", "mean_rel_abund", colorset, current_tax,
                     x_label = "Sample", y_label = "Relative Abundance (%)") +
        ggtitle("Relative Abundance") +
        facet_add(present_factors) +
        scale_y_continuous(expand = c(0, 0))

      n_samples <- length(unique(plot_data_rel$Sample))
      fig.width <- max(12, n_samples * 0.40)

      figure_file_path = paste0(tax_folder_png, project_name, "_barplot_relative_", tax, "_level.png")
      ggsave(filename = figure_file_path, plot = barplot_relative, width = fig.width, height = 10, dpi = 300)
      log_message(paste("Relative barplot saved as .png object in", figure_file_path), log_file)

      figure_file_path = paste0(tax_folder_pdf, project_name, "_barplot_relative_", tax, "_level.pdf")
      ggsave(filename = figure_file_path, plot = barplot_relative, width = fig.width, height = 10)
      log_message(paste("Relative barplot saved as .pdf object in", figure_file_path), log_file)

    } else if (length(na_types) == 2) {
      plot_data_rel_dna =
        plot_data_rel %>% filter(na_type == "dna")

      barplot_relative_dna =
        base_barplot(plot_data_rel_dna, "Sample", "mean_rel_abund", colorset, current_tax,
                     x_label = "Sample", y_label = "Relative Abundance (%)") +
        ggtitle("Relative Abundance - DNA") +
        facet_add(present_factors) +
        scale_y_continuous(expand = c(0, 0))

      plot_data_rel_rna =
        plot_data_rel %>% filter(na_type == "rna")

      barplot_relative_rna =
        base_barplot(plot_data_rel_rna, "Sample", "mean_rel_abund", colorset, current_tax,
                     x_label = "Sample", y_label = "Relative Abundance (%)") +
        ggtitle("Relative Abundance - RNA") +
        facet_add(present_factors) +
        scale_y_continuous(expand = c(0, 0))

      n_samples_dna <- length(unique(plot_data_rel_dna$Sample))
      fig.width_dna <- max(12, n_samples_dna * 0.40)

      n_samples_rna <- length(unique(plot_data_rel_rna$Sample))
      fig.width_rna <- max(12, n_samples_rna * 0.40)

      # save figure as png
      figure_file_path = paste0(tax_folder_png, project_name, "_barplot_relative_dna_", tax, "_level.png")
      ggsave(filename = figure_file_path, plot = barplot_relative_dna, width = fig.width_dna, height = 10, dpi = 600)
      log_message(paste("Relative barplot saved as .png object in", figure_file_path), log_file)

      figure_file_path = paste0(tax_folder_png, project_name, "_barplot_relative_rna_", tax, "_level.png")
      ggsave(filename = figure_file_path, plot = barplot_relative_rna, width = fig.width_rna, height = 10, dpi = 600)
      log_message(paste("Relative barplot saved as .png object in", figure_file_path), log_file)

      # save figure as pdf
      figure_file_path = paste0(tax_folder_pdf, project_name, "_barplot_relative_dna_", tax, "_level.pdf")
      ggsave(filename = figure_file_path, plot = barplot_relative_dna, width = fig.width_dna, height = 10)
      log_message(paste("Relative barplot saved as .pdf object in", figure_file_path), log_file)

      figure_file_path = paste0(tax_folder_pdf, project_name, "_barplot_relative_rna_", tax, "_level.pdf")
      ggsave(filename = figure_file_path, plot = barplot_relative_rna, width = fig.width_rna, height = 10)
      log_message(paste("Relative barplot saved as .pdf object in", figure_file_path), log_file)

    }

    # Absolute data plots
    # if normalized data is present, an absolute barplot is created
    if (is.null(norm_method)) {
      log_message("No absolute data provided; skipping absolute barplot creation", log_file)
    } else if (norm_method == "fcm" || norm_method == "qpcr") {
      if (norm_method == "fcm") {
        absolute_data = physeq[[paste0("psmelt_fcm_norm_rarefied_", tax)]]
      } else if (norm_method == "qpcr") {
        absolute_data = physeq[[paste0("psmelt_qpcr_norm_rarefied_", tax)]]
      }

      taxonomy_order <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
      current_tax <- tax
      current_index <- match(current_tax, taxonomy_order)
      higher_levels <- taxonomy_order[1:(current_index)]

      # modify taxonomic labels
      genus_abund_norm =
        absolute_data %>%
        group_by(Sample, !!!syms(higher_levels), na_type, !!!syms(present_factors)) %>%
        summarise(norm_abund = sum(Abundance), .groups = "drop") %>%
        group_by(Sample, !!sym(current_tax)) %>%
        mutate(!!sym(current_tax) := str_replace(!!sym(current_tax), "(.*)_unclassified", "Unclassified *\\1*")) %>%
        mutate(!!sym(current_tax) := case_when(
          str_detect(!!sym(current_tax), "Genus of") ~ str_replace(!!sym(current_tax), "Genus of (\\S+)", "Genus of *\\1*"),
          str_detect(!!sym(current_tax), "(\\S+)\\s+(\\S+)") ~ str_replace(!!sym(current_tax), "(\\S+)\\s+(\\S+)", "*\\1* (*\\2*)"),
          TRUE ~ str_replace(!!sym(current_tax), "^(\\S*)$", "*\\1*")
        ))

      # create a wide table with tax rank and absolute abundance
      genus_abund_wide_norm =
        genus_abund_norm %>%
        select(!!!syms(higher_levels), Sample, norm_abund) %>%
        pivot_wider(names_from = Sample, values_from = norm_abund, values_fill = 0) %>%
        group_by_at(vars(!!!syms(higher_levels))) %>%
        rowwise() %>%
        mutate(mean = mean(c_across(where(is.numeric)))) %>%
        arrange(desc(mean)) %>%
        select(-mean) %>%
        ungroup() %>%
        mutate(across(where(is.numeric), ~round(.x, digits = 2))) %>%
        select(!!!syms(higher_levels), everything())

      # create a folder for each taxonomic rank table
      barplot_csv_folder = paste0(output_folder_csv_files, "Barplot/")
      if(!dir.exists(barplot_csv_folder)) { dir.create(barplot_csv_folder) }
      tax_csv_folder = paste0(barplot_csv_folder, tax, "/")
      if(!dir.exists(tax_csv_folder)) { dir.create(tax_csv_folder) }

      # save the wide table
      output_file_path = paste0(tax_csv_folder, project_name, "_barplot_absolute_data_", tax, "_level.csv")
      write_csv(genus_abund_wide_norm, file = output_file_path, col_names = TRUE)
      log_message(paste("Absolute data saved as .csv object in", output_file_path), log_file)

      # calculate the top ntaxa shown in the legend
      legend_cutoff_norm =
        genus_abund_norm %>%
        group_by(!!sym(current_tax)) %>%
        summarise(mean_norm_abund = mean(norm_abund, na.rm = TRUE), .groups = "drop") %>%
        arrange(desc(mean_norm_abund)) %>%
        slice_head(n = ntaxa)

      legend_cutoff_names_norm = legend_cutoff_norm %>% pull(!!sym(current_tax))
      legend_cutoff_value_norm = legend_cutoff_norm %>% pull(mean_norm_abund) %>% min()

      scale_legend_cutoff_value_norm = 10^floor(log10(legend_cutoff_value_norm))
      scaled_legend_cutoff_value_norm = round(legend_cutoff_value_norm / scale_legend_cutoff_value_norm, 1)

      # assign a color to all top ntaxa
      genus_pool_norm =
        genus_abund_norm %>%
        group_by(!!sym(current_tax)) %>%
        summarise(pool = max(norm_abund) < legend_cutoff_value_norm,
                  mean = mean(norm_abund), .groups = "drop")

      norm_taxa = inner_join(genus_abund_norm, genus_pool_norm, by = current_tax) %>%
        mutate(!!sym(current_tax) := if_else(!!sym(current_tax) %in% legend_cutoff_names_norm,
                                             !!sym(current_tax),
                                             glue("Other max.<{scaled_legend_cutoff_value_norm}x10^{log10(scale_legend_cutoff_value_norm)}\ncell equivalents"))) %>%
        group_by(Sample, !!!syms(higher_levels), na_type, !!!syms(present_factors)) %>%
        summarise(norm_abund = sum(norm_abund),
                  median = median(mean), .groups = "drop") %>%
        mutate(!!sym(current_tax) := factor(!!sym(current_tax)),
               !!sym(current_tax) := fct_reorder(!!sym(current_tax), median, .desc = TRUE)) %>%
        arrange(desc(median)) %>%
        pull(!!sym(current_tax)) %>%  # use the dynamic column
        unique()

      new_taxa <- setdiff(norm_taxa, names(colorset))
      num_new_taxa = length(new_taxa)
      if (num_new_taxa > 0) {
      assigned_taxa <- names(colorset)
      available_colors <- setdiff(scales::hue_pal()(num_new_taxa + length(colorset)), colorset)
      colorset[new_taxa] <- available_colors[1:num_new_taxa]
      names(colorset)[(length(colorset) - num_new_taxa + 1):length(colorset)] <- new_taxa
      }
      colorset[glue("Other max.<{scaled_legend_cutoff_value_norm}x10^{log10(scale_legend_cutoff_value_norm)}\ncell equivalents")] <- "#D3D3D3"

      # creating the plot data for absolute barplots
      plot_data_norm =
        inner_join(genus_abund_norm, genus_pool_norm, by = current_tax) %>%
        mutate(!!sym(current_tax) := if_else(!!sym(current_tax) %in% legend_cutoff_names_norm,
                                             !!sym(current_tax),
                                             glue("Other max.<{scaled_legend_cutoff_value_norm}x10^{log10(scale_legend_cutoff_value_norm)}\ncell equivalents"))) %>%
        group_by(Sample, !!!syms(higher_levels), na_type, !!!syms(present_factors)) %>%
        summarise(norm_abund = sum(norm_abund),
                  median = median(mean), .groups = "drop") %>%
        mutate(!!sym(current_tax) := factor(!!sym(current_tax)),
               !!sym(current_tax) := fct_reorder(!!sym(current_tax), median, .desc = TRUE))

      tax_rds_folder = paste0(output_folder_rds_files, tax, "/")
      if(!dir.exists(tax_rds_folder)) { dir.create(tax_rds_folder) }

      output_file_path = paste0(tax_rds_folder, project_name, "_pstibble_absolute_data_", tax, "_level.rds")
      saveRDS(plot_data_norm, file = output_file_path)
      log_message(paste("Absolute pstibble data saved as .rds object in", output_file_path), log_file)

      barplot_folder_png = paste0(figure_folder_png, "Barplot/")
      if(!dir.exists(barplot_folder_png)) { dir.create(barplot_folder_png) }
      tax_folder_png = paste0(barplot_folder_png, tax, "/")
      if(!dir.exists(tax_folder_png)) { dir.create(tax_folder_png) }

      barplot_folder_pdf = paste0(figure_folder_pdf, "Barplot/")
      if(!dir.exists(barplot_folder_pdf)) { dir.create(barplot_folder_pdf) }
      tax_folder_pdf = paste0(barplot_folder_pdf, tax, "/")
      if(!dir.exists(tax_folder_pdf)) { dir.create(tax_folder_pdf) }

      na_types = unique(plot_data_norm$na_type)

      if (length(na_types) == 1) {
        # choice of `liquid` for cells per ml annotation on the y as of the barplot or `solid` for cells per gram
        if (sample_matrix == "liquid") {
          ylabel = "Cell equivalents (Cells/ml) sample"
        } else if (sample_matrix == "solid") {
          ylabel = "Cell equivalents (Cells/gram) sample"
        }

        if (!is.null(group_by_factor)) {
          all_plots = list()
          factors = unique(plot_data_norm[[group_by_factor]])
          for (x in factors) {
            data_filtered = plot_data_norm %>% filter(.data[[group_by_factor]] == x)
            plot <- base_barplot(data_filtered, "Sample", "norm_abund", colorset, current_tax,
                                 x_label = "Sample", y_label = ylabel) +
              ggtitle("Absolute Abundance") +
              facet_add(present_factors) +
              scale_y_continuous(labels = function(x) {
                ifelse(x == 0, "0", sapply(x, function(num) {
                  base <- floor(log10(abs(num)))
                  mantissa <- num / 10^base
                  ifelse(base == 0, as.character(mantissa),
                         as.expression(bquote(.(round(mantissa, 1)) ~ "×" ~ 10^.(base))))
                }))
              }, expand = c(0, 0), limits = c(0, NA)) +
              theme(legend.position = "none")
            all_plots[[x]] <- plot
          }
          legend = get_legend(all_plots[[1]] + theme(legend.position = "right"))
          nplots = length(all_plots)
          final_plot = wrap_plots(all_plots, ncol = nplots) +
            plot_annotation(title = "Absolute Abundance - DNA")
          barplot_absolute = plot_grid(final_plot, legend, ncol = 1, rel_heights = c(3, 1))
        } else {
          barplot_absolute =
            base_barplot(plot_data_norm, "Sample", "norm_abund", colorset, current_tax,
                         x_label = "Sample", y_label = ylabel) +
            ggtitle("Absolute Abundance") +
            facet_add(present_factors) +
            scale_y_continuous(labels = function(x) {
              ifelse(x == 0, "0", sapply(x, function(num) {
                base <- floor(log10(abs(num)))
                mantissa <- num / 10^base
                ifelse(base == 0, as.character(mantissa),
                       as.expression(bquote(.(round(mantissa, 1)) ~ "×" ~ 10^.(base))))
              }))
            }, expand = c(0, 0), limits = c(0, NA))
        }

        n_samples <- length(unique(plot_data_norm$Sample))
        fig.width <- max(12, n_samples * 0.40)

        # figures saved as png
        figure_file_path = paste0(tax_folder_png, project_name, "_barplot_absolute_", tax, "_level.png")
        ggsave(filename = figure_file_path, plot = barplot_absolute, width = fig.width, height = 10, dpi = 600)
        log_message(paste("Absolute barplot saved as .png object in", figure_file_path), log_file)

        # figure saved as pdf
        figure_file_path = paste0(tax_folder_pdf, project_name, "_barplot_absolute_", tax, "_level.pdf")
        ggsave(filename = figure_file_path, plot = barplot_absolute, width = fig.width, height = 10)
        log_message(paste("Absolute barplot saved as .pdf object in", figure_file_path), log_file)

      } else if (length(na_types) == 2) {
        if (sample_matrix == "liquid") {
          ylabel = "Cell equivalents (Cells/ml) sample"
        } else if (sample_matrix == "solid") {
          ylabel = "Cell equivalents (Cells/gram) sample"
        }
        # DNA
        plot_data_norm_dna =
          plot_data_norm %>%
          filter(na_type == "dna")
        if (!is.null(group_by_factor)) {
          all_plots = list()
          factors = unique(plot_data_norm_dna[[group_by_factor]])
          for (x in factors) {
            data_filtered = plot_data_norm_dna %>% filter(.data[[group_by_factor]] == x)
            plot <- base_barplot(data_filtered, "Sample", "norm_abund", colorset, current_tax,
                                 x_label = "Sample", y_label = ylabel) +
              ggtitle("Absolute Abundance - DNA") +
              facet_add(present_factors) +
              scale_y_continuous(labels = function(x) {
                ifelse(x == 0, "0", sapply(x, function(num) {
                  base <- floor(log10(abs(num)))
                  mantissa <- num / 10^base
                  ifelse(base == 0, as.character(mantissa),
                         as.expression(bquote(.(round(mantissa, 1)) ~ "×" ~ 10^.(base))))
                }))
              }, expand = c(0, 0), limits = c(0, NA)) +
              theme(legend.position = "none")
            all_plots[[x]] <- plot
          }
          legend = get_legend(all_plots[[1]] + theme(legend.position = "right"))
          nplots = length(all_plots)
          final_plot = wrap_plots(all_plots, ncol = nplots) +
            plot_annotation(title = "Absolute Abundance - DNA")
          barplot_absolute_dna = plot_grid(final_plot, legend, ncol = 1, rel_heights = c(3, 1))
        } else {
          barplot_absolute_dna =
            base_barplot(plot_data_norm_dna, "Sample", "norm_abund", colorset, current_tax,
                         x_label = "Sample", y_label = ylabel) +
            ggtitle("Absolute Abundance - DNA") +
            facet_add(present_factors) +
            scale_y_continuous(labels = function(x) {
              ifelse(x == 0, "0", sapply(x, function(num) {
                base <- floor(log10(abs(num)))
                mantissa <- num / 10^base
                ifelse(base == 0, as.character(mantissa),
                       as.expression(bquote(.(round(mantissa, 1)) ~ "×" ~ 10^.(base))))
              }))
            }, expand = c(0, 0), limits = c(0, NA))
        }
        # RNA
        plot_data_norm_rna =
          plot_data_norm %>%
          filter(na_type == "rna")
        if (!is.null(group_by_factor)) {
          all_plots = list()
          factors = unique(plot_data_norm_rna[[group_by_factor]])
          for (x in factors) {
            data_filtered = plot_data_norm_rna %>% filter(.data[[group_by_factor]] == x)
            plot <- base_barplot(data_filtered, "Sample", "norm_abund", colorset, current_tax,
                                 x_label = "Sample", y_label = ylabel) +
              ggtitle("Absolute Abundance - RNA") +
              facet_add(present_factors) +
              scale_y_continuous(labels = function(x) {
                ifelse(x == 0, "0", sapply(x, function(num) {
                  base <- floor(log10(abs(num)))
                  mantissa <- num / 10^base
                  ifelse(base == 0, as.character(mantissa),
                         as.expression(bquote(.(round(mantissa, 1)) ~ "×" ~ 10^.(base))))
                }))
              }, expand = c(0, 0), limits = c(0, NA)) +
              theme(legend.position = "none")
            all_plots[[x]] <- plot
          }
          legend = get_legend(all_plots[[1]] + theme(legend.position = "right"))
          nplots = length(all_plots)
          final_plot = wrap_plots(all_plots, ncol = nplots) +
            plot_annotation(title = "Absolute Abundance - DNA")
          barplot_absolute_rna = plot_grid(final_plot, legend, ncol = 1, rel_heights = c(3, 1))
        } else {
          barplot_absolute_rna =
            base_barplot(plot_data_norm_rna, "Sample", "norm_abund", colorset, current_tax,
                         x_label = "Sample", y_label = ylabel) +
            ggtitle("Absolute Abundance - RNA") +
            facet_add(present_factors) +
            scale_y_continuous(labels = function(x) {
              ifelse(x == 0, "0", sapply(x, function(num) {
                base <- floor(log10(abs(num)))
                mantissa <- num / 10^base
                ifelse(base == 0, as.character(mantissa),
                       as.expression(bquote(.(round(mantissa, 1)) ~ "×" ~ 10^.(base))))
              }))
            }, expand = c(0, 0), limits = c(0, NA))
        }

        n_samples_dna <- length(unique(plot_data_norm_dna$Sample))
        fig.width_dna <- max(12, n_samples_dna * 0.40)

        n_samples_rna <- length(unique(plot_data_norm_rna$Sample))
        fig.width_rna <- max(12, n_samples_rna * 0.40)

        # figures saved as png
        figure_file_path = paste0(tax_folder_png, project_name, "_barplot_absolute_dna_", tax, "_level.png")
        ggsave(filename = figure_file_path, plot = barplot_absolute_dna, width = fig.width_dna, height = 10, dpi = 600)
        log_message(paste("Absolute barplot saved as .png object in", figure_file_path), log_file)

        figure_file_path = paste0(tax_folder_png, project_name, "_barplot_absolute_rna_", tax, "_level.png")
        ggsave(filename = figure_file_path, plot = barplot_absolute_rna, width = fig.width_rna, height = 8, dpi = 600)
        log_message(paste("Absolute barplot saved as .png object in", figure_file_path), log_file)

        # figures saved as pdf
        figure_file_path = paste0(tax_folder_pdf, project_name, "_barplot_absolute_dna_", tax, "_level.pdf")
        ggsave(filename = figure_file_path, plot = barplot_absolute_dna, width = fig.width_dna, height = 8)
        log_message(paste("Absolute barplot saved as .pdf object in", figure_file_path), log_file)

        figure_file_path = paste0(tax_folder_pdf, project_name, "_barplot_absolute_rna_", tax, "_level.pdf")
        ggsave(filename = figure_file_path, plot = barplot_absolute_rna, width = fig.width_rna, height = 8)
        log_message(paste("Absolute barplot saved as .pdf object in", figure_file_path), log_file)
      }
    }
  }

  log_message("Barplot successfully plotted.", log_file)
}
