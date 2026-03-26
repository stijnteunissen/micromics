### heatmap clr ####

heatmap_clr = function(physeq = rarefied_genus_psmelt,
                       ntaxa = NULL,
                       norm_method = NULL,
                       taxrank = c("Phylum", "Class", "Order", "Family", "Genus"),
                       date_factor = NULL,
                       SampleID_xas = FALSE) {

  log_message(paste("Step 15: Creating CLR heatmap.", paste(projects, collapse = ", ")), log_file)

  base_heatmap = function(plot_data, x_value, clr_value, legend_name, x_label = "Sample", tax_column) {
    p <- ggplot(plot_data, aes(x = Sample,
                               y = !!sym(tax_column))) +
      geom_tile(aes(fill = !!sym(clr_value)), color = NA) +
      #scale_fill_gradient(low = "white", high = "darkred", name = legend_name) +
      scale_fill_gradient2(low = "darkgreen", mid = "white", high = "darkred", midpoint = 0, name = "CLR", limits = c(-10, 40)) +
      labs(x = x_label,
           y = tax_column) +
      theme_classic() +
      theme(axis.text.y = element_markdown(size = 10),
            axis.ticks.x = element_blank(),
            strip.placement = "outside",
            strip.background = element_blank(),
            #strip.text = element_text(face = "bold", angle = 90, vjust = 0.5, hjust = 0),
            strip.text = element_text(face = "bold"),
            ggh4x.facet.nestline = element_line(colour = "black")) +
      scale_x_discrete(expand = c(0, 0)) +
      geom_text(aes(label = ifelse(clr_value > 1, paste0(round(clr_value, 1)), ifelse(clr_value == 0, ".", "")),
                    color = ifelse(clr_value > 50, "#D3D3D3", "black")),
                size = 3)

    if (!is.null(present_factors) && !isTRUE(SampleID_xas)) {
      p = p + theme(axis.text.x = element_blank())
    } else {
      p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0))
    }
    return(p)
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
  output_folder_csv_files = paste0(project_folder, "/output_data/csv_files/")
  if(!dir.exists(output_folder_csv_files)) { dir.create(output_folder_csv_files, recursive = TRUE) }

  if (is.null(ntaxa)) {
    ntaxa = 23
  }

  for (tax in taxrank) {
    log_message(paste("Processing taxonomic level:", tax), log_file)

    if (is.null(norm_method)) {
      copy_number_corrected_data <- physeq[[paste0("psdata_copy_number_corrected_", tax)]]
    } else if (norm_method == "fcm" || norm_method == "qpcr") {
      copy_number_corrected_data <- physeq[[paste0("psdata_copy_number_corrected_", tax)]]
    }

    copy_number_corrected_data = physeq_IN_list$psdata_copy_number_corrected_Genus
    tax = "Genus"

    # clr transform
    copy_number_corrected_data_clr <- microbiome::transform(copy_number_corrected_data, "clr") %>%
      psmelt() %>%
      as_tibble()

    # relative abundance for selection taxa
    copy_number_corrected_data_rel <- copy_number_corrected_data %>%
      psmelt() %>%
      as_tibble()

    # Create output folders (barplot folder and tax folder)
    heatmap_folder_png = paste0(figure_folder_png, "Heatmap/")
    if(!dir.exists(heatmap_folder_png)) { dir.create(heatmap_folder_png) }
    tax_folder_png = paste0(heatmap_folder_png, tax, "/")
    if(!dir.exists(tax_folder_png)) { dir.create(tax_folder_png) }

    heatmap_folder_pdf = paste0(figure_folder_pdf, "Heatmap/")
    if(!dir.exists(heatmap_folder_pdf)) { dir.create(heatmap_folder_pdf) }
    tax_folder_pdf = paste0(heatmap_folder_pdf, tax, "/")
    if(!dir.exists(tax_folder_pdf)) { dir.create(tax_folder_pdf) }

    variable_columns = intersect(present_variable_factors, colnames(copy_number_corrected_data_clr))
    factor_columns = unique(c(variable_columns))
    present_factors = if (length(factor_columns) > 0) factor_columns else NULL

    # taxa selection
    # mean relative abundace
    taxa_mean <- copy_number_corrected_data_rel %>%
      group_by(Sample, !!sym(tax), na_type, !!!syms(present_factors)) %>%
      summarise(abund = sum(Abundance), .groups = "drop") %>%
      group_by(Sample) %>%
      mutate(mean_rel_abund = abund/sum(abund) * 100) %>%
      ungroup() %>%
      group_by(!!sym(tax)) %>%
      summarise(mean_abund_sample = mean(mean_rel_abund, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(mean_abund_sample)) %>%
      slice_head(n = ntaxa)

    top_taxa_mean <- taxa_mean[[tax]]

    # Define the complete taxonomic order and extract all higher levels than the current tax
    taxonomy_order <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
    current_tax <- tax   # current taxonomic column (e.g., "Family" or "Tax_label")
    current_index <- match(current_tax, taxonomy_order)
    higher_levels <- taxonomy_order[1:(current_index)]

    # CLR data
    genus_abund_clr <- copy_number_corrected_data_clr %>%
      group_by(Sample, !!!syms(higher_levels), na_type, !!!syms(present_factors)) %>%
      summarise(clr_value = sum(Abundance), .groups = "drop") %>%
      ungroup()

    # filter taxa
    plot_data_clr_mean <- genus_abund_clr %>%
      filter((!!sym(tax)) %in% top_taxa_mean) %>%
      mutate(!!sym(tax) := factor(!!sym(tax),
                                  levels = taxa_mean %>%
                                    arrange(desc(mean_abund_sample)) %>%
                                    pull(!!sym(tax))))

    plot_data_clr_mean <- plot_data_clr_mean %>%
      group_by(Sample, !!sym(tax)) %>%
      mutate(!!sym(tax) := str_replace(!!sym(tax), "(.*)_unclassified", "Unclassified *\\1*")) %>%
      mutate(!!sym(tax) := case_when(
        str_detect(!!sym(tax), "Genus of") ~ str_replace(!!sym(tax), "Genus of (\\S+)", "Genus of *\\1*"),
        str_detect(!!sym(tax), "(\\S+)\\s+(\\S+)") ~ str_replace(!!sym(tax), "(\\S+)\\s+(\\S+)", "*\\1* (*\\2*)"),
        TRUE ~ str_replace(!!sym(tax), "^(\\S*)$", "*\\1*")
      ))

    if (!is.null(date_factor) && date_factor %in% present_factors) {
      genus_abund_clr <- genus_abund_clr %>%
        mutate(!!sym(date_factor) := as.Date(!!sym(date_factor), format = "%d/%m/%Y")) %>%
        arrange(!!sym(date_factor))
    }

    genus_abund_wide_clr <- genus_abund_clr %>%
      group_by(Sample, !!sym(current_tax)) %>%
      mutate(!!sym(current_tax) := str_replace(!!sym(current_tax), "(.*)_unclassified", "Unclassified \\1")) %>%
      mutate(!!sym(current_tax) := case_when(
        str_detect(!!sym(current_tax), "Genus of") ~ str_replace(!!sym(current_tax), "Genus of (\\S+)", "Genus of \\1"),
        str_detect(!!sym(current_tax), "(\\S+)\\s+(\\S+)") ~ str_replace(!!sym(current_tax), "(\\S+)\\s+(\\S+)", "\\1 (\\2)"),
        TRUE ~ str_replace(!!sym(current_tax), "^(\\S*)$", "\\1")
      )) %>%
      select(!!!syms(higher_levels), Sample, clr_value) %>%
      pivot_wider(names_from = Sample, values_from = clr_value, values_fill = 0) %>%
      group_by_at(vars(!!!syms(higher_levels))) %>%
      rowwise() %>%
      mutate(mean = mean(c_across(where(is.numeric)))) %>%
      arrange(desc(mean)) %>%
      select(-mean) %>%
      ungroup() %>%
      mutate(across(where(is.numeric), ~round(.x, digits = 2))) %>%
      select(!!!syms(higher_levels), everything())

    # create a folder for each taxonomic rank table
    heatmap_csv_folder = paste0(output_folder_csv_files, "Heatmap/")
    if(!dir.exists(heatmap_csv_folder)) { dir.create(heatmap_csv_folder) }
    tax_csv_folder = paste0(heatmap_csv_folder, tax, "/")
    if(!dir.exists(tax_csv_folder)) { dir.create(tax_csv_folder) }

    # save the wide table
    output_file_path = paste0(tax_csv_folder, project_name, "_heatmap_clr_data_", tax, "_level.csv")
    write_csv(genus_abund_wide_clr, file = output_file_path, col_names = TRUE)
    log_message(paste("CLR data saved as .csv object in", output_file_path), log_file)

    na_types = unique(plot_data_clr_mean$na_type)

    if (length(na_types) == 1) {
      heatmap_clr_mean =
        base_heatmap(plot_data_clr_mean, "Sample", "clr_value", legend_name = "CLR", x_label = "Sample", tax) +
        scale_color_identity() +
        facet_add(present_factors)

      n_samples <- length(unique(plot_data_clr_mean$Sample))
      fig.width <- max(12, n_samples * 0.40)

      # save mean
      figure_file_path = paste0(tax_folder_png, project_name, "_heatmap_CLR_mean_", tax, "_level.png")
      ggsave(filename = figure_file_path, plot = heatmap_clr_mean, width = fig.width, height = 10, limitsize = FALSE, dpi = 600)
      log_message(paste("Relative heatmap saved as .png object in", figure_file_path), log_file)

      figure_file_path = paste0(tax_folder_pdf, project_name, "_heatmap_CLR_mean_", tax, "_level.pdf")
      ggsave(filename = figure_file_path, plot = heatmap_clr_mean, width = fig.width, height = 10, limitsize = FALSE)
      log_message(paste("Relative heatmap saved as .pdf object in", figure_file_path), log_file)

    } else if (length(na_types) == 2) {
      plot_data_dna = plot_data_clr_mean %>% filter(na_type == "dna")

      heatmap_clr_dna =
        base_heatmap(plot_data_dna, "Sample", "clr_value", legend_name = "CLR", x_label = "Sample", tax) +
        scale_color_identity() +
        facet_add(present_factors, include_na_type = TRUE) +
        theme(legend.position = "right")

      plot_data_rna = plot_data_clr_mean %>% filter(na_type == "rna")

      heatmap_clr_rna =
        base_heatmap(plot_data_rna, "Sample", "clr_value", legend_name = "CLR", x_label = "Sample", tax) +
        scale_color_identity() +
        facet_add(present_factors, include_na_type = TRUE) +
        theme(legend.position = "right")

      num_samples_dna = plot_data_dna %>% distinct(Sample) %>% nrow()
      num_samples_rna = plot_data_rna %>% distinct(Sample) %>% nrow()

      dna_width = num_samples_dna
      rna_width = num_samples_rna + 1
      total_width = dna_width + rna_width

      rel_widths = c(dna_width / total_width, rna_width / total_width)

      legend = get_legend(heatmap_clr_dna + theme(legend.position = "right"))

      combined_heatmap_rel = plot_grid(heatmap_clr_dna, heatmap_clr_rna,
                                       ncol = 2, labels = c("A", "B"), rel_widths = rel_widths)

      heatmap_relative = plot_grid(combined_heatmap_rel, legend, ncol = 2, rel_widths = c(3, 0.5))

      n_samples_dna <- length(unique(plot_data_dna$Sample))
      fig.width_dna <- max(12, n_samples_dna * 0.40)

      n_samples_rna <- length(unique(plot_data_rna$Sample))
      fig.width_rna <- max(12, n_samples_rna * 0.40)

      figure_file_path = paste0(tax_folder_png, project_name, "_heatmap_clr_dna_", tax, "_level.png")
      ggsave(filename = figure_file_path, plot = heatmap_clr_dna, width = fig.width_dna, height = 10, limitsize = FALSE, dpi = 600)
      log_message(paste("CLR heatmap saved as .png object in", figure_file_path), log_file)

      figure_file_path = paste0(tax_folder_pdf, project_name, "_heatmap_clr_dna_", tax, "_level.pdf")
      ggsave(filename = figure_file_path, plot = heatmap_clr_dna, width = fig.width_dna, height = 10, limitsize = FALSE)
      log_message(paste("CLR heatmap saved as .pdf object in", figure_file_path), log_file)

      figure_file_path = paste0(tax_folder_png, project_name, "_heatmap_clr_rna_", tax, "_level.png")
      ggsave(filename = figure_file_path, plot = heatmap_clr_rna, width = fig.width_rna, height = 10, limitsize = FALSE, dpi = 600)
      log_message(paste("CLR heatmap saved as .png object in", figure_file_path), log_file)

      figure_file_path = paste0(tax_folder_pdf, project_name, "_heatmap_clr_rna_", tax, "_level.pdf")
      ggsave(filename = figure_file_path, plot = heatmap_clr_rna, width = fig.width_rna, height = 10, limitsize = FALSE)
      log_message(paste("CLR heatmap saved as .pdf object in", figure_file_path), log_file)
    }
  }
  log_message("Heatmap CLR successfully plotted.", log_file)
}

