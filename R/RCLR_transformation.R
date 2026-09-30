RCLR_transformation <- function(physeq, taxrank, facet_var, color_var, shape_var, time_var, project_id, base_path, log_file) {

  # Setup folder paths
  clean_rds_folder <- file.path(base_path, "01_r_objects/clean")
  figures_folder <- file.path(base_path, "03_figures")
  rclr_folder <- file.path(figures_folder, "RCLR_plots")

  if (!dir.exists(rclr_folder)) {
    dir.create(rclr_folder, recursive = TRUE)
  }

  beta_div_folder <- file.path(rclr_folder, "beta_diversity")
  if (!dir.exists(beta_div_folder)) {dir.create(beta_div_folder, recursive = TRUE)}

  rclr_load_folder <- file.path(rclr_folder, "RCLR_loadings")
  if (!dir.exists(rclr_load_folder)) {dir.create(rclr_load_folder, recursive = TRUE)}

  # loop through each specified taxonomic rank
  for (tax in taxrank) {
    physeq_rmp <- physeq$physeq_rmp_rarefied

    # Agglomerate taxa if not at ASV level
    if (tax != "ASV") {
      physeq_rmp_glom <- phyloseq::tax_glom(physeq_rmp, taxrank = tax)
    } else {
      physeq_rmp_glom <- physeq_rmp
    }

    # Extract OTU table and convert to matrix
    otu_matrix <- as(otu_table(physeq_rmp_glom), "matrix")

    # Ensure samples are rows, taxa are columns
    if (taxa_are_rows(physeq_rmp_glom)) {
      otu_matrix <- t(otu_matrix)
    }

    # perform RCLR transformation and PCA
    x_rclr <- vegan::decostand(otu_matrix, method = "rclr")
    pca_rclr <- prcomp(x_rclr, center = TRUE, scale. = FALSE)

    # Save transformed phyloseq object
    new_otu_table <- otu_table(x_rclr, taxa_are_rows = FALSE)
    physeq_rclr <- physeq_rmp_glom
    otu_table(physeq_rclr) <- new_otu_table

    if (tax == "ASV") {
      tax_folder <- clean_rds_folder
    } else {
      tax_folder <- file.path(clean_rds_folder, tax)
      if (!dir.exists(tax_folder)) dir.create(tax_folder, recursive = TRUE)
    }

    saveRDS(physeq_rclr, file = file.path(tax_folder, glue::glue("{project_id}_phyloseq_rcrl_{tax}.rds")))

    # Build coordinates data frame
    df_rclr <- dplyr::tibble(
      Sample = rownames(x_rclr),
      PC1 = pca_rclr$x[, 1],
      PC2 = pca_rclr$x[, 2],
      PC3 = pca_rclr$x[, 3]
    )

    # Extract an clean sample metadata
    meta_df <- phyloseq::sample_data(physeq_rclr) %>%
      data.frame() %>%
      tibble::rownames_to_column(var = "SampleID") %>%
      dplyr::as_tibble()

    # Combine PCA coordinates with metadata
    df_rclr_meta <- df_rclr %>%
      inner_join(., meta_df, by = c("Sample" = "SampleID"))

    # Calculate variance percentage for axes labels
    # Calculate percentage of variance explanied for each PC axis
    pca_variance <- (pca_rclr$sdev^2) / sum(pca_rclr$sdev^2) * 100
    pc1_label <- paste0("PC1 (", round(pca_variance[1], 1), "%)")
    pc2_label <- paste0("PC2 (", round(pca_variance[2], 1), "%)")
    pc3_label <- paste0("PC3 (", round(pca_variance[3], 1), "%)")

    mapping_list <- list()
    if (!is.null(color_var)) { mapping_list$color <- color_var }
    if (!is.null(shape_var)) { mapping_list$shape <- shape_var }

    rclr_beta_div <- ggplot(df_rclr_meta, aes(x = PC1, y = PC2)) +
      geom_point(mapping = aes(!!!lapply(mapping_list, ggplot2::sym)), alpha = 0.8, size = 2.5) +
      scale_color_gradient(low = "#2C5282", high = "#C05621", name = "Time (days)") +
      labs(
        title = "Beta Diversity PCA on RCLR",
        x = pc1_label,
        y = pc2_label
      ) +
      theme_light() +
      theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13, margin = margin(b = 4)),
        axis.title = element_text(color = "gray20", face = "plain"),
        axis.text = element_text(color = "gray40"),
        legend.position = "right",
        strip.background = element_rect(fill = "gray95", color = "gray85"),
        strip.text = element_text(face = "bold", color = "gray20")
      )

    if (!is.null(facet_var)) {
      rclr_beta_div <- rclr_beta_div + facet_wrap(vars(!!sym(facet_var)))
    }

    ggsave(filename = file.path(beta_div_folder, glue::glue("rclr_beta_diversity_{tax}.pdf")), plot = rclr_beta_div, width = 12, height = 8)
    ggsave(filename = file.path(beta_div_folder, glue::glue("rclr_beta_diversity_{tax}.png")), plot = rclr_beta_div, width = 12, height = 8, dpi = 600)

    # RCLR Loading --------------------------------------------------------
    # Function to generate loadings plot for a specific PC
    create_loadings_plot <- function(pc_axis, data, current_tax) {

      label_col <- if (current_tax == "ASV") "Taxon" else current_tax

      # Select top taxa based on absolute impact on the specified PC axis
      top_taxa <- data %>%
        arrange(desc(abs(.data[[pc_axis]]))) %>%
        slice_head(n = 20) %>%
        mutate(
          direction = if_else(.data[[pc_axis]] > 0, "Positive Loading", "Negative Loading"), # Identify if the effect is positive or negative
          taxon_unique = make.unique(as.character(!!sym(tax))), # Creates unique names
          taxon_unique = reorder(taxon_unique, .data[[pc_axis]]) # Orders the axis by absolute impact
        )

      # Generate the plot
      ggplot(top_taxa, aes(x = taxon_unique, y = .data[[pc_axis]], fill = direction)) +
        geom_col(alpha = 0.8, width = 0.75) +
        geom_hline(yintercept = 0, linetype = "solid", color = "gray40", linewidth = 0.4) +
        scale_fill_manual(values = c("Positive Loading" = "#4A90E2", "Negative Loading" = "#D0021B")) +
        coord_flip() +
        theme_classic() +
        scale_x_discrete(labels = function(x) gsub("\\.\\d+$", "", x)) +
        scale_y_continuous(expand = expansion(mult = c(0.08, 0.08))) +
        labs(title = pc_axis, x = tax, y = pc_axis) +
        theme(
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_line(color = "gray92", linewidth = 0.5),
          plot.title = element_text(face = "bold", size = 13, margin = margin(b = 4)),
          axis.title.x = element_text(margin = margin(t = 10), face = "plain", color = "gray20"),
          axis.title.y = element_text(margin = margin(r = 10), face = "plain", color = "gray20"),
          axis.text = element_text(color = "gray30"),
          legend.position = "none"
        )
    }

    # Clean taxonomy table for loadings analysis
    tax_table_clean <- physeq_rmp_glom %>%
      psmelt() %>%
      as_tibble() %>%
      select(OTU, any_of(tax)) %>%
      distinct(OTU, .keep_all = TRUE)

    # Join PCA rotations with taxonomic names
    loadings_rclr <- as_tibble(pca_rclr$rotation, rownames = "Taxon") %>%
      dplyr::left_join(tax_table_clean, by = c("Taxon" = "OTU"))

    # Generate separate plots for PC1, PC2 and PC3
    if (tax != "ASV") {
      plot_pc1 <- create_loadings_plot(pc_axis = "PC1", data = loadings_rclr, tax)
      plot_pc2 <- create_loadings_plot(pc_axis = "PC2", data = loadings_rclr, tax)
      plot_pc3 <- create_loadings_plot(pc_axis = "PC3", data = loadings_rclr, tax)

      # Assemble the three loadings plots side-by-side into a single row
      combined_loadings <- cowplot::plot_grid(
        plot_pc1,
        plot_pc2,
        plot_pc3,
        ncol = 3,
        labels = c("A", "B", "C")
      )

      # Create a clean main title
      main_title <- cowplot::ggdraw() +
        cowplot::draw_label(
          label = paste0("Top ", tax, " Loadings Driving RCLR"),
          fontface = "bold",
          size = 15,
          hjust = 0.5
        )

      # Combine the main title and the loadings
      final_plot <- cowplot::plot_grid(
        main_title,
        combined_loadings,
        ncol = 1,
        rel_heights = c(0.1, 0.9)
      )

      ggsave(filename = file.path(rclr_load_folder, glue::glue("rclr_loadings_{tax}.pdf")), plot = final_plot, width = 16, height = 8)
      ggsave(filename = file.path(rclr_load_folder, glue::glue("rclr_loadings_{tax}.png")), plot = final_plot, width = 16, height = 8, dpi = 600)
    }

    # Optional plots if time var is present
    if (!is.null(time_var)) {

      if (tax != "ASV") {
        # Select top n taxa wiht highest absolute loadings impact
        top_load_taxa <- loadings_rclr %>%
          dplyr::mutate(max_abs_load = pmax(abs(PC1), abs(PC2), abs(PC3), na.rm = TRUE)) %>%
          dplyr::arrange(desc(max_abs_load)) %>%
          dplyr::slice_head(n = 20) %>%
          dplyr::pull(Taxon)

        psdata_rclr <- physeq_rclr %>%
          phyloseq::psmelt() %>%
          tibble::as_tibble() %>%
          dplyr::filter(OTU %in% top_load_taxa)

        rclr_time_folder <- file.path(rclr_folder, "RCLR_Time")
        if (!dir.exists(rclr_time_folder)) {dir.create(rclr_time_folder, recursive = TRUE)}

        # RCLR Abundance vs Time
        for (load_tax in top_load_taxa) {
          psdata_rclr_tax <- psdata_rclr %>%
            dplyr::filter(OTU == load_tax)

          tax_name <- psdata_rclr_tax %>%
            dplyr::select(OTU, !!sym(tax)) %>%
            BiocGenerics::unique() %>%
            dplyr::pull(!!sym(tax))

          tax_name_clean <- gsub("[[:space:]/]+", "_", tax_name)
          tax_name_clean <- gsub("[^[:alnum:]_]", "", tax_name_clean)
          tax_name_clean <- gsub("_{2,}", "_", tax_name_clean)

          plot_rclr_time <- ggplot(psdata_rclr_tax, aes(x = !!sym(time_var), y = Abundance)) +
            geom_point(alpha = 0.8, size = 2) +
            geom_smooth(method = "loess", color = "black", se = FALSE, linewidth = 0.8, formula = y ~ x) +
            labs(
              title = glue::glue("RCLR over Time ({tax}: {tax_name_clean})"),
              x = time_var,
              y = "RCLR"
            ) +
            theme_light() +
            theme(
              panel.grid.minor = element_blank(),
              strip.background = element_rect(fill = "gray95", color = "gray85"),
              strip.text = element_text(face = "bold", color = "gray20"),
              plot.title = element_text(face = "bold", size = 13),
              axis.text = element_text(color = "gray30"),
              axis.title = element_text(color = "gray20")
            )

          if (!is.null(facet_var)) {
            plot_rclr_time <- plot_rclr_time + facet_wrap(vars(!!sym(facet_var)))
          }

          ggsave(filename = file.path(rclr_time_folder, glue::glue("rclr_over_time_{tax}_{tax_name_clean}.pdf")), plot = plot_rclr_time, width = 16, height = 8)
          ggsave(filename = file.path(rclr_time_folder, glue::glue("rclr_over_time_{tax}_{tax_name_clean}.png")), plot = plot_rclr_time, width = 16, height = 8, dpi = 600)
        }
      }

      # PC scores vs Time
      plot_pc_time <- function(df, pc) {
        p <- ggplot(df, aes(x = !!sym(time_var), y = .data[[pc]])) +
          geom_point(alpha = 0.8, size = 2) +
          geom_smooth(method = "loess", color = "steelblue", se = FALSE, linewidth = 0.8, formula = y ~ x) +
          labs(
            title = glue::glue("{pc} Scores over Time ({tax})"),
            x = time_var,
            y = glue::glue("{pc}")
          ) +
          theme_light() +
          theme(
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = "gray95", color = "gray85"),
            strip.text = element_text(face = "bold", color = "gray20"),
            plot.title = element_text(face = "bold", size = 13),
            axis.text = element_text(color = "gray30"),
            axis.title = element_text(color = "gray20")
          )

        if (!is.null(facet_var)) {
          p <- p + facet_wrap(vars(!!sym(facet_var)))
        }
        return(p)
      }

      pc_time_folder <- file.path(rclr_folder, "PC_Time")
      if (!dir.exists(pc_time_folder)) {dir.create(pc_time_folder, recursive = TRUE)}

      rclr_p1 <- plot_pc_time(df_rclr_meta, "PC1")
      rclr_p2 <- plot_pc_time(df_rclr_meta, "PC2")
      rclr_p3 <- plot_pc_time(df_rclr_meta, "PC3")

      combined_pcs_time <- cowplot::plot_grid(
        rclr_p1,
        rclr_p2,
        rclr_p3,
        ncol = 1,
        align = "v",
        labels = c("A", "B", "C")
      )

      ggsave(filename = file.path(pc_time_folder, glue::glue("pc_over_time_{tax}.pdf")), plot = combined_pcs_time, width = 16, height = 12)
      ggsave(filename = file.path(pc_time_folder, glue::glue("pc_over_time_{tax}.png")), plot = combined_pcs_time, width = 16, height = 12, dpi = 600)

    }
  }
}
