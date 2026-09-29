#' Generate Alpha Diversity Plots
#'
#' This function calculates alpha diversity metrics from a phyloseq object and
#' generates alpha diversity plots at various taxonomic levels or at the ASV
#' level. The alpha diversity measures (Observed, Chao1, Shannon, Simpson) are
#' computed using the \code{estimate_richness} function from the phyloseq
#' package. Depending on the input parameters, the function can handle
#' normalized data (using flow cytometry or qPCR methods) and can separate plots
#' by DNA and RNA types.
#'
#' @param physeq A phyloseq object containing microbial community data.
#' @param norm_method A character string specifying the normalization method for the data. Options include:
#'   \itemize{
#'     \item \code{"fcm"}: Use flow cytometry-normalized data.
#'     \item \code{"qpcr"}: Use qPCR-normalized data.
#'     \item \code{NULL}: Use only copy number corrected data (default).
#'   }
#' @param taxrank A character vector specifying the taxonomic levels for which
#'   alpha diversity is to be calculated. If the first element is (taxrank = `asv`),
#'   ASV-level data is processed. Otherwise, the function
#'   processes data for each taxonomic level provided default is taxrank = c('Phylum', 'Class', 'Order', 'Family', 'Genus').
#' @param date_factor An optional character string indicating the name of the
#'   date column in the sample metadata. If provided, the column is converted to
#'   a Date object ("%d/%m/%Y") and used to order the data.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Extracts the appropriate data object from \code{physeq} based on the chosen \code{norm_method} and taxonomic level.
#'   \item Estimates alpha diversity metrics (Observed, Chao1, Shannon, Simpson) using \code{estimate_richness}.
#'   \item Merges the alpha diversity estimates with sample metadata.
#'   \item Optionally orders the data by a date factor if provided.
#'   \item Appends dummy data to the dataset for visualization purposes.
#'   \item Exports the combined alpha diversity data as CSV files.
#'   \item Generates bar plots for the diversity metrics, optionally separating DNA and RNA data if both are present.
#'   \item Saves the resulting plots as PDF files in the project's figures folder.
#' }
#'
#' @return A combined ggplot object containing the generated alpha diversity plots.
#'
#' @examples
#' \dontrun{
#'   # Generate alpha diversity plots at the ASV level without normalization
#'   alpha_plot <- alpha_diversity(physeq = my_physeq, taxrank = "asv")
#'
#'   # Generate alpha diversity plots at the Phylum level using flow cytometry-normalized data
#'   alpha_plot <- alpha_diversity(physeq = my_physeq, norm_method = "fcm", taxrank = "Phylum")
#'
#'   # Generate alpha diversity plots with a specified date factor for ordering samples
#'   alpha_plot <- alpha_diversity(physeq = my_physeq, date_factor = "Sample_Date")
#' }
#'
#' @export
alpha_diversity = function(physeq, norm_method, taxrank, facet_vars, plot_width, project_id, base_path, log_file) {

  base_alpha_plot = function(alpha_data, x_value, y_value, x_label, y_label, facet_vars) {

    p = ggplot(alpha_data, aes(x = !!sym(x_value), y = !!sym(y_value))) +
      geom_point(color = "black", size = 3, alpha = 0.8) +
      #geom_col(fill = "steelblue", color = "steelblue", show.legend = FALSE) +
      theme_classic() +
      labs(x = x_label, y = y_label) +
      theme(
        legend.position = "none",
        legend.text = element_markdown(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(face = "bold", angle = 90, vjust = 0.5, hjust = 0),
        strip.placement = "outside",
        strip.text = element_text(face = "bold"),
        strip.background = element_blank(),
        ggh4x.facet.nestline = element_line(colour = "black")
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +  # Correct placement outside theme with 5%
      expand_limits(y = c(min(alpha_data[[y_value]]) - 1, max(alpha_data[[y_value]]) + 1))  # Correct placement outside theme

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

  alpha_div_folder <- file.path(figures_folder, "alpha_diversity")
  if(!dir.exists(alpha_div_folder)) { dir.create(alpha_div_folder, recursive = TRUE) }

  alpha_div_table_folder <- file.path(export_folder, "alpha_diversity")
  if(!dir.exists(alpha_div_table_folder)) { dir.create(alpha_div_table_folder, recursive = TRUE) }

  # ASV and taxa
  for (tax in taxrank) {
    physeq_rmp <- physeq$physeq_rmp_rarefied

    if (tax != "ASV") {
      physeq_rmp_glom <- phyloseq::tax_glom(physeq_rmp, taxrank = tax)
    } else {
      physeq_rmp_glom <- physeq_rmp
    }

    alpha_div = phyloseq::estimate_richness(physeq_rmp_glom, measures = c("Observed", "Chao1", "Shannon", "Simpson"))
    alpha_div_df = alpha_div %>% tibble::rownames_to_column(var = "SampleID")
    metadata = phyloseq::sample_data(physeq_rmp_glom) %>%
      data.frame() %>%
      tibble::rownames_to_column(var = "SampleID") %>%
      dplyr::as_tibble()
    alpha_div_df_meta = inner_join(metadata, alpha_div_df, by = "SampleID")

    alpha_div_df_meta_export <- alpha_div_df_meta %>%
      dplyr::mutate(
        Observed = round(Observed, 2),
        Chao1 = round(Chao1, 2),
        Shannon = round(Shannon, 2),
        Simpson = round(Simpson, 2)) %>%
      select(SampleID, read_count, Observed, Chao1, Shannon, Simpson)

    readr::write_csv(alpha_div_df_meta_export, file.path(alpha_div_table_folder, glue::glue("alpha_diversity_data_{tax}.csv")))

    Chao1_plot <- base_alpha_plot(alpha_div_df_meta, "SampleID", "Chao1", x_label = "Sample", y_label = "Choa1 Index", facet_vars)
    Shannon_plot <- base_alpha_plot(alpha_div_df_meta, "SampleID", "Shannon", x_label = "Sample", y_label = "Shannon Index", facet_vars)

    combined_plot <- cowplot::plot_grid(Chao1_plot, Shannon_plot, align = "v", labels = c("A", "B"), ncol = 1)

    ggsave(filename = file.path(alpha_div_folder, glue::glue("alpha_diversity_rmp_{tax}.png")), plot = combined_plot, width = plot_width, height = 10, dpi = 600)
    ggsave(filename = file.path(alpha_div_folder, glue::glue("alpha_diversity_rmp_{tax}.pdf")), plot = combined_plot, width = plot_width, height = 10)
  }
}
