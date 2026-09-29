#' Generate Beta Diversity Ordination Plots
#'
#' This function computes beta diversity using ordination methods (e.g., PCoA) on a phyloseq object and generates corresponding ordination plots. It supports both ASV-level data and data aggregated at various taxonomic levels, and it can handle both relative and, if available, absolute abundance data. Four distance metrics are used for relative abundance plots (Jaccard, Bray-Curtis, Unweighted UniFrac, and Weighted UniFrac), and Manhattan distance is used for absolute abundance plots when a normalization method is provided.
#'
#' @param physeq A phyloseq object containing microbial community data.
#' @param taxrank A character vector specifying the taxonomic levels to process. If the first element (case-insensitive) is \code{"asv"}, ASV-level beta diversity is computed; otherwise, beta diversity is computed for each taxonomic level provided (default: \code{c("Phylum", "Class", "Order", "Family", "Genus")}).
#' @param norm_method A character string indicating the normalization method to be used for generating absolute abundance data. Options include:
#'   \itemize{
#'     \item \code{"fcm"}: Flow cytometry normalization.
#'     \item \code{"qpcr"}: qPCR normalization.
#'     \item \code{NULL}: No absolute abundance processing (default).
#'   }
#' @param ordination_method A character string specifying the ordination method to use (default is \code{"PCoA"}).
#' @param color_factor An optional character string specifying the sample metadata column used to color the points in the ordination plot.
#' @param color_continuous A logical value indicating whether the color scale should be continuous (\code{TRUE}) or discrete (\code{FALSE}). Default is \code{TRUE}.
#' @param shape_factor An optional character string specifying the sample metadata column used to assign point shapes in the ordination plot.
#' @param size_factor An optional character string specifying the sample metadata column used to assign point sizes in the ordination plot.
#' @param alpha_factor An optional character string specifying the sample metadata column used to assign point transparency in the ordination plot.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Sets up project directories and folder paths for saving output files.
#'   \item Defines an internal helper function, \code{base_beta_plot}, which performs the ordination (using the specified \code{ordination_method} and distance metric), removes default point layers, and then adds a customized geom_point layer with aesthetics defined by the provided factors.
#'   \item For ASV-level data (when \code{taxrank[1]} is \code{"asv"}):
#'     \itemize{
#'       \item Processes relative abundance data (transformed to percentages) and, if a normalization method is specified, also processes absolute abundance data.
#'       \item Generates ordination plots using multiple distance metrics:
#'         \itemize{
#'           \item \strong{Jaccard}: Based on binary presence/absence.
#'           \item \strong{Bray-Curtis}: Incorporates both presence and abundance.
#'           \item \strong{Unweighted UniFrac}: Considers lineage presence only.
#'           \item \strong{Weighted UniFrac}: Considers both lineage presence and abundance.
#'         }
#'       \item If both DNA and RNA data are present, separate plots are generated for each.
#'     }
#'   \item For other taxonomic levels:
#'     \itemize{
#'       \item Similar processing is performed for each taxonomic rank, with output saved in dedicated subfolders.
#'     }
#'   \item Aesthetic scales (colors, shapes, sizes, alpha) are defined based on the unique levels in the sample data.
#' }
#'
#' @return A ggplot object representing the combined beta diversity ordination plot for the relative (and, if applicable, absolute) data.
#'
#' @examples
#' \dontrun{
#'   # Example: Generate beta diversity plots at the ASV level using PCoA ordination and flow cytometry normalization
#'   beta_plot <- beta_diversity(
#'     physeq = my_physeq,
#'     taxrank = "asv",
#'     norm_method = "fcm",
#'     ordination_method = "PCoA",
#'     color_factor = "Treatment",
#'     shape_factor = "Replica",
#'     size_factor = "Timepoint"
#'   )
#'
#'   # Example: Generate beta diversity plots for Phylum and Class levels without absolute data processing
#'   beta_plot <- beta_diversity(
#'     physeq = my_physeq,
#'     taxrank = c("Phylum", "Class"),
#'     norm_method = NULL,
#'     ordination_method = "PCoA",
#'     color_factor = "Soil_Type"
#'   )
#' }
#'
#' @export
beta_diversity <- function(physeq, taxrank, norm_method, ordination_method,
                           color_var, shape_var, size_var, alpha_var, facet_var,
                           project_id, base_path, log_file) {

  # Function for creating the base beta-diversity plot
  base_beta_plot <- function(physeq, ordination_method, distance_method, title,
                             color_var, shape_var, size_var, alpha_var, facet_var) {

    # Perform the ordination
    ordination_res <- ordinate(physeq, method = ordination_method, distance = distance_method)

    # Build a dynamic mapping list to force phyloseq to load ALL required metadata columns
    mapping_list <- list()
    if (!is.null(color_var)) mapping_list$color <- color_var
    if (!is.null(shape_var)) mapping_list$shape <- shape_var
    if (!is.null(size_var))  mapping_list$size  <- size_var
    if (!is.null(alpha_var)) mapping_list$alpha <- alpha_var
    if (!is.null(facet_var)) mapping_list$facet <- facet_var

    base_plot <- plot_ordination(physeq, ordination = ordination_res, axes = c(1, 2))
    base_plot$mapping <- utils::modifyList(base_plot$mapping, ggplot2::aes(!!!lapply(mapping_list, ggplot2::sym)))

    # Safely clear default phyloseq point layers to prevent doubling
    base_plot$layers <- list()

    # Add points with the specified aesthetics
    base_plot <- base_plot +
      geom_point(stroke = 1, size = if (is.null(size_var)) 3 else NULL) +
      ggtitle(title) +
      labs(color = color_var, shape = shape_var, size = size_var, alpha = alpha_var) +
      theme_classic() +
      theme(
        panel.background = element_rect(fill = "transparent"),
        panel.grid = element_line(colour = "grey90"),
        strip.text = element_text(face = "bold"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
        axis.line.y = element_line(color = "black", linewidth = 0.3),
        axis.line.x = element_line(color = "black", linewidth = 0.3)) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)))  # 5% expansion for y-axis

    if (!is.null(facet_var)) {
      base_plot <- base_plot + facet_wrap(vars(!!sym(facet_var)))
    }

    return(base_plot)
  }

  # Set up project and folder paths
  figures_folder <- file.path(base_path, "03_figures")

  beta_div_folder <- file.path(figures_folder, "beta_diversity")
  if(!dir.exists(beta_div_folder)) { dir.create(beta_div_folder, recursive = TRUE) }

  for (tax in taxrank) {
    physeq_rmp <- physeq$physeq_rmp_rarefied

    if (tax != "ASV") {
      physeq_rmp_glom <- phyloseq::tax_glom(physeq_rmp, taxrank = tax)
    } else {
      physeq_rmp_glom <- physeq_rmp
    }

    # Check if a phylogenetic tree is present in the object
    has_tree <- !is.null(phyloseq::phy_tree(physeq_rmp_glom, errorIfNULL = FALSE))

    # Create beta-diversity plots
    plot_Jac <- base_beta_plot(physeq_rmp_glom, ordination_method, "jaccard", "Jaccard\n(binary presence)",
                               color_var, shape_var, size_var, alpha_var, facet_var)
    plot_BC <- base_beta_plot(physeq_rmp_glom, ordination_method, "bray", "Bray-Curtis\n(presence + abundance)",
                              color_var, shape_var, size_var, alpha_var, facet_var)

    # Extract shared legend safely before stripping panel themes
    shared_legend <- cowplot::get_legend(plot_Jac + ggplot2::theme(legend.position = "right"))

    # Strip legends from standard panels
    plot_Jac <- plot_Jac + ggplot2::theme(legend.position = "none")
    plot_BC  <- plot_BC  + ggplot2::theme(legend.position = "none")

    # Conditionally generate UniFrac plots only if a tree exists
    if (has_tree) {
      log_message(paste("Phylogenetic tree detected. Generating UniFrac plots for", tax), log_file)

      plot_uu <- base_beta_plot(physeq_rmp_glom, ordination_method, "uunifrac", "Unweighted UniFrac\n(lineage presence)",
                                color_var, shape_var, size_var, alpha_var, facet_var) + ggplot2::theme(legend.position = "none")
      plot_wu <- base_beta_plot(physeq_rmp_glom, ordination_method, "wunifrac", "Weighted UniFrac\n(lineage abundance)",
                                color_var, shape_var, size_var, alpha_var, facet_var) + ggplot2::theme(legend.position = "none")

      # Assemble 4-panel grid layout
      combined_panels <- cowplot::plot_grid(plot_Jac, plot_BC, plot_uu, plot_wu, ncol = 2, labels = c("A", "B", "C", "D"))
      plot_ncol <- 2
    } else {
      log_message(paste("No phylogenetic tree detected. Skipping UniFrac plots for", tax), log_file)

      # Assemble 2-panel grid layout (only Jaccard and Bray-Curtis side-by-side)
      combined_panels <- cowplot::plot_grid(plot_Jac, plot_BC, ncol = 2, labels = c("A", "B"))
      plot_ncol <- 1
    }

    # Assemble final grid visualization with the shared legend
    final_grid_plot <- cowplot::plot_grid(combined_panels, shared_legend, ncol = 2, rel_widths = c(3, 0.6))

    ggsave(filename = file.path(beta_div_folder, glue::glue("beta_diversity_rmp_{tax}.png")), plot = final_grid_plot, width = 12, height = 10, dpi = 600)
    ggsave(filename = file.path(beta_div_folder, glue::glue("beta_diversity_rmp_{tax}.pdf")), plot = final_grid_plot, width = 12, height = 10)

  }
}
