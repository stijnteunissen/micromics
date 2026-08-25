#' Generate and Save Rarefaction Curve
#'
#' Generates a sequencing depth rarefaction curve from a phyloseq object
#' and exports the resulting plot as a PDF.
#'
#' @param physeq A phyloseq object to be analyzed.
#' @param color Character string specifying the metadata column used for coloring lines.
#'   Defaults to \code{"sample_or_control"}.
#' @param base_path Character string specifying the root project directory.
#' @param project_id Character string specifying the unique project name.
#' @param log_file Character string specifying the path to the log file.
#'
#' @return None. This function is called for its side effects (saving a PDF file).
#' @export
rarefaction_curve <- function(physeq, color = "sample_or_control", project_id, base_path, log_file) {

  log_message("Generating rarefaction curve plot on raw ASV levels.", status = "start", log_file)

  # Define storage directory path cleanly
  figure_folder <- file.path(base_path, "figures")

  # Validate mandatory core coloring column presence
  if (!"sample_or_control" %in% colnames(phyloseq::sample_data(physeq))) {
    error_message <- "Error: 'sample_or_control' column is missing from the sample metadata."
    log_message(error_message, status = "error", log_file)
    stop(error_message, call. = FALSE)
  }

  log_message(glue::glue("Using '{color}' column for line coloring in rarefaction plot."), status = "info", log_file)

  # Internal helper function to plot the underlying vegan curve
  plot_rarecurve <- function(data, step = 100, ylim_val = NULL, xlim_val = NULL, color_col = NULL) {
    abund_table <- as.data.frame(phyloseq::otu_table(data))

    if (!is.null(color_col)) {
      # Extract metadata vector safely using S4 slot extraction conventions
      metadata_df <- as.data.frame(phyloseq::sample_data(data))
      group_vector <- as.character(metadata_df[[color_col]])
      unique_groups <- unique(group_vector)

      # Generate standardized color hues
      hues <- seq(15, 375, length = length(unique_groups) + 1)
      color_palette <- hcl(h = hues, l = 65, c = 100)[seq_along(unique_groups)]

      # Map samples to their corresponding group color
      sample_colors <- color_palette[match(group_vector, unique_groups)]
    } else {
      sample_colors <- "black"
    }

    # Setup flexible plot limits based on per-sample ASV maximum
    if (is.null(xlim_val)) {
      # Open a silent graphics device to calculate max depth per sample
      pdf(NULL)
      captured_curve <- vegan::rarecurve(t(abund_table), step = 100, label = FALSE)
      dev.off()
      xlim_val <- c(0, max(sapply(captured_curve, function(x) max(attr(x, "Subsample")))))
    }

    if (is.null(ylim_val)) {
      # Dynamically scale y-axis to the maximum number of ASVs found within a single sample
      pdf(NULL)
      captured_curve <- vegan::rarecurve(t(abund_table), step = 100, label = FALSE)
      dev.off()
      ylim_val <- c(0, max(sapply(captured_curve, function(x) max(x))))
    }

    # Draw the standardized base graphic curve
    vegan::rarecurve(
      x = t(abund_table),
      step = step,
      xlim = xlim_val,
      ylim = ylim_val,
      label = FALSE,
      col = sample_colors,
      xlab = "Sequencing depth",
      ylab = "Total ASV"
    )

    # Add legend if coloring variables are present
    if (!is.null(color_col)) {
      legend("bottomright", legend = unique_groups, fill = color_palette, bty = "n")
    }
  }

  # Construct output destination file path
  output_pdf_path <- file.path(figure_folder, glue::glue("{project_id}_rarefaction_curve.pdf"))

  # Open graphics engine and write plot stream
  grDevices::cairo_pdf(file = output_pdf_path, width = 7, height = 5)
  plot_rarecurve(physeq, color_col = color, xlim_val = c(0, max(phyloseq::sample_sums(physeq))))
  grDevices::dev.off()

  # Log successful termination boundaries
  log_message(glue::glue("Rarefaction curve PDF stored successfully at: {output_pdf_path}"), status = "info", log_file)
  log_message("Rarefaction curve visualization complete.", status = "success", log_file)
}
