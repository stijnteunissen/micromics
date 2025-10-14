amp_rarecurve <- function(data, step = 100, ylim = NULL, xlim = NULL, label = FALSE, color = NULL, legend = TRUE, color.vector = NULL, legend.position = "topleft") {
  abund = otu_table(data)
  abund = as.data.frame(abund)

  max_depth <- max(rowSums(abund))

  if (!is.null(color)) {
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    group_vector <- sample_data(data)[, color]@.Data %>% as.data.frame()
    names(group_vector) <- "color_variable"
    group_vector <- as.character(group_vector$color_variable)
    groups <- unique(group_vector)
    n = length(groups)
    cols = gg_color_hue(n)
    if (!is.null(color.vector)) {
      cols <- color.vector
    }
    col_vector <- rep("black", length(group_vector))
    for (i in seq_along(group_vector)) {
      col_vector[i] <- cols[match(group_vector[i], groups)]
    }
  } else {
    col_vector = "black"
  }
  if (is.null(ylim) & is.null(xlim)) {
    vegan::rarecurve(t(abund), step = step, label = label, col = col_vector, xlab = "Sequencing depth", ylab = "Total ASV")
  }
  if (!is.null(ylim) & !is.null(xlim)) {
    vegan::rarecurve(t(abund), step = step, ylim = ylim, xlim = c(0, max_depth), label = label, col = col_vector, xlab = "Sequencing depth", ylab = "Total ASV")
  }
  if (!is.null(ylim) & is.null(xlim)) {
    vegan::rarecurve(t(abund), step = step, ylim = ylim, label = label, col = col_vector, xlab = "Sequencing depth", ylab = "Total ASV")
  }
  if (is.null(ylim) & !is.null(xlim)) {
    vegan::rarecurve(t(abund), step = step, xlim = c(0, max_depth), label = label, col = col_vector, xlab = "Sequencing depth", ylab = "Total ASV")
  }
  if (!is.null(color) & legend) {
    legend(legend.position, legend = groups, fill = cols, bty = "n")
  }
}
#' @title Generate and Save Rarefaction Curve
#' @description This function creates a rarefaction curve for a given phyloseq object and
#'              saves the plot as a PDF.
#'
#' @inheritParams tax_clean
#'
#' @param color A character string specifying the column in the sample metadata to use for
#'              coloring the samples. Default is `NULL`, which automatically sets the color
#'              to `"sample_or_control"`.
#'
#' @return The rarefaction curve plot object.
#' @details This function first checks whether the `sample_or_control` column exists in the
#'          sample metadata. It then generates a rarefaction curve using the `amp_rarecurve`
#'          function and saves the plot as a PDF file in the specified directory.
#' @examples
#' rarefaction_curve(physeq = physeq)
#'
#' @export
rarefaction_curve <- function(physeq = resolved_tree_physeq, color = NULL) {

  log_message(paste("Step 6: Creating rarefaction curve: creating rarefaction curve before cleaning on ASV level.", paste(projects, collapse = ", ")), log_file)

  if (!("sample_or_control" %in% colnames(sample_data(physeq)))) {
    stop("error: 'sample_or_control' column is missing from the sample data.")
  }

  if (is.null(color)) {
    color = "sample_or_control"
  } else {
    message(paste("message: using", color, "for color in rarefaction plot."))
  }

  project_name = projects
  project_folder = paste0(base_path, project_name)
  figure_folder_pdf = paste0(project_folder, "/figures/PDF_figures/")
  if(!dir.exists(figure_folder_pdf)) { dir.create(figure_folder_pdf, recursive = TRUE) }

  figure_path <- file.path(figure_folder_pdf, paste0(project_name, "_rarefaction_curve.pdf"))
  cairo_pdf(file = figure_path, width = 7, height = 5)
  amp_rarecurve(physeq, color = color, legend.position = "bottomright", xlim = c(0, max(sample_sums(physeq))))
  dev.off()
  log_message("Created successfully rarefaction curve.", log_file)
}
