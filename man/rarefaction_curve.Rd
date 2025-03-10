% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rarefaction_curve.R
\name{rarefaction_curve}
\alias{rarefaction_curve}
\title{Generate and Save Rarefaction Curve}
\usage{
rarefaction_curve(
  physeq = resolved_tree_physeq,
  rank = "ASV",
  processing_stage,
  color = NULL,
  base_path = "./"
)
}
\arguments{
\item{physeq}{A phyloseq object containing OTU data and associated metadata.}

\item{rank}{The rank level at which taxa should be aggregated. Default is \code{rank}.}

\item{processing_stage}{A character string describing the processing stage (e.g., "raw" or "filtered").
This string will be included in the saved PDF file name.}

\item{color}{A character string specifying the column in the sample metadata to use for
coloring the samples. Default is \code{NULL}, which automatically sets the color
to \code{"sample_type"}.}

\item{base_path}{A character string specifying the base directory for saving output files.
Default is \code{"./"}.}
}
\value{
The rarefaction curve plot object.
}
\description{
This function creates a rarefaction curve for a given phyloseq object and
saves the plot as a PDF.
}
\details{
This function first checks whether the \code{sample_type} column exists in the
sample metadata. It then generates a rarefaction curve using the \code{amp_rarecurve}
function and saves the plot as a PDF file in the specified directory.
}
\examples{
rarefaction_curve(physeq = physeq, rank = "Genus", processing_stage = "filtered", color = "Treatment")
}
