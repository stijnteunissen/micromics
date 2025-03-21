% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/creating_physeq_object.R
\name{creating_physeq_object}
\alias{creating_physeq_object}
\title{Create a Phyloseq Object}
\usage{
creating_physeq_object(projects)
}
\arguments{
\item{projects}{A character vector containing the name of the project folders.
Each project folder must already exist within the specified
base path. The \code{projects} argument represents the name of the
project folder, which must contain a subfolder named \code{qiime2_output}
with the output files generated from QIIME2. This folder structure
is necessary for setting up downstream analysis.
\code{projects} is typically reference in a configuration file.}
}
\value{
A \code{phyloseq} object that integrates feature tables, taxonomy, phylogenetic trees, and metadata.
This object is also saved as an RDS file for further usage in downstream analyses.
}
\description{
This function creates a \code{phyloseq} object using input data files such as
feature tables, taxonomic assignments, phylogenetic tree, and unified
metadata into a single \code{phyloseq} object for downstream analysis.
}
\details{
This function performs the following steps:
\itemize{
\item Defines the paths to the required input files (feature table, rooted tree, taxonomy, and metadata).
\item Searches for and retrieves these files from the \code{input_data} directory.
\item Calls the \code{qza_to_phyloseq()} function to create a \code{phyloseq} object using the input files.
\item Adds read count information to the sample metadata within the \code{phyloseq} object.
}

The function assumes the following files are present in the \code{input_data} directory:
\itemize{
\item \code{table.qza}: Feature table containing sample feature data.
\item \code{rooted-tree.qza}: Phylogenetic tree.
\item \code{classifier.qza}: Taxonomic classification file.
\item \code{metadata_formatted.tsv}: Unified sample metadata.
}

The resulting \code{phyloseq} object is essential for downstream analyses and integrates all input files into a single, structured object.
The created \code{phyloseq} object is saved as an RDS file named \verb{<project_name>_uncleaned_phyloseq.rds} in the \code{output_data/rds_files} directory.
}
\examples{
\dontrun{
# Example usage:
physeq_object <- creating_physeq_object(projects)
}

}
