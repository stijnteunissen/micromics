% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/creating_folders_test_script.R
\name{creating_folders}
\alias{creating_folders}
\title{Create Project Folders and Process Files}
\usage{
creating_folders(projecten)
}
\arguments{
\item{projecten}{A character vector containing the names of projects. Each project will have its own folder structure created.}
}
\value{
No return value. The function creates directories and copies files as a side effect.
}
\description{
This function creates the necessary directory structure for specified projects,
checks for required and optional files in the source directory, and copies
relevant files to the appropriate locations.
}
\details{
The function performs the following steps for each project:
\itemize{
\item Creates directories if they do not already exist.
\item Checks for required files (\code{table.qza}, \code{rooted-tree.qza}, \code{classifier.qza}, \code{metadata.tsv}).
\item Logs errors if required files are missing and stops execution.
\item Checks for optional files (\code{dna-sequences.csv}, \code{qPCR.csv}, \code{fcm.csv}, \code{metadata_extra.tsv}).
\item Logs warnings if optional files are missing.
\item Copies relevant files to the \code{input_R_data} directory.
}
}
\examples{
# Example usage:
projecten <- c("Project1", "Project2")
creating_folders(projecten)
}
