#' Resolve Multichotomies in Phyloseq Object's Tree
#'
#' This function processes the phylogenetic tree within a `phyloseq` object by
#' resolving polytomous branching in the QIIME2 FastTree2 phylogeny, converting
#' it into a fully bifurcated (binary) tree using the `ape` package for
#' phylogenetic analysis. This step is essential for accurate evolutionary
#' analysis, and the updated `phyloseq` object is saved as an RDS file.
#'
#' @inheritParams tax_clean
#'
#' @details
#' This function performs the following steps:
#' \itemize{
#'   \item Checks if the `phyloseq` object contains a binary phylogenetic tree.
#'   \item Resolves polychotomous nodes (if present) using `ape::multi2di()`.
#'   \item Ensures that the tree is binary after resolution; raises an error if unresolved.
#'   \item Merges the resolved tree back with the `otu_table`, `sample_data`, and `tax_table` in the `phyloseq` object.
#'   \item Saves the updated `phyloseq` object with the resolved tree as an RDS file.
#' }
#'
#' @return A `phyloseq` object with a binary phylogenetic tree.
#' The `phyloseq` with a binary phylogenetic tree is saved as an RDS file named `<project_name>_phyloseq_resolved_tree.rds` in the `output_data/rds_files/Before_cleaning_rds_files` directory.
#'
#' @examples
#' \dontrun{
#' # Resolve tree and save the updated phyloseq object
#' resolved_physeq <- resolve_tree(physeq = cleaned_physeq)
#' }
#'
#' @export
resolve_tree = function(physeq, project_id, base_path, log_file) {

  log_message("Starting phylogenetic tree validation and node resolution.", status = "start", log_file)

  # Define storage directory path
  raw_rds_folder <- file.path(base_path, "raw_rds")

  # Extract tree slot
  current_tree <- phyloseq::phy_tree(physeq, errorIfNULL = FALSE)

  # Skip tree resolution if slot is empty
  if (is.null(current_tree)) {
    log_message("No phylogenetic tree found in the phyloseq object. Skipping tree resolution.", status = "info", log_file)
    return(physeq)
  } else {
    # Check if the tree is binary
    if (!ape::is.binary(current_tree)) {
      # Resolve polychotomous nodes
      phy_tree_resolved <- ape::multi2di(current_tree)

      # Check if resolved
      if (!ape::is.binary(phy_tree_resolved)) {
        error_message <- "Error: Unable to resolve polychotomous tree nodes."
        log_message(error_message, status = "error", log_file)
        stop(error_message, call. = FALSE)
      }
      binary_tree <- phy_tree_resolved
    } else {
      log_message("Phylogenetic tree is already binary. No resolution needed.", status = "info", log_file)
      binary_tree <- current_tree
    }
    # Reconstruct phyloseq object with the updated binary tree
    phyloseq_tree_resolved <- phyloseq::merge_phyloseq(
      phyloseq::otu_table(physeq),
      phyloseq::sample_data(physeq),
      phyloseq::tax_table(physeq),
      binary_tree
    )

    # Add sampleIDs
    phyloseq::sample_data(phyloseq_tree_resolved)$sample_id <- phyloseq::sample_names(phyloseq_tree_resolved)

    # Save the resolved object
    output_file_path = file.path(raw_rds_folder, glue::glue("{project_id}_phyloseq_tree_resolved.rds"))
    saveRDS(phyloseq_tree_resolved, file = output_file_path)

    # Log successful completion boundary
    log_message(glue::glue("Phyloseq object with resolved tree saved at: {output_file_path}"), status = "info", log_file)
    log_message("Phylogenetic tree resolution complete.", status = "success", log_file)

    return(phyloseq_tree_resolved)
  }
}
