#' Rarefy counts
#'
#' This function performs rarefaction on a normalized `phyloseq` object.
#' Rarefaction is based on the biomass of each sample to identify the minimum
#' sampling depth across the dataset. By rarefying all samples to this sampling
#' depth, we ensure that the data are normalized for sequencing effort. This
#' step prevents an overestimation of genus abundance in samples with higher
#' sequencing depths, allowing for more accurate and comparable results.
#'
#' @inheritParams normalise_data
#'
#' @param norm_method A character string specifying the normalization method. Acceptable values are:
#'   \itemize{
#'     \item `NULL`: Use this option if no FCM or qPCR data is available, or if you wish to retain only relative abundances.
#'     \item `"fcm"`: Use this option if the data have been normalized using flow cytometry (FCM).
#'     \item `"qpcr"`: Use this option if the data have been normalized using quantitative PCR (qPCR).
#'   }
#'
#' @details
#' - For `"fcm"` normalization:
#'   - Rarefies only the `fcm` normalized data, while the `copy number corrected` data remains unchanged.
#'   - Rarefies based on the calculated sampling depth, derived from the ratio of total reads per sample to the
#'     estimated cell counts.
#'   - Uses a custom function (`avgrarefy`) to perform multiple iterations of rarefaction and averages the results.
#' - For `"qpcr"` normalization:
#'   - Rarefies only the `qpcr` normalized data, while the `copy number corrected` data remains unchanged.
#'   - Rarefies based on the calculated sampling depth, derived from the ratio of total reads per sample to
#'     predicted 16S rRNA gene copy numbers.
#'   - Uses the same `avgrarefy` function for averaging rarefied counts.
#'
#' - The input values are subject to a scaling limit of 1e7. If input values exceed this limit due to the prior normalization,
#'   all values are scaled down by a calculated scaling factor. The rarefied counts are rescaled back to their original scale
#'   after rarefaction. However, due to this scaling and the rounding of scaled values, slight variations in the rarefied counts
#'   may occur.
#'
#' The function uses parallel processing to improve the efficiency of rarefaction. It saves the rarefied `phyloseq`
#' object as an `.rds` file in the appropriate output folder.
#'
#' @return The rarefied `phyloseq` object is returned and also saved as an `.rds` file. The file name and location
#' depend on the specified normalization method:
#' - For `"fcm"`: `"project_name_phyloseq_asv_level_fcm_normalised_cell_concentration_rarefied.rds"`
#' - For `"qpcr"`: `"project_name_phyloseq_asv_level_qpcr_normalised_cell_concentration_rarefied.rds"`
#'
#' @examples
#' # Rarefy using FCM normalization
#' rarefied_physeq <- rarefy_counts(physeq = normalised_physeq, norm_method = "fcm")
#'
#' # Rarefy using qPCR normalization
#' rarefied_physeq <- rarefy_counts(physeq = normalised_physeq, norm_method = "qpcr")
#'
#' @note Ensure that the `phyloseq` object is properly normalized before applying this function.
#' Missing or invalid `rarefy_to` values will result in warnings and skipped samples.
#'
#' @export
rarefy_counts = function(physeq, norm_method = NULL, copy_correction = TRUE, iteration = 100, project_id, base_path, log_file) {

  log_message("Starting rarefying data", status = "start", log_file)

  # Define directory paths
  clean_rds_folder <- file.path(base_path, "clean_rds")

  # Ensure 'vegan' is installed on the master node; install if missing
  if (!requireNamespace("vegan", quietly = TRUE)) {
    install.packages("vegan", repos = "https://cloud.r-project.org")
  }

  # Determine number of workers
  ncores <- parallel::detectCores()
  nworkers <- max(1, ncores - 2)
  cl <- parallel::makeCluster(nworkers)

  # Prepare all workers once at the beginning
  parallel::clusterEvalQ(cl, {
    if (!requireNamespace("vegan", quietly = TRUE)) {
      install.packages("vegan", repos = "https://cloud.r-project.org", quiet = TRUE)
    }
    library(vegan)
  })

  # Internal function to calculate averaged rarefactions
  avgrarefy <- function(cl_object, x, rarefy_to, iterations, seed = 711) {
    set.seed(seed)

    # Export data and target depth to workers
    parallel::clusterExport(cl_object, varlist = c("x", "rarefy_to"), envir = environment())

    # Perform parallel rarefactions
    tablist <- parallel::parLapply(cl_object, seq_len(iterations), function(i) {
      suppressWarnings(vegan::rrarefy(x, sample = rarefy_to))
    })

    # Average the results
    afunc  <- array(unlist(tablist), c(dim(tablist[[1]]), iterations))
    output <- apply(afunc, 1:2, mean)
    return(round(output, 0))
  }

  # Process standard rarefaction
  log_message("Performing standard rarefaction.", status = "info", log_file)

  # Extract the copy number corrected phyloseq object or use direct object input
  physeq_rmp <- if (is.list(physeq)) physeq[["physeq_copy_corrected"]] else physeq

  # Convert phyloseq otu table to matrix format
  otu_matrix = as(phyloseq::otu_table(physeq_rmp), "matrix")

  if (phyloseq::taxa_are_rows(physeq_rmp)) {
    otu_matrix <- t(otu_matrix)
  }

  # Determine minimal sampling depth
  min_sample <- min(phyloseq::sample_sums(physeq_rmp))

  # rarefaction taking mean of n iterations
  rarefied_matrix <- avgrarefy(cl_object = cl, x = otu_matrix, rarefy_to = min_sample, iterations = iteration, seed = 711)

  rownames(rarefied_matrix) <- rownames(otu_matrix)  # samples
  colnames(rarefied_matrix) <- colnames(otu_matrix)  # taxa

  # Reconstruct phyloseq object with rarefied counts
  rarefied_otu_table <- phyloseq::otu_table(rarefied_matrix, taxa_are_rows = FALSE)

  if (!phyloseq::taxa_are_rows(rarefied_otu_table)) {
    rarefied_otu_table <- t(rarefied_otu_table)
  }

  physeq_rmp_rarefied <- physeq_rmp
  phyloseq::otu_table(physeq_rmp_rarefied) <- rarefied_otu_table

  # Define path names based on the copy correction status
  if (copy_correction == FALSE) {
    saveRDS(file = file.path(clean_rds_folder, glue::glue("{project_id}_phyloseq_without_copy_corrected_counts_rarefied.rds")), object = physeq_rmp_rarefied)
  } else {
    saveRDS(file = file.path(clean_rds_folder, glue::glue("{project_id}_phyloseq_copy_corrected_counts_rarefied.rds")), object = physeq_rmp_rarefied)
  }

  # Process dynamic biomass rarefaction
  if (!is.null(norm_method)) {

    log_message(glue::glue("Starting biomass rarefaction"), status = "info", log_file)

    physeq_qmp <- physeq[["physeq_biomass_normalised"]]

    # Clean empty profiles form the phyloseq object
    physeq_qmp <- phyloseq::prune_samples(phyloseq::sample_sums(physeq_qmp) > 0, physeq_qmp)        # Remove samples with zero counts
    physeq_qmp <- phyloseq::prune_taxa(rowSums(phyloseq::otu_table(physeq_qmp)) > 0, physeq_qmp)  # Remove taxa with zero counts across all samples

    # Convert phyloseq sample data to data frame
    sample_data <- data.frame(phyloseq::sample_data(physeq_qmp))

    # Convert phyloseq otu table to matrix
    otu_matrix <- as(phyloseq::otu_table(physeq_qmp), "matrix")

    if (phyloseq::taxa_are_rows(physeq_qmp)) {
      otu_matrix <- t(otu_matrix)
    }

    # calculate sample sizes and rarefy depths
    sample_size <- rowSums(otu_matrix) # total reads per sample

    if (norm_method == "fcm") {
      cell_count_table <- sample_data$scaled_cells_per_ml # extract cell counts
    } else if (norm_method == "qpcr") {
      cell_count_table <- sample_data$scaled_copies_per_ml # extract cell counts
    }

    sampling_depths <- sample_size / cell_count_table # sampling depths (total reads divided by cell count)
    minimum_sampling_depth <- min(sampling_depths) # minimum sampling depth across all samples

    # Caclusate target dephts per sample
    rarefy_to <- round(cell_count_table * minimum_sampling_depth, digits = 0) # number of reads to rarefy for each sample

    # rarefy each sample based on the calculated rarefying targets (rarfy_to)
    rarefied_matrix <- matrix(
      nrow = nrow(otu_matrix),
      ncol = ncol(otu_matrix),
      dimnames = list(rownames(otu_matrix), colnames(otu_matrix))
    )

    for (i in seq_len(nrow(otu_matrix))) {
      sample_name <- rownames(otu_matrix)[i]
      sample_counts <- otu_matrix[sample_name, , drop = FALSE]

      if (!is.na(rarefy_to[i]) && rarefy_to[i] > 0) {
        rarefied_sample <- avgrarefy(cl_object = cl, x = sample_counts, rarefy_to = rarefy_to[i], iterations = iteration, seed = 711)
        rarefied_matrix[sample_name, ] <- as.numeric(rarefied_sample)
      }
    }

    # transpose
    rarefied_matrix_t <- t(rarefied_matrix)
    colnames(rarefied_matrix_t) <- phyloseq::sample_names(physeq_qmp)
    rownames(rarefied_matrix_t) <- phyloseq::taxa_names(physeq_qmp)

    # Extract scale factor
    scale_factor_df <- data.frame(phyloseq::sample_data(physeq_qmp))
    scale_factor_df$SampleID <- rownames(scale_factor_df)

    for (i in seq_len(nrow(scale_factor_df))) {
      sample_id <- scale_factor_df$SampleID[i]
      scale_factor <- scale_factor_df$scale_factor[i]
      if (!is.na(scale_factor) && scale_factor != 1) {
        rarefied_matrix_t[sample_id, ] <- rarefied_matrix_t[sample_id, ] * scale_factor
      }
    }

    # Reconstruct biomass normalised phyloseq object
    otu_rescaled <- phyloseq::otu_table(rarefied_matrix_t, taxa_are_rows = TRUE)
    physeq_qmp_rarefied <- physeq_qmp
    phyloseq::otu_table(physeq_qmp_rarefied) <- otu_rescaled

    saveRDS(file = file.path(clean_rds_folder, glue::glue("{project_id}_phyloseq_biomass_normalised_counts_rarefied.rds")), object = physeq_qmp_rarefied)

    return(list(physeq_rmp_rarefied = physeq_rmp_rarefied, physeq_qmp_rarefied = physeq_qmp_rarefied))

  } else {
    return(list(physeq_rmp_rarefied = physeq_rmp_rarefied))
  }
  # Stop the global running cluster safely
  parallel::stopCluster(cl)
}
