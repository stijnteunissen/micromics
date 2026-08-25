#' Clean Taxonomy Table in Phyloseq Object
#'
#' This function cleans and filters the taxonomy table within a `phyloseq` object.
#' It removes unclassified or ambiguous names, and replaces these and missing taxon names at genus level with placeholders derived from higher
#' taxonomic ranks. Optionally  this function filters specific taxa (e.g., Eukaryota, chloroplasts, mitochondria, and ASVs unclassified at Kingdom and Phylum levels).
#' The cleaned taxonomy enables taxonomic agglomeration (`tax_glom` form the package `phyloseq`) to genus level without merging unclassified taxa from diverse ancestry.
#'
#' @param physeq A `phyloseq` object containing the microbiome data.
#'               This is the input object that the function processes.
#'
#' @param tax_filter Iindicating whether to apply additional
#'                   filtering to remove unwanted taxa (e.g., chloroplasts, mitochondria).
#'                   \itemize{
#'                   \item `TRUE`: Apply filtering to remove taxa such as Eukaryota, chloroplasts, mitochondria, and unclassified taxa.
#'                   \item `FALSE`: Skip filtering; the taxonomy table will only be cleaned but no taxa will be removed.
#'                   }
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#' \item Cleans the taxonomy table:
#'   \itemize{
#'   \item Replaces ambiguous or placeholder taxa names (e.g., "uncultured organism") with `NA`.
#'   \item Assigns missing taxon names based on their higher-level ranks (e.g., "Phylum of Kingdom").
#'   \item Cleans the taxonomic ranks to remove overly nested placeholders.
#'   \item Ensures consistent naming for `Kingdom` (e.g., replacing "d__Bacteria" with "Bacteria").
#'   }
#' \item If `tax_filter = TRUE`, filters out unwanted taxa such as:
#'   \itemize{
#'   \item Chloroplasts (at `Class` or `Order` level).
#'   \item Mitochondria (at `Family` level).
#'   \item Unassigned taxa or taxa classified as `Eukaryota`.
#'   }
#' \item Saves the cleaned `phyloseq` object as an RDS file.
#' \item Logs the number of ASVs (amplicon sequence variants) removed during cleaning.
#' }
#'
#' The cleaned `phyloseq` object is returned for further analysis.
#'
#' @return A cleaned `phyloseq` object with ambiguous and unwanted taxa removed or replaced.
#' The cleaned `phyloseq` object is saved as an RDS file named `<project_name>_phyloseq_cleaned.rds` in the `output_data/rds_files/Before_cleaning_rds_files` directory.
#'
#' @examples
#' \dontrun{
#' # Clean taxonomy table and apply filtering
#' physeq_cleaned <- tax_clean(physeq = physeq, tax_filter = TRUE)
#'
#' # Clean taxonomy table without filtering
#' physeq_cleaned_no_filter <- tax_clean(physeq = physeq, tax_filter = FALSE)
#' }
#'
#' @export
tax_clean = function(physeq, tax_filter = TRUE, project_id, base_path, log_file) {

  log_message("Starting taxonomy table cleaning", status = "start", log_file)

  # Define storage directory path
  raw_rds_folder <- file.path(base_path, "raw_rds")

  # Track initial ASV counts
  ntaxa_in <- phyloseq::ntaxa(physeq)

  # Extract taxonomy table into a clean data frame
  tax_table_clean <- data.frame(phyloseq::tax_table(physeq))

  # Replace ambiguous strings with explicit NA values and roll forward taxonomy names
  tax_table_clean2 <- tax_table_clean %>%
    dplyr::mutate(dplyr::across(where(is.factor), as.character)) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ stringr::str_replace_all(
      .x,
      "Incertae_Sedis|Ambiguous_taxa|metagenome|uncultured archeaon|uncultured bacterium|uncultured prokaryote|uncultured soil bacterium|uncultured rumen bacterium|uncultured compost bacterium|uncultured organism|^uncultured|uncultured$",
      NA_character_
    ))) %>%
    dplyr::mutate(
      Phylum  = dplyr::if_else(is.na(Phylum), paste0("Phylum of ", Kingdom), Phylum),
      Class   = dplyr::if_else(is.na(Class), paste0("Class of ", Phylum), Class),
      Order   = dplyr::if_else(is.na(Order), paste0("Order of ", Class), Order),
      Family  = dplyr::if_else(is.na(Family), paste0("Family of ", Order), Family),
      Genus   = dplyr::if_else(is.na(Genus), paste0("Genus of ", Family), Genus),
      Species = dplyr::if_else(is.na(Species), paste0("Species of ", Genus), Species)
    ) %>%
    dplyr::mutate(dplyr::across(
      .cols = Kingdom:Species,
      .fns = ~ dplyr::if_else(stringr::str_detect(., "\\bof\\b.*\\bof\\b"), paste(stringr::word(., 1), stringr::word(., 2), stringr::word(., -1)), .)
    )) %>%
    dplyr::mutate(
      Kingdom = dplyr::if_else(Kingdom == "d__Bacteria", "Bacteria", Kingdom),
      Kingdom = dplyr::if_else(Kingdom == "d__Archaea", "Archaea", Kingdom)
    )

  # Handle numeric/placeholder Genus annotations conditionally
  tax_table_clean3 <- tax_table_clean2 %>%
    dplyr::rowwise() %>%
    dplyr::mutate(Genus = dplyr::case_when(
      grepl("^\\d+$|^[A-Z\\d_-]+$", Genus) & !stringr::str_detect(Genus, "[A-Za-z]{4,}") ~ {
        dplyr::case_when(
          !grepl("\\d", Family) ~ paste0("Genus of ", Family, " (", Genus, ")"),
          !grepl("\\d", Order)  ~ paste0("Genus of ", Order, " (", Genus, ")"),
          !grepl("\\d", Class)  ~ paste0("Genus of ", Class, " (", Genus, ")"),
          !grepl("\\d", Phylum) ~ paste0("Genus of ", Phylum, " (", Genus, ")"),
          TRUE                  ~ paste0("Genus (", Genus, ")")
        )
      },
      TRUE ~ Genus
    )) %>%
    dplyr::ungroup()

  # Reconvert table back into a matrix format and overwrite taxonomy slot
  tax_matrix <- as.matrix(tax_table_clean3)
  rownames(tax_matrix) <- phyloseq::taxa_names(physeq)
  phyloseq::tax_table(physeq) <- phyloseq::tax_table(tax_matrix)

  # Apply taxonomic filters to remove off-target host or unassigned reads
  if (tax_filter) {
    physeq <- physeq %>%
      phyloseq::subset_taxa(
        Class != "Chloroplast" &
          Order != "Chloroplast" &
          Family != "Mitochondria" &
          Kingdom != "d__Eukaryota" &
          Kingdom != "Unassigned" &
          Phylum != "Phylum of d__Bacteria" &
          Phylum != "Phylum of d__Archaea" &
          Phylum != "Phylum of Bacteria"
      )
  }

  # Calculate the absolute number of filtered ASVs
  removed_asv_count <- ntaxa_in - phyloseq::ntaxa(physeq)

  # Save the clean taxonomy object
  output_file_path <- file.path(raw_rds_folder, glue::glue("{project_id}_phyloseq_cleaned.rds"))
  saveRDS(physeq, file = output_file_path)

  # Log messages
  log_message(glue::glue("Taxonomy table cleaned: {removed_asv_count} ASVs removed."), status = "info", log_file)
  log_message(glue::glue("Cleaned phyloseq object stored at: {output_file_path}"), status = "info", log_file)
  log_message("Taxonomy table standardization complete.", status = "success", log_file)

  return(physeq)
}
