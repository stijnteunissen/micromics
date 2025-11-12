# Clean Taxonomy Table in Phyloseq Object

This function cleans and filters the taxonomy table within a `phyloseq`
object. It removes unclassified or ambiguous names, and replaces these
and missing taxon names at genus level with placeholders derived from
higher taxonomic ranks. Optionally this function filters specific taxa
(e.g., Eukaryota, chloroplasts, mitochondria, and ASVs unclassified at
Kingdom and Phylum levels). The cleaned taxonomy enables taxonomic
agglomeration (`tax_glom` form the package `phyloseq`) to genus level
without merging unclassified taxa from diverse ancestry.

## Usage

``` r
tax_clean(physeq = physeq, tax_filter = TRUE)
```

## Arguments

- physeq:

  A `phyloseq` object containing the microbiome data. This is the input
  object that the function processes.

- tax_filter:

  Iindicating whether to apply additional filtering to remove unwanted
  taxa (e.g., chloroplasts, mitochondria).

  - `TRUE`: Apply filtering to remove taxa such as Eukaryota,
    chloroplasts, mitochondria, and unclassified taxa.

  - `FALSE`: Skip filtering; the taxonomy table will only be cleaned but
    no taxa will be removed.

## Value

A cleaned `phyloseq` object with ambiguous and unwanted taxa removed or
replaced. The cleaned `phyloseq` object is saved as an RDS file named
`<project_name>_phyloseq_cleaned.rds` in the
`output_data/rds_files/Before_cleaning_rds_files` directory.

## Details

The function performs the following steps:

- Cleans the taxonomy table:

  - Replaces ambiguous or placeholder taxa names (e.g., "uncultured
    organism") with `NA`.

  - Assigns missing taxon names based on their higher-level ranks (e.g.,
    "Phylum of Kingdom").

  - Cleans the taxonomic ranks to remove overly nested placeholders.

  - Ensures consistent naming for `Kingdom` (e.g., replacing
    "d\_\_Bacteria" with "Bacteria").

- If `tax_filter = TRUE`, filters out unwanted taxa such as:

  - Chloroplasts (at `Class` or `Order` level).

  - Mitochondria (at `Family` level).

  - Unassigned taxa or taxa classified as `Eukaryota`.

- Saves the cleaned `phyloseq` object as an RDS file.

- Logs the number of ASVs (amplicon sequence variants) removed during
  cleaning.

The cleaned `phyloseq` object is returned for further analysis.

## Examples

``` r
if (FALSE) { # \dontrun{
# Clean taxonomy table and apply filtering
physeq_cleaned <- tax_clean(physeq = physeq, tax_filter = TRUE)

# Clean taxonomy table without filtering
physeq_cleaned_no_filter <- tax_clean(physeq = physeq, tax_filter = FALSE)
} # }
```
