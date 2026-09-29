# Aggregate Phyloseq Data by Taxonomic Level

This function aggregates ASV data at specified taxonomic levels (e.g.,
Phylum, Class, Order, Family, or Genus) using the `tax_glom` function
from the phyloseq package.

## Usage

``` r
physeq_to_taxa_tibbles(
  physeq,
  norm_method = NULL,
  copy_correction = TRUE,
  taxrank = c("Phylum", "Class", "Order", "Family", "Genus"),
  project_id,
  base_path,
  log_file
)
```

## Arguments

- taxrank:

  A character vector indicating the taxonomic levels at which to group
  the data.

## Value

The function saves multiple `phyloseq` objects as RDS files. The
aggregated objects are saved in the output directory
`output_data/rds_files/After_cleaning_rds_files/`.

## Details

The function applies the `tax_glom` function to group ASVs at each
specified taxonomic level. It creates a dedicated folder for each
taxonomic level under the output directory and saves the aggregated data
as RDS files.

## Examples

``` r
if (FALSE) { # \dontrun{
# Aggregate data using flow cytometry normalization
result <- group_tax(physeq = rarefied_asv_physeq, norm_method = "fcm")

# Aggregate data using qPCR normalization
result <- group_tax(physeq = rarefied_asv_physeq, norm_method = "qpcr")
} # }
```
