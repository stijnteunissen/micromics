# Aggregate Phyloseq Data by Taxonomic Level

This function aggregates ASV data at specified taxonomic levels (e.g.,
Phylum, Class, Order, Family, or Genus) using the `tax_glom` function
from the phyloseq package.

## Usage

``` r
group_tax(
  physeq = rarefied_asv_physeq,
  norm_method = NULL,
  taxrank = c("Phylum", "Class", "Order", "Family", "Genus"),
  copy_correction = TRUE
)
```

## Arguments

- physeq:

  A `phyloseq` object containing the microbiome data. This is the input
  object that the function processes.

- norm_method:

  A character string specifying the normalization method. Acceptable
  values are:

  - `NULL`: Use this option if no FCM or qPCR data is available, or if
    you wish to retain only relative abundances.

  - `"fcm"`: Use this option if the data have been normalized using flow
    cytometry (FCM).

  - `"qpcr"`: Use this option if the data have been normalized using
    quantitative PCR (qPCR).

- taxrank:

  A character vector indicating the taxonomic levels at which to group
  the data.

- copy_correction:

  A logical value indicating whether the data should be corrected for
  the predicted 16S rRNA copy numbers prior to biomass normalization.
  Options are:

  - `TRUE`: Both relative and absolute abundances are corrected using
    the predicted copy numbers.

  - `FALSE`: Abundances are not corrected by copy number. Note that qPCR
    normalization requires copy number correction to provide absolute
    data.

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
