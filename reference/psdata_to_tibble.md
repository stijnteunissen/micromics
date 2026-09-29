# Convert Phyloseq Data to Tibble and Export by Taxonomic Level

This function transforms a `phyloseq` object into a tibble for each
specified taxonomic level using the `psmelt` function from the phyloseq
package.

## Usage

``` r
psdata_to_tibble(
  physeq = rarefied_genus_physeq,
  norm_method = NULL,
  taxrank = c("Phylum", "Class", "Order", "Family", "Genus")
)
```

## Value

The function saves multiple `phyloseq`-derived tibbles as RDS files and
returns a list of tibbles:

- **If `norm_method` is `NULL`:** A tibble of copy number–corrected
  counts is saved as
  `<project_name>_psmelt_<tax>_level_copy_number_corrected_counts.rds`
  for each taxonomic level.

      \item **If \code{norm_method = "fcm"}:**
            Two tibbles are saved for each taxonomic level:
            \itemize{
              \item A tibble of copy number–corrected counts.
              \item A tibble of FCM-normalized, rarefied counts saved as
                    \code{<project_name>_psmelt_<tax>_level_fcm_normalised_cell_concentration_rarefied.rds}.
            }

      \item **If \code{norm_method = "qpcr"}:**
            Two tibbles are saved for each taxonomic level:
            \itemize{
              \item A tibble of copy number–corrected counts.
              \item A tibble of qPCR-normalized, rarefied counts saved as
                    \code{<project_name>_psmelt_<tax>_level_qpcr_normalised_cell_concentration_rarefied.rds}.
            }

## Details

The primary task of this function is to convert a `phyloseq` object into
a tibble at each specified taxonomic level using the `psmelt` function.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Export data without normalization (copy number–corrected counts only)
  result <- psdata_to_tibble(physeq = rarefied_genus_physeq)

  # Export data with flow cytometry normalization
  result <- psdata_to_tibble(physeq = rarefied_genus_physeq, norm_method = "fcm")

  # Export data with qPCR normalization
  result <- psdata_to_tibble(physeq = rarefied_genus_physeq, norm_method = "qpcr")
} # }
```
