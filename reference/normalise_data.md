# Normalize Phyloseq Data

This function applies normalization to a `phyloseq` object, converting
data to absolute values based on 16S rRNA copy numbers and sample
biomass, using either flow cytometry (FCM) data or qPCR data. The
function can apply copy number correction prior to biomass
normalization.

## Usage

``` r
normalise_data(
  physeq = without_mock_physeq,
  norm_method = NULL,
  copy_correction = TRUE
)
```

## Arguments

- physeq:

  A `phyloseq` object containing the microbiome data. This is the input
  object that the function processes.

- norm_method:

  A character string specifying the normalization method. Options are:

  - `"fcm"`: Normalize based on flow cytometry data, converting
    abundances to cell concentrations (cells/mL or per gram sample).

  - `"qpcr"`: Normalize based on qPCR data, converting abundances to
    cell equivalents (cells/mL or per gram sample).

  - `NULL`: Apply only copy number correction without further
    normalization.

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

The function saves multiple `phyloseq` objects as RDS files:

- `<project_name>_phyloseq_asv_level_without_copy_number_corrected_counts.rds`:
  if `copy_correction = FALSE`; a phyloseq object with uncorrected
  counts.

- `<project_name>_phyloseq_asv_level_copy_number_corrected_counts.rds`:
  if `copy_correction = TRUE`; a phyloseq object with counts corrected
  by predicted 16S copy numbers.

- `<project_name>_phyloseq_asv_level_fcm_normalised_cell_concentration.rds`:
  if `norm_method = "fcm"` and `copy_correction = TRUE`; a phyloseq
  object with abundances normalized based on FCM data and copy number
  correction.

- `<project_name>_phyloseq_asv_level_fcm_normalised_cell_concentration_without_copy_number_corrected_count.rds`:
  if `norm_method = "fcm"` and `copy_correction = FALSE`; a phyloseq
  object with FCM normalization applied without copy number correction.

- `<project_name>_phyloseq_asv_level_qpcr_normalised_cell_concentration.rds`:
  if `norm_method = "qpcr"`; a phyloseq object with abundances
  normalized to cell equivalents using qPCR data.

The relative phyloseq object (without biomass normalization) is saved in
the `output_data/rds_files/After_cleaning_rds_files/ASV` directory, and
all biomass-normalized phyloseq objects are saved in the
`output_data/rds_files/Before_cleaning_rds_files` directory.

## Details

The function follows these steps based on the chosen parameters:

1.  **Copy Number Correction (if `copy_correction = TRUE`):**

    - Correct ASV abundances by dividing each count by its predicted 16S
      rRNA copy number. The prediction is based on the method described
      in ["Accounting for 16S rRNA copy number prediction uncertainty
      and its implications in bacterial diversity
      analyses"](https://dx.doi.org/10.1038/s43705-023-00266-0).

    - This correction adjusts for variability in 16S rRNA gene copy
      numbers across taxa, enabling the calculation of cell equivalents.

2.  **FCM Normalization (`norm_method = "fcm"`):**

    - When `copy_correction = TRUE`: The copy number–corrected
      abundances are multiplied by the FCM data, where FCM data (cells
      per mL or per gram) is included in the metadata or a file with
      "fcm" in the name, with the column `cells_per_ml`.

    - When `copy_correction = FALSE`: The raw abundances are multiplied
      by the FCM data without prior copy number correction.

3.  **qPCR Normalization (`norm_method = "qpcr"`):**

    - The qPCR data, provided in 16S copies per mL or per gram sample
      (included in the metadata or a file with "qpcr" in the name and
      column `sq_calc_mean`), is used together with copy number
      predictions to calculate absolute abundances.

### DNA vs. RNA Normalization

The interpretation of normalized data depends on the nucleic acid type:

- **DNA:** Normalized abundances usually represent **cells per mL (or
  per gram)**, assuming one genome copy per cell.

- **RNA:** Normalized abundances often represent **copies per cell
  equivalent per mL (or per gram)**; RNA reflects transcriptional
  activity and may vary considerably with cell condition.

## References

Gao, Y., & Wu, M. (2023). Accounting for 16S rRNA copy number prediction
uncertainty and its implications in bacterial diversity analyses. *ISME
Communications, 3*(1), 59.
doi:[10.1038/s43705-023-00266-0](https://dx.doi.org/10.1038/s43705-023-00266-0)

## Examples

``` r
if (FALSE) { # \dontrun{
# Apply only copy number correction
result <- normalize_data(physeq = physeq, norm_method = NULL)

# Normalize using flow cytometry (FCM) data
result <- normalize_data(physeq = physeq, norm_method = "fcm", copy_correction = TRUE)

# Normalize using qPCR data
result <- normalize_data(physeq = physeq, norm_method = "qpcr", copy_correction = TURE)
} # }
```
