# Rarefying Phyloseq Data

This function performs rarefaction on a normalized `phyloseq` object.
Rarefaction is based on the biomass of each sample to identify the
minimum sampling depth across the dataset. By rarefying all samples to
this sampling depth, we ensure that the data are normalized for
sequencing effort. This step prevents an overestimation of genus
abundance in samples with higher sequencing depths, allowing for more
accurate and comparable results.

## Usage

``` r
rarefying(
  physeq = physeq,
  norm_method = NULL,
  iteration = 100,
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

The rarefied `phyloseq` object is returned and also saved as an `.rds`
file. The file name and location depend on the specified normalization
method:

- For `"fcm"`:
  `"project_name_phyloseq_asv_level_fcm_normalised_cell_concentration_rarefied.rds"`

- For `"qpcr"`:
  `"project_name_phyloseq_asv_level_qpcr_normalised_cell_concentration_rarefied.rds"`

## Details

- For `"fcm"` normalization:

  - Rarefies only the `fcm` normalized data, while the
    `copy number corrected` data remains unchanged.

  - Rarefies based on the calculated sampling depth, derived from the
    ratio of total reads per sample to the estimated cell counts.

  - Uses a custom function (`avgrarefy`) to perform multiple iterations
    of rarefaction and averages the results.

- For `"qpcr"` normalization:

  - Rarefies only the `qpcr` normalized data, while the
    `copy number corrected` data remains unchanged.

  - Rarefies based on the calculated sampling depth, derived from the
    ratio of total reads per sample to predicted 16S rRNA gene copy
    numbers.

  - Uses the same `avgrarefy` function for averaging rarefied counts.

- The input values are subject to a scaling limit of 1e7. If input
  values exceed this limit due to the prior normalization, all values
  are scaled down by a calculated scaling factor. The rarefied counts
  are rescaled back to their original scale after rarefaction. However,
  due to this scaling and the rounding of scaled values, slight
  variations in the rarefied counts may occur.

The function uses parallel processing to improve the efficiency of
rarefaction. It saves the rarefied `phyloseq` object as an `.rds` file
in the appropriate output folder.

## Note

Ensure that the `phyloseq` object is properly normalized before applying
this function. Missing or invalid `rarefy_to` values will result in
warnings and skipped samples.

## Examples

``` r
# Rarefy using FCM normalization
rarefied_physeq <- rarefying(physeq = normalised_physeq, norm_method = "fcm")
#> Error in log_message(paste("Step 10: rarefied data", paste(projects, collapse = ", ")),     log_file): could not find function "log_message"

# Rarefy using qPCR normalization
rarefied_physeq <- rarefying(physeq = normalised_physeq, norm_method = "qpcr")
#> Error in log_message(paste("Step 10: rarefied data", paste(projects, collapse = ", ")),     log_file): could not find function "log_message"
```
