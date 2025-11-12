# Generate Alpha Diversity Plots

This function calculates alpha diversity metrics from a phyloseq object
and generates alpha diversity plots at various taxonomic levels or at
the ASV level. The alpha diversity measures (Observed, Chao1, Shannon,
Simpson) are computed using the `estimate_richness` function from the
phyloseq package. Depending on the input parameters, the function can
handle normalized data (using flow cytometry or qPCR methods) and can
separate plots by DNA and RNA types.

## Usage

``` r
alpha_diversity(
  physeq = physeq,
  norm_method = NULL,
  taxrank = c("Phylum", "Class", "Order", "Family", "Genus"),
  date_factor = NULL
)
```

## Arguments

- physeq:

  A phyloseq object containing microbial community data.

- norm_method:

  A character string specifying the normalization method for the data.
  Options include:

  - `"fcm"`: Use flow cytometry-normalized data.

  - `"qpcr"`: Use qPCR-normalized data.

  - `NULL`: Use only copy number corrected data (default).

- taxrank:

  A character vector specifying the taxonomic levels for which alpha
  diversity is to be calculated. If the first element is (taxrank =
  `asv`), ASV-level data is processed. Otherwise, the function processes
  data for each taxonomic level provided default is taxrank =
  c('Phylum', 'Class', 'Order', 'Family', 'Genus').

- date_factor:

  An optional character string indicating the name of the date column in
  the sample metadata. If provided, the column is converted to a Date
  object ("%d/%m/%Y") and used to order the data.

## Value

A combined ggplot object containing the generated alpha diversity plots.

## Details

The function performs the following steps:

1.  Extracts the appropriate data object from `physeq` based on the
    chosen `norm_method` and taxonomic level.

2.  Estimates alpha diversity metrics (Observed, Chao1, Shannon,
    Simpson) using `estimate_richness`.

3.  Merges the alpha diversity estimates with sample metadata.

4.  Optionally orders the data by a date factor if provided.

5.  Appends dummy data to the dataset for visualization purposes.

6.  Exports the combined alpha diversity data as CSV files.

7.  Generates bar plots for the diversity metrics, optionally separating
    DNA and RNA data if both are present.

8.  Saves the resulting plots as PDF files in the project's figures
    folder.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Generate alpha diversity plots at the ASV level without normalization
  alpha_plot <- alpha_diversity(physeq = my_physeq, taxrank = "asv")

  # Generate alpha diversity plots at the Phylum level using flow cytometry-normalized data
  alpha_plot <- alpha_diversity(physeq = my_physeq, norm_method = "fcm", taxrank = "Phylum")

  # Generate alpha diversity plots with a specified date factor for ordering samples
  alpha_plot <- alpha_diversity(physeq = my_physeq, date_factor = "Sample_Date")
} # }
```
