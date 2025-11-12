# Create Barplots for Relative and Absolute Abundance

This function generates barplots for microbial data at the genus level.
It supports both relative and absolute abundance data and can include
facets based on available metadata factors. The resulting plots can be
saved as PDF files, and the underlying data can be exported as CSV and
RDS files.

## Usage

``` r
barplot(
  physeq = rarefied_genus_psmelt,
  ntaxa = NULL,
  norm_method = NULL,
  sample_matrix = NULL,
  group_by_factor = NULL,
  taxrank = c("Phylum", "Class", "Order", "Family", "Genus"),
  date_factor = NULL
)
```

## Arguments

- physeq:

  A `phyloseq` object containing the microbiome data. This is the input
  object that the function processes.

- ntaxa:

  An integer specifying the maximum number of taxa to display in the
  barplot. Default is 23.

- norm_method:

  A string indicating the normalization method used for absolute
  abundance data. Options are `"fcm"` (flow cytometry) or `"qpcr"`
  (quantitative PCR). (relative abundance only).

- sample_matrix:

  An optional matrix specifying the sample structure or metadata.

- group_by_factor:

  with this option you can separtate de barplot for factors

- taxrank:

  A character vector indicating the taxonomic levels at which to group
  the data.

## Value

The function generates and saves barplots as PDF files in the project’s
`figures/` folder. It also saves the processed data as CSV and RDS files
in the corresponding `output_data/` folders. Additionally, the function
outputs the plot object for further customization if needed.

## Details

- Relative abundance plots show proportions of taxa in each sample, with
  taxa having a mean relative abundance below 1% grouped as "Other".

- Absolute abundance plots use normalized cell equivalents
  (`norm_method = "fcm"` or `"qpcr"`) to display the number of cells per
  mL for each taxon.

- Facets are added based on metadata factors present in the phyloseq
  object.

- Taxa labels are styled to include genus and species names, if
  available.

## Examples

``` r
# Example usage
barplot(
  physeq = rarefied_genus_psmelt,
  ntaxa = 20,
  colorset = my_colors,
  norm_method = "fcm",
  sample_matrix = sample_metadata
)
#> Error in barplot(physeq = rarefied_genus_psmelt, ntaxa = 20, colorset = my_colors,     norm_method = "fcm", sample_matrix = sample_metadata): unused argument (colorset = my_colors)
```
