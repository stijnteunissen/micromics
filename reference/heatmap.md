# Generate a Heatmap of Relative Abundance

This function creates a heatmap of relative abundance data from a
phyloseq object (or similar data frame) at the genus level. It
calculates the relative abundance of each taxon per sample and groups
taxa with low relative abundance (below a defined threshold) into an
"Other" category. The heatmap is then facetted based on additional
sample metadata if available and saved as a PDF.

## Usage

``` r
heatmap(
  physeq = rarefied_genus_psmelt,
  ntaxa = NULL,
  norm_method = NULL,
  taxrank = c("Phylum", "Class", "Order", "Family", "Genus"),
  date_factor = NULL,
  SampleID_xas = FALSE
)
```

## Arguments

- physeq:

  A phyloseq object containing normalized genus-level data. The default
  is `rarefied_genus_psmelt`.

- ntaxa:

  An integer specifying the maximum number of taxa to display
  individually. Taxa below the threshold are grouped into "Other". If
  `NULL`, `ntaxa` is set to 23.

- norm_method:

  A character string specifying the normalization method. If `NULL`, the
  function uses the provided `physeq` directly. If set to `"fcm"` or
  `"qpcr"`, the function extracts the corresponding
  `psmelt_copy_number_corrected_` data based on the taxonomic rank.

- taxrank:

  A character string indicating the taxonomic rank to use for grouping
  taxa. The default is `"Genus"`.

## Value

A `ggplot` object representing the heatmap of relative abundance.

## Details

The function performs the following steps:

1.  Sets up project folder paths for figures and output data.

2.  Extracts and processes the input data to compute the relative
    abundance (in percentage) of each taxon per sample.

3.  Groups taxa with a mean relative abundance below a defined cutoff
    into an "Other" category.

4.  Optionally orders the data by `Sample_Date` if that factor is
    present in the metadata.

5.  Creates a base heatmap using `ggplot2`, with samples on the x-axis
    and taxa on the y-axis. The fill color reflects the relative
    abundance, and text labels are added for values exceeding a
    threshold.

6.  If more than one `na_type` is present (e.g., both DNA and RNA),
    separate heatmaps are generated for each and then combined.

7.  Saves the final heatmap as a PDF file in the project's figures
    folder.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Generate a heatmap using default parameters
  heatmap_plot <- heatmap(physeq = rarefied_genus_psmelt)

  # Generate a heatmap with a specified number of taxa and a normalization method
  heatmap_plot <- heatmap(
    physeq = rarefied_genus_psmelt,
    ntaxa = 20,
    norm_method = "fcm",
    taxrank = c("Phylum", "Class", "Order", "Family", "Genus")
  )
} # }
```
