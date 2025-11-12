# Generate Beta Diversity Ordination Plots

This function computes beta diversity using ordination methods (e.g.,
PCoA) on a phyloseq object and generates corresponding ordination plots.
It supports both ASV-level data and data aggregated at various taxonomic
levels, and it can handle both relative and, if available, absolute
abundance data. Four distance metrics are used for relative abundance
plots (Jaccard, Bray-Curtis, Unweighted UniFrac, and Weighted UniFrac),
and Manhattan distance is used for absolute abundance plots when a
normalization method is provided.

## Usage

``` r
beta_diversity(
  physeq = physeq,
  taxrank = c("Phylum", "Class", "Order", "Family", "Genus"),
  norm_method = NULL,
  ordination_method = "PCoA",
  color_factor = NULL,
  color_continuous = TRUE,
  shape_factor = NULL,
  size_factor = NULL,
  alpha_factor = NULL
)
```

## Arguments

- physeq:

  A phyloseq object containing microbial community data.

- taxrank:

  A character vector specifying the taxonomic levels to process. If the
  first element (case-insensitive) is `"asv"`, ASV-level beta diversity
  is computed; otherwise, beta diversity is computed for each taxonomic
  level provided (default:
  `c("Phylum", "Class", "Order", "Family", "Genus")`).

- norm_method:

  A character string indicating the normalization method to be used for
  generating absolute abundance data. Options include:

  - `"fcm"`: Flow cytometry normalization.

  - `"qpcr"`: qPCR normalization.

  - `NULL`: No absolute abundance processing (default).

- ordination_method:

  A character string specifying the ordination method to use (default is
  `"PCoA"`).

- color_factor:

  An optional character string specifying the sample metadata column
  used to color the points in the ordination plot.

- color_continuous:

  A logical value indicating whether the color scale should be
  continuous (`TRUE`) or discrete (`FALSE`). Default is `TRUE`.

- shape_factor:

  An optional character string specifying the sample metadata column
  used to assign point shapes in the ordination plot.

- size_factor:

  An optional character string specifying the sample metadata column
  used to assign point sizes in the ordination plot.

- alpha_factor:

  An optional character string specifying the sample metadata column
  used to assign point transparency in the ordination plot.

## Value

A ggplot object representing the combined beta diversity ordination plot
for the relative (and, if applicable, absolute) data.

## Details

The function performs the following steps:

1.  Sets up project directories and folder paths for saving output
    files.

2.  Defines an internal helper function, `base_beta_plot`, which
    performs the ordination (using the specified `ordination_method` and
    distance metric), removes default point layers, and then adds a
    customized geom_point layer with aesthetics defined by the provided
    factors.

3.  For ASV-level data (when `taxrank[1]` is `"asv"`):

    - Processes relative abundance data (transformed to percentages)
      and, if a normalization method is specified, also processes
      absolute abundance data.

    - Generates ordination plots using multiple distance metrics:

      - **Jaccard**: Based on binary presence/absence.

      - **Bray-Curtis**: Incorporates both presence and abundance.

      - **Unweighted UniFrac**: Considers lineage presence only.

      - **Weighted UniFrac**: Considers both lineage presence and
        abundance.

    - If both DNA and RNA data are present, separate plots are generated
      for each.

4.  For other taxonomic levels:

    - Similar processing is performed for each taxonomic rank, with
      output saved in dedicated subfolders.

5.  Aesthetic scales (colors, shapes, sizes, alpha) are defined based on
    the unique levels in the sample data.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Example: Generate beta diversity plots at the ASV level using PCoA ordination and flow cytometry normalization
  beta_plot <- beta_diversity(
    physeq = my_physeq,
    taxrank = "asv",
    norm_method = "fcm",
    ordination_method = "PCoA",
    color_factor = "Treatment",
    shape_factor = "Replica",
    size_factor = "Timepoint"
  )

  # Example: Generate beta diversity plots for Phylum and Class levels without absolute data processing
  beta_plot <- beta_diversity(
    physeq = my_physeq,
    taxrank = c("Phylum", "Class"),
    norm_method = NULL,
    ordination_method = "PCoA",
    color_factor = "Soil_Type"
  )
} # }
```
