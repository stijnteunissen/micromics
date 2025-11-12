# Generate and Save Rarefaction Curve

This function creates a rarefaction curve for a given phyloseq object
and saves the plot as a PDF.

## Usage

``` r
rarefaction_curve(physeq = resolved_tree_physeq, color = NULL)
```

## Arguments

- physeq:

  A `phyloseq` object containing the microbiome data. This is the input
  object that the function processes.

- color:

  A character string specifying the column in the sample metadata to use
  for coloring the samples. Default is `NULL`, which automatically sets
  the color to `"sample_or_control"`.

## Value

The rarefaction curve plot object.

## Details

This function first checks whether the `sample_or_control` column exists
in the sample metadata. It then generates a rarefaction curve using the
`amp_rarecurve` function and saves the plot as a PDF file in the
specified directory.

## Examples

``` r
rarefaction_curve(physeq = physeq)
#> Error in log_message(paste("Step 6: Creating rarefaction curve: creating rarefaction curve before cleaning on ASV level.",     paste(projects, collapse = ", ")), log_file): could not find function "log_message"
```
