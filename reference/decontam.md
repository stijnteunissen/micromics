# Decontaminate a Phyloseq Object Using Specified Methods

This function removes contamination from a phyloseq object using the
[`decontam`](https://benjjneb.github.io/decontam/vignettes/decontam_intro.html)
package. It supports the frequency, prevalence, or both methods for
contaminant identification. Note that the prevalence method relies on
the presence of blank samples; if blank samples are not available, only
the frequency method can be performed.

## Usage

``` r
decontam(
  physeq = resolved_tree_physeq,
  decon_method = c("frequency", "prevalence", "both"),
  blank = TRUE
)
```

## Arguments

- physeq:

  A `phyloseq` object containing the microbiome data. This is the input
  object that the function processes.

- decon_method:

  A character string specifying the contamination removal method.
  Possible values are:

  - `frequency`: Identifies contaminants by examining the distribution
    of sequence feature frequencies as a function of the input DNA
    concentration.

  - `prevalence`: Identifies contaminants by comparing the prevalence
    (presence/absence across samples) of sequence features in true
    samples versus negative controls (blanks).

  - `both`: Applies both frequency and prevalence methods sequentially.
    Taxa flagged by either method are considered contaminants. This
    option requires that blank samples are available.

- blank:

  A logical value indicating whether blank samples were included in the
  dataset.

  - `TRUE`: Blank samples are present, allowing the use of the `both`
    method.

  - `FALSE`: No blank samples are available; in this case, only the
    `frequency` method can be applied.

## Value

A phyloseq object with contaminants removed. The decontaminated object
is saved as an RDS file named
`<project_name>_phyloseq_asv_level_decontam.rds` in the
`output_data/rds_files/Before_cleaning_rds_files` directory.

## Details

The function uses the `decontam` package to remove contaminants based on
the specified method. If `both` is chosen, the function applies the
frequency method first and then the prevalence method, flagging any taxa
identified by either method as contaminants. In addition, diagnostic
plots (showing read counts and contaminant prevalence) are generated and
saved as PDF files.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example usage:
decontam(physeq = phyloseq_data, decon_method = "both", blank = TRUE)
} # }
```
