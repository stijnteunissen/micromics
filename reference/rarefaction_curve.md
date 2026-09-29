# Generate and Save Rarefaction Curve

Generates a sequencing depth rarefaction curve from a phyloseq object
and exports the resulting plot as a PDF.

## Usage

``` r
rarefaction_curve(
  physeq,
  color = "sample_or_control",
  project_id,
  base_path,
  log_file
)
```

## Arguments

- physeq:

  A phyloseq object to be analyzed.

- color:

  Character string specifying the metadata column used for coloring
  lines. Defaults to `"sample_or_control"`.

- project_id:

  Character string specifying the unique project name.

- base_path:

  Character string specifying the root project directory.

- log_file:

  Character string specifying the path to the log file.

## Value

None. This function is called for its side effects (saving a PDF file).
