# Unify and Format Metadata

This function merges and formats metadata from various sources,
including QIIME metadata, experimental sample metadata, and qPCR or FCM
data, to create a unified metadata file for downstream analyses.

## Usage

``` r
unify_metadata(projects)
```

## Arguments

- projects:

  A character vector containing the names of the project (folders).

## Value

A data frame containing the unified metadata. The data frame is also
saved as a file in the project's `input_data` folder.

## Details

The function ensures that the metadata is unified and correctly
formatted for further analysis by:

- Reading and processing the `metadata_extra.tsv` file, which must
  contain at least:

  - `SampleID`: A unique identifier for each sample.

  - `sample_type`: Indicates whether the sample is a `sample`, `mock`,
    or `blank`.

  - `DNA_Concentration`: The DNA concentration (in ng/µl).

- Optionally processing and integrating qPCR or FCM data:

  - If qPCR data is available, calculating the mean for duplicates and
    merging it with `metadata_extra`.

  - If FCM data is available, calculating the mean for duplicates and
    merging it with `metadata_extra`.

- Reading the `metadata.tsv` file (QIIME metadata), ensuring it contains
  the `SampleID` column, and combining it with the processed
  `metadata_extra`.

- Writing the final combined metadata to a file. Note that the output
  file is named by concatenating the project name with
  `_metadata_formatted.tsv` and is saved in the project's `input_data`
  folder.

All metadata files (QIIME metadata, experimental sample metadata, and
qPCR/FCM data) must include the `SampleID` column for proper merging.
This column serves as the key to align data from multiple sources.

## Note

This function requires that the folder structure has been set up (using
the `create_folders` function) before running.

## Examples

``` r
if (FALSE) { # \dontrun{
# Process and unify metadata for a project
unified_metadata <- unify_metadata(projects)
} # }
```
