# Remove Mock Features from a Phyloseq Object

This function removes mock features from a `phyloseq` object. These mock
features, which can appear in other samples due to cross-contamination,
are removed to minimize their impact on the analysis samples. In
addition, the function filters the dataset to retain only samples
without controls. Users can choose whether to remove the mock features
by setting the `mock` parameter.

## Usage

``` r
remove_mock(physeq = decontam_physeq, mock_genera = mock_genera, mock = TRUE)
```

## Arguments

- physeq:

  A `phyloseq` object containing the microbiome data. This is the input
  object that the function processes.

- mock_genera:

  A vector of genera representing the taxa that make up the mock
  community. These taxa are used to identify mock features to be removed
  from the `phyloseq` object.

- mock:

  A logical value determining whether to filter out mock features.

  - `TRUE`: Remove mock features from the `phyloseq` object and retain
    only samples for downstream analysis.

  - `FALSE`: Retain mock features. Use this option if no mock community
    is present in the dataset.

## Value

A filtered `phyloseq` object is returned and saved as an RDS file named
`project_name_phyloseq_asv_level_without_mock.rds` in the
`output_data/rds_files/Before_cleaning_rds_files/` directory.

## Details

The function performs the following steps:

- If `mock = FALSE`, the function filters the dataset to retain only
  samples without controls, leaving the mock features intact. This
  option is suitable for datasets where no mock community is included.

- If `mock = TRUE`, the function:

  - Identifies mock features based on the provided `mock_genera`.

  - Removes the mock features from the dataset.

  - Retains only samples without controls.

## Examples

``` r
if (FALSE) { # \dontrun{
# Remove mock ASVs from the phyloseq object
physeq_no_mock <- remove_mock(physeq = physeq, mock_genera = c("Mock_Genus1", "Mock_Genus2"), mock = TRUE)

# Retain mock ASVs but filter to only samples (or when no mock community is present)
physeq_no_filter <- remove_mock(physeq = physeq, mock_genera = c("Mock_Genus1", "Mock_Genus2"), mock = FALSE)
} # }
```
