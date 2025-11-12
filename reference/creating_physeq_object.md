# Create a Phyloseq Object

This function creates a `phyloseq` object using input data files such as
feature tables, taxonomic assignments, phylogenetic tree, and unified
metadata into a single `phyloseq` object for downstream analysis.

## Usage

``` r
creating_physeq_object(projects)
```

## Arguments

- projects:

  A character vector containing the names of the project (folders).

## Value

A `phyloseq` object that integrates feature tables, taxonomy,
phylogenetic trees, and metadata. This object is also saved as an RDS
file for further usage in downstream analyses.

## Details

This function performs the following steps:

- Defines the paths to the required input files (feature table, rooted
  tree, taxonomy, and metadata).

- Searches for and retrieves these files from the `input_data`
  directory.

- Calls the `qza_to_phyloseq()` function from the `qiime2R` package to
  generate a `phyloseq` object based on the provided input files.

- Adds read count information to the sample metadata within the
  `phyloseq` object.

The function assumes the following files are present in the `input_data`
directory:

- `table.qza`: Feature table containing sample feature data.

- `rooted-tree.qza`: Phylogenetic tree.

- `classifier.qza`: Taxonomic classification file.

- `metadata_formatted.tsv`: Unified sample metadata.

The resulting `phyloseq` object is essential for downstream analyses and
integrates all input files into a single, structured object. The created
`phyloseq` object is saved as an RDS file named
`<project_name>_phyloseq_uncleaned.rds` in the
`output_data/rds_files/Before_cleaning_rds_files` directory.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example usage:
physeq_object <- creating_physeq_object(projects)
} # }
```
