# Resolve Multichotomies in Phyloseq Object's Tree

This function processes the phylogenetic tree within a `phyloseq` object
by resolving polytomous branching in the QIIME2 FastTree2 phylogeny,
converting it into a fully bifurcated (binary) tree using the `ape`
package for phylogenetic analysis. This step is essential for accurate
evolutionary analysis, and the updated `phyloseq` object is saved as an
RDS file.

## Usage

``` r
resolve_tree(physeq = cleaned_physeq)
```

## Arguments

- physeq:

  A `phyloseq` object containing the microbiome data. This is the input
  object that the function processes.

## Value

A `phyloseq` object with a binary phylogenetic tree. The `phyloseq` with
a binary phylogenetic tree is saved as an RDS file named
`<project_name>_phyloseq_resolved_tree.rds` in the
`output_data/rds_files/Before_cleaning_rds_files` directory.

## Details

This function performs the following steps:

- Checks if the `phyloseq` object contains a binary phylogenetic tree.

- Resolves polychotomous nodes (if present) using `ape::multi2di()`.

- Ensures that the tree is binary after resolution; raises an error if
  unresolved.

- Merges the resolved tree back with the `otu_table`, `sample_data`, and
  `tax_table` in the `phyloseq` object.

- Saves the updated `phyloseq` object with the resolved tree as an RDS
  file.

## Examples

``` r
if (FALSE) { # \dontrun{
# Resolve tree and save the updated phyloseq object
resolved_physeq <- resolve_tree(physeq = cleaned_physeq)
} # }
```
