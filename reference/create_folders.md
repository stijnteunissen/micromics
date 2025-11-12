# Create Project Folder Structure

This function creates a folder structure for projects, ensuring that all
necessary directories exist and that specific files required for
downstream analyses are present.

## Usage

``` r
create_folders(projects)
```

## Arguments

- projects:

  A character vector containing the names of the project (folders).

- base_path:

  A character string indicating the base directory where the project
  folders are located. This path is used to locate the project folder
  and create the required subfolders.

- log_file:

  A character string specifying the path to the log file where warnings
  and errors will be recorded.

## Value

None. This function is called for its side effects.

## Details

This function facilitates the setup of downstream analyses by:

- Creating a consistent directory structure for each project, including
  subfolders such as `input_data`, `output_data`, `figures`, and
  `messages`.

- Copying essential files from the `qiime2_output` folder into the
  `input_data` folder. These essential files include:

  - `table.qza`: The feature table from QIIME2.

  - `rooted-tree.qza`: The phylogenetic tree used for diversity
    analysis.

  - `classifier.qza`: The classifier used for taxonomy assignment.

  - `metadata.tsv`: The QIIME2 sample metadata file required for
    analyses.

  - `metadata_extra.tsv`: Any additional metadata provided for the
    samples.

  - `pantaxa_stats_NCBI.tsv`: The reference database for copy number
    correction from
    [rrndb](https://rrndb.umms.med.umich.edu/downloads/).

  - `prediction.RDS`: The predicted 16S copy numbers for each feature.

- Checking for optional files that enhance analyses, such as:

  - `qPCR.csv`: Contains quantitative PCR data.

  - `fcm.csv`: Contains flow cytometry data.

- Logging warnings for missing optional files and errors for missing
  required files.

The function ensures that downstream analyses—which rely on specific
input files (e.g., `table.qza`, `rooted-tree.qza`, etc.) have access to
these files in the correct directory structure. If any required files
are missing from the `qiime2_output` folder, the function stops
execution and logs an error message.

## Note

Each project folder must already exist within the `base_path` directory
and must contain a subfolder named `qiime2_output`, which holds the
outputs of QIIME2 analysis and other necessary files. The function sets
up the folder structure for downstream analysis within this project
folder.

## Examples

``` r
if (FALSE) { # \dontrun{
# Define the base path and log file location
projects <- "project_name"
base_path <- "path/to/projects"
log_file <- "path/to/log_file.log"

# Create folder structures for projects
create_folders(projects)
} # }
```
