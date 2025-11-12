# Export Project Data to a Zip Archive

This function exports project-related figures and output data by
creating a dedicated export folder within the project directory. It
organizes the export folder into subdirectories for figures, CSV files,
and RDS files, copies the relevant files from their original locations,
and finally compresses the export folder into a zip archive containing
only relative file paths.

## Usage

``` r
export_data()
```

## Value

This function does not return a value. It performs file operations and
creates a zip archive in the project folder.

## Details

The function performs the following steps:

1.  Constructs paths for the project folder, figures, and output data
    using global variables (e.g., `projects` and `base_path`).

2.  Creates an export folder named `export_data_<project_name>` within
    the project folder.

3.  Creates subdirectories within the export folder for figures, CSV
    files, and RDS files.

4.  Copies figures from the project's figure folder to the export
    folder.

5.  Copies CSV files from the output CSV folder to the corresponding
    subfolder in the export directory.

6.  Copies RDS files from the `After_cleaning_rds_files` subfolder (if
    present) to the export folder.

7.  Zips the contents of the export folder into a zip archive containing
    only the relative file paths.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Export project data to a zip archive
  export_data()
} # }
```
