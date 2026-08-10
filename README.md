# micromics

*micromics*is an R package designed to analyse and visualise microbiome data derived from *QIIME2* output files in a streamlined and reproducible manner. It simplifies the import and processing of *QIIME2* artefacts and includes a normalisation step using biomass data to convert read counts into biological estimates, such as cell equivalents per millilitre of sample. This conversion allows for more accurate interpretation and comparison of microbial abundance across different samples. 

## Installation

*micromics* is implemented in R (4.5.1). You can install *micromics* using *devtools*:

``` r
# install.packages
devootls::install_github("stijnteunissen/micromics")
```

## License

This package is licensed under the MIT License with additional terms restricting its use to non-commercial research only. See the [LICENSE](LICENSE.md) file for details.

## Citing micromics
Teunissen, S., Silvius, J., van Veelen, H. P. J. (2025). micromics: an R package to integrate QIIME2 output and biomass quantitation for pseudo-quantitative microbiome analysis, R package version 1.0.0, GitHub. https://github.com/stijnteunissen/micromics.
