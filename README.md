
<!-- README.md is generated from README.Rmd. Please edit that file -->

# regutools

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/ComunidadBioInfo/regutools.svg?branch=master)](https://travis-ci.org/ComunidadBioInfo/regutools)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Codecov test
coverage](https://codecov.io/gh/ComunidadBioInfo/regutools/branch/master/graphs/badge.svg)](https://codecov.io/gh/ComunidadBioInfo/regutools?branch=master)
<!-- badges: end -->

The goal of `regutools` is to provide an R interface for extracting and
processing data from [RegulonDB](http://regulondb.ccg.unam.mx/).

## Installation

You can install the released version of `regutools` from
[Bioconductor](http://bioconductor.org/) with:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("regutools")

## Check that you have a valid Bioconductor installation
BiocManager::valid()
```

And the development version from [GitHub](https://github.com/) with:

``` r
BiocManager::install("ComunidadBioinfo/regutools")
```

## Example

This is a basic example which shows you how to use `regutools`. For more
details, please check the vignette.

``` r
library('regutools')
## basic example code

## Connect to the RegulonDB database if necessary
if(!exists('regulondb_conn')) regulondb_conn <- connect_database()

## Build a regulondb object
e_coli_regulondb <-
    regulondb(
        database_conn = regulondb_conn,
        organism = "E.coli",
        database_version = "1",
        genome_version = "1"
    )

## Get the araC regulators
araC_regulation <-
    get_gene_regulators(
        e_coli_regulondb,
        genes = c("araC"),
        format = "multirow",
        output.type = "TF"
    )

## Summarize the araC regulation
get_regulatory_summary(e_coli_regulondb, araC_regulation)
#> regulondb_result with 3 rows and 7 columns
#>         TF Regulated_genes_per_TF          Percent Activator Repressor     Dual
#>   <factor>               <factor>         <factor>  <factor>  <factor> <factor>
#> 1     AraC                      1 33.3333333333333         0         0        1
#> 2      CRP                      1 33.3333333333333         1         0        0
#> 3     XylR                      1 33.3333333333333         0         1        0
#>   Regulated_genes
#>          <factor>
#> 1            araC
#> 2            araC
#> 3            araC
```

# Citation

Below is the citation output from using `citation('regutools')` in R.
Please run this yourself to check for any updates on how to cite
**regutools**.

``` r
citation('regutools')
#> 
#> To cite package 'regutools' in publications use:
#> 
#>   Joselyn Chavez, Carmina Barberena-Jonas, Jesus E. Sotelo-Fonseca,
#>   Jose Alquicira-Hernandez, Heladia Salgado, Alejandro Reyes and
#>   Leonardo Collado-Torres (2020). regutools: regutools: an R package
#>   for data extraction from RegulonDB. R package version 0.99.0.
#>   https://github.com/ComunidadBioInfo/regutools
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {regutools: regutools: an R package for data extraction from RegulonDB},
#>     author = {Joselyn Chavez and Carmina Barberena-Jonas and Jesus E. Sotelo-Fonseca and Jose Alquicira-Hernandez and Heladia Salgado and Alejandro Reyes and Leonardo Collado-Torres},
#>     year = {2020},
#>     note = {R package version 0.99.0},
#>     url = {https://github.com/ComunidadBioInfo/regutools},
#>   }
```

# Testing

Testing on Bioc-devel is feasible thanks to [R
Travis](http://docs.travis-ci.com/user/languages/r/).
