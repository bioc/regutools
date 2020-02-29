
<!-- README.md is generated from README.Rmd. Please edit that file -->

# regutools <img src="man/figures/logo.png" align="right" width="150px"/>

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/ComunidadBioInfo/regutools.svg?branch=master)](https://travis-ci.org/ComunidadBioInfo/regutools)
[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Codecov test
coverage](https://codecov.io/gh/ComunidadBioInfo/regutools/branch/master/graphs/badge.svg)](https://codecov.io/gh/ComunidadBioInfo/regutools?branch=master)
[![BioC
status](http://www.bioconductor.org/shields/build/release/bioc/regutools.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/regutools)
<!-- badges: end -->

The goal of `regutools` is to provide an R interface for extracting and
processing data from [RegulonDB](http://regulondb.ccg.unam.mx/). This
package was created as a collaboration by members of the [Community of
Bioinformatics Software Developers](https://comunidadbioinfo.github.io/)
(CDSB in Spanish).

For more details, please check the [documentation
website](http://comunidadbioinfo.github.io/regutools) or the
Bioconductor package landing page
[here](https://bioconductor.org/packages/regutools).

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
#> snapshotDate(): 2019-10-29

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
#> Barberena-Jonas C, Chavez J, Sotelo-Fonseca JE, Alquicira-Hernandez J,
#> Salgado H, Collado-Torres L, Reyes A (2020). _regutools: an R package
#> for data extraction from RegulonDB_. doi: 10.18129/B9.bioc.regutools
#> (URL: https://doi.org/10.18129/B9.bioc.regutools),
#> https://github.com/comunidadbioinfo/regutools - R package version
#> 0.99.0, <URL: http://www.bioconductor.org/packages/regutools>.
#> 
#> Barberena-Jonas C, Chavez J, Sotelo-Fonseca JE, Alquicira-Hernandez J,
#> Salgado H, Collado-Torres L, Reyes A (2020). "Programmatic access to
#> bacterial regulatory networks with regutools." _bioRxiv_. doi:
#> 10.1101/xxxyyy (URL: https://doi.org/10.1101/xxxyyy), <URL:
#> https://doi.org/10.1101/xxxyyy>.
#> 
#> To see these entries in BibTeX format, use 'print(<citation>,
#> bibtex=TRUE)', 'toBibtex(.)', or set
#> 'options(citation.bibtex.max=999)'.
```

# Development tools

  - Testing on Bioc-devel is possible thanks to [R
    Travis](http://docs.travis-ci.com/user/languages/r/).
  - Code coverage assessment is possible thanks to
    [codecov](https://codecov.io/gh).
  - The [documentation
    website](http://comunidadbioinfo.github.io/regutools) is
    automatically updated thanks to
    *[pkgdown](https://CRAN.R-project.org/package=pkgdown)* and
    *[travis](https://github.com/ropenscilabs/travis)*.

# RegulonDB

<a href="http://regulondb.ccg.unam.mx/"><img src="http://regulondb.ccg.unam.mx/img/logo.jpg"></a>

# A CDSB community project

<a href="https://comunidadbioinfo.github.io/"><img src="https://comunidadbioinfo.github.io/img/Logo_texto-768x107.png"></a>
