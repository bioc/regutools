
<!-- README.md is generated from README.Rmd. Please edit that file -->

# regutools <img src="man/figures/logo.png" align="right" width="150px"/>

<!-- badges: start -->

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![BioC
status](http://www.bioconductor.org/shields/build/release/bioc/regutools.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/regutools)
[![Codecov test
coverage](https://codecov.io/gh/ComunidadBioInfo/regutools/branch/master/graph/badge.svg)](https://codecov.io/gh/ComunidadBioInfo/regutools?branch=master)
[![R build
status](https://github.com/ComunidadBioInfo/regutools/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/ComunidadBioInfo/regutools/actions)
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
#> snapshotDate(): 2020-03-31

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
#>            TF Regulated_genes_per_TF          Percent   Activator   Repressor
#>   <character>            <character>      <character> <character> <character>
#> 1        AraC                      1 33.3333333333333           0           0
#> 2         CRP                      1 33.3333333333333           1           0
#> 3        XylR                      1 33.3333333333333           0           1
#>          Dual Regulated_genes
#>   <character>     <character>
#> 1           1            araC
#> 2           0            araC
#> 3           0            araC
```

## Documentation

For more information about `regutools` check the vignettes [through
Bioconductor](http://bioconductor.org/packages/regutools) or at the
[documentation website](http://comunidadbioinfo.github.io/regutools).

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `regutools` using from
[Bioconductor](http://bioconductor.org/) the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("regutools")
```

## Citation

Below is the citation output from using `citation('regutools')` in R.
Please run this yourself to check for any updates on how to cite
**regutools**.

``` r
print(citation('regutools'), bibtex = TRUE)
#> 
#> Chávez J, Barberena-Jonas C, Sotelo-Fonseca JE, Alquicira-Hernandez J,
#> Salgado H, Collado-Torres L, Reyes A (2020). _regutools: an R package
#> for data extraction from RegulonDB_. doi: 10.18129/B9.bioc.regutools
#> (URL: https://doi.org/10.18129/B9.bioc.regutools),
#> https://github.com/comunidadbioinfo/regutools - R package version
#> 0.99.13, <URL: http://www.bioconductor.org/packages/regutools>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {regutools: an R package for data extraction from RegulonDB},
#>     author = {Joselyn Chávez and Carmina Barberena-Jonas and Jesus Emiliano Sotelo-Fonseca and Jose Alquicira-Hernandez and Heladia Salgado and Leonardo Collado-Torres and Alejandro Reyes},
#>     year = {2020},
#>     url = {http://www.bioconductor.org/packages/regutools},
#>     note = {https://github.com/comunidadbioinfo/regutools - R package version 0.99.13},
#>     doi = {10.18129/B9.bioc.regutools},
#>   }
#> 
#> Chávez J, Barberena-Jonas C, Sotelo-Fonseca JE, Alquicira-Hernandez J,
#> Salgado H, Collado-Torres L, Reyes A (2020). "Programmatic access to
#> bacterial regulatory networks with regutools." _bioRxiv_. doi:
#> 10.1101/xxxyyy (URL: https://doi.org/10.1101/xxxyyy), <URL:
#> https://doi.org/10.1101/xxxyyy>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {Programmatic access to bacterial regulatory networks with regutools},
#>     author = {Joselyn Chávez and Carmina Barberena-Jonas and Jesus Emiliano Sotelo-Fonseca and Jose Alquicira-Hernandez and Heladia Salgado and Leonardo Collado-Torres and Alejandro Reyes},
#>     year = {2020},
#>     journal = {bioRxiv},
#>     doi = {10.1101/xxxyyy},
#>     url = {https://doi.org/10.1101/xxxyyy},
#>   }
```

Please note that the `regutools` was only made possible thanks to many
other R and bioinformatics software authors, which are cited either in
the vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the derfinderPlot project is released with a
[Contributor Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

## Development tools

  - Continuous code testing is possible thanks to [GitHub
    actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
    through *[usethis](https://CRAN.R-project.org/package=usethis)*,
    *[remotes](https://CRAN.R-project.org/package=remotes)*,
    *[sysreqs](https://github.com/r-hub/sysreqs)* and
    *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)*
    customized to use [Bioconductor’s docker
    containers](https://www.bioconductor.org/help/docker/) and
    *[BiocCheck](https://bioconductor.org/packages/3.11/BiocCheck)*.
  - Code coverage assessment is possible thanks to
    [codecov](https://codecov.io/gh) and
    *[covr](https://CRAN.R-project.org/package=covr)*.
  - The [documentation
    website](http://comunidadbioinfo.github.io/regutools) is
    automatically updated thanks to
    *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
  - The code is styled automatically thanks to
    *[styler](https://CRAN.R-project.org/package=styler)*.
  - The documentation is formatted thanks to
    *[devtools](https://CRAN.R-project.org/package=devtools)* and
    *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

For more details, check the `dev` directory.

# RegulonDB

<a href="http://regulondb.ccg.unam.mx/"><img src="http://regulondb.ccg.unam.mx/img/logo.jpg"></a>

# A CDSB community project

<a href="https://comunidadbioinfo.github.io/"><img src="https://comunidadbioinfo.github.io/img/Logo_texto-768x107.png"></a>
