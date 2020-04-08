# regutools 0.99.8

* Change biocViews from AnnotationData to Software 


# regutools 0.99.7

BUG FIXES

* Fix a bug in pkgdown that is currently resolved by
`Rscript -e 'remotes::install_github("r-lib/pkgdown#1276")'`. 
Details at https://github.com/r-lib/pkgdown/pull/1276 and
related issues.

# regutools 0.99.6

* Manual adjust for some lines > 80 characters

# regutools 0.99.5

* Reindent lines and try line-wrapping for lines > 80 characters

# regutools 0.99.4

NEW FEATURES

* Added a `NEWS.md` file to track changes to the package.

# regutools 0.99.3

SIGNIFICANT USER-VISIBLE CHANGES

* DESCRIPTION file now contains an article-like summary in the field `regutools` package Description.
* DESCRIPTION was updated to show only one maintainer as required by http://bioconductor.org/developers/package-guidelines/#description

BUG FIXES

* Fix a bug in the function `get_regulatory_network()` related with Cytoscape conection.

# regutools 0.99.2

SIGNIFICANT CHANGES

* Removed the regutools.Rproj file (you can still have it, but it's not
version controlled).
* Bump R depends to 4.0 as required by BiocCheck. We'll see if it breaks.

# regutools 0.99.1

BUG FIXES

* Address as many errors, warnings and notes as possible from
http://bioconductor.org/spb_reports/regutools_buildreport_20200302121040.html#malbec2_check_anchor. For the Windows error, it's likely related to `mode = "w"` versus
`mode = "wb"` for `utils::download.file()` as in 
https://github.com/LieberInstitute/spatialLIBD/commit/7ac24dc72cc75e42c9af8382c661517f4a64f08d.

# regutools 0.99.0

* Submitted to Bioconductor
