# regutools 0.99.2

* Removed the regutools.Rproj file (you can still have it, but it's not
version controlled).
* Bump R depends to 4.0 as required by BiocCheck. We'll see if it breaks.

# regutools 0.99.1

* Address as many errors, warnings and notes as possible from
http://bioconductor.org/spb_reports/regutools_buildreport_20200302121040.html#malbec2_check_anchor. For the Windows error, it's likely related to `mode = "w"` versus
`mode = "wb"` for `utils::download.file()` as in 
https://github.com/LieberInstitute/spatialLIBD/commit/7ac24dc72cc75e42c9af8382c661517f4a64f08d.

# regutools 0.99.0

* Added a `NEWS.md` file to track changes to the package.
