#' Methods for regulondb objects
#' @name show
#' @aliases show,regulondb-method
#' @param object A regulondb object
#' @docType methods
#' @rdname regulondb-methods
#' @return A [regulondb][regutools::regulondb] object.
setMethod("show",
        signature(object = "regulondb"),
        function(object) {
            cat(
                sprintf(
                    "regulondb object\n  organism: %s\n  genome_build: %s\n  database_version: %s\n  database_conn: %s\n",
                    object@organism,
                    object@genome_version,
                    object@database_version,
                    object@dbname
                )
            )
        })
