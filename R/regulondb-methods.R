#' @exportMethod
setMethod(
  "show",
  signature(object = "regulondb"),
  function(object) {
    cat(sprintf("regulondb object\n  organism: %s\n  genome_build: %s\n  database_version: %s\n",
                object@organism, object@genome_version,
                object@database_version))
  }
)
