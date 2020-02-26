#' @title Guess gene id type
#' @description Given a gene identifier, return the most likely gene_id type.
#' @param gene Character vector of gene identifiers (id, name, bnumber or gi).
#' @param regulondb A [regulondb()] object.
#' @keywords geneid, bnumber, gi, synonyms
#' @author Jesús Emiliano Sotelo Fonseca
#' @examples
#' ## Download the database if necessary
#' if(!file.exists(file.path(tempdir(), 'regulondb_sqlite3.db'))) {
#'     download_database(tempdir())
#' }
#'
#' ## Build the regulon db object
#' e_coli_regulondb <-
#'     regulondb(
#'         database_path = file.path(tempdir(), "regulondb_sqlite3.db"),
#'         organism = "E.coli",
#'         database_version = "1",
#'         genome_version = "1"
#'     )
#'
#' ## Lists all available identifiers for "araC"
#' ## Guess name
#' guess_id("araC", e_coli_regulondb)
#'
#' ## Guess id
#' guess_id("ECK120000050", e_coli_regulondb)
#'
#' ## Guess bnumber
#' guess_id("b0064", e_coli_regulondb)
#'
#' @export

guess_id <- function(gene, regulondb) {
    # Function checks
    stopifnot(validObject(regulondb))

    if (class(gene) != "character") {
        stop("'genes' should be a character vector of gene identifiers.",
            call. = FALSE)
    }
    if (length(gene) > 1) {
        stop("'genes' should be a character vector of length one", call. = FALSE)
    }

    gene_id_type <- "name"

    # Guess bnumber
    if (grepl("^b[0-9]", gene)) {
        gene_id_type <- "bnumber"
        return(gene_id_type)
    }
    # Guess id
    if (grepl("(eck|ECK)12", gene)) {
        gene_id_type <- "id"
        # Guess name
    } else if (grepl("[a-z]", gene)) {
        gene_id_type <- "name"
    }

    return(gene_id_type)
}
