#' @title Retrieve gene synonyms
#' @description Given a list of genes (id, name, bnumber or gi),
#' get the gene synonyms (name, bnumber of gi).
#' @param regulondb A [regulondb()] object.
#' @param genes Character vector of gene identifiers (id, name, bnumber or gi).
#' @param from A `character()` specifying one of: id, name, bnumber of gi
#' @param to A `character()` specifying one or more of: id, name, bnumber of gi
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
#' get_gene_synonyms(e_coli_regulondb, "araC", from = "name")
#'
#' ## Retrieve only the ID
#' get_gene_synonyms(e_coli_regulondb, "araC", from = "name", to = "id")
#'
#' ## Use an ID to retrieve the synonyms
#' get_gene_synonyms(e_coli_regulondb, "ECK120000998", from = "id")
#'
#' @export

get_gene_synonyms <-
    function(regulondb,
        genes,
        from = "name",
        to = c("id", "name", "bnumber", "gi")) {
        # Function checks
        stopifnot(validObject(regulondb))

        if (length(from) > 1) {
            stop("'from' should be only one of name, bnumber or GI.", call. = FALSE)
        }

        if (class(genes) != "character") {
            stop("'genes' should be a character vector of gene identifiers.",
                call. = FALSE)
        }

        if (!all(to %in% c("id", "name", "bnumber", "gi"))) {
            stop("'to' should be a character vector with one or more of name, bnumber or GI.",
                call. = FALSE)
        }

        if (!from %in% c("id", "name", "bnumber", "gi")) {
            stop("'from' should be one or more of name, bnumber or GI.", call. = FALSE)
        }

        gene_filter <- list(genes)
        names(gene_filter) <- from

        #Convert GIs to gene names
        gene_synonyms <-
            get_dataset(
                regulondb,
                dataset = "GENE",
                attributes = to,
                filters = gene_filter
            )

        return(gene_synonyms)
    }
