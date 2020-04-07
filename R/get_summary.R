#' @title Return summary of gene regulation.
#' @description This function takes the output of [get_gene_regulators()] with
#' format multirow,
#' onerow or table, or a vector with genes and retrieves information about the
#' TFs and their regulated genes
#' @param regulondb A [regulondb()] object.
#' @param gene_regulators Result from [get_gene_regulators()] or vector of genes
#' @return A data frame with the following columns:
#' \itemize{
#' \item The name or gene of TF
#' \item Regulated Genes per TF
#' \item Percent of regulated genes per TF
#' \item positive, negative or dual regulation
#' \item Name(s) of regulated genes
#' }
#' @keywords regulation retrieval, summary, networks,
#' @author Carmina Barberena Jonas, Jesús Emiliano Sotelo Fonseca,
#' José Alquicira Hernández, Joselyn Chávez
#' @examples
#' ## Connect to the RegulonDB database if necessary
#' if(!exists('regulondb_conn')) regulondb_conn <- connect_database()
#'
#' ## Build the regulon db object
#' e_coli_regulondb <-
#'     regulondb(
#'         database_conn = regulondb_conn,
#'         organism = "E.coli",
#'         database_version = "1",
#'         genome_version = "1"
#'     )
#'
#' ## Get the araC regulators
#' araC_regulation <-
#'     get_gene_regulators(
#'         e_coli_regulondb,
#'         genes = c("araC"),
#'         format = "multirow",
#'         output.type = "TF"
#'     )
#'
#' ## Summarize the araC regulation
#' get_regulatory_summary(e_coli_regulondb, araC_regulation)
#'
#' ## Retrieve summary of genes 'araC' and 'modB'
#' get_regulatory_summary(e_coli_regulondb,
#'     gene_regulators = c("araC", "modB"))
#'
#' ## Obtain the summary for 'ECK120000050' and 'modB'
#' get_regulatory_summary(e_coli_regulondb,
#'     gene_regulators = c("ECK120000050", "modB"))
#'
#' @export
#' @importFrom stats complete.cases

get_regulatory_summary <- function(regulondb, gene_regulators) {
    regulation <- gene_regulators
    # Check that class is regulation or character

    #Regulation of a gene list
    if (is(regulation,  "character")) {
        regulation <- get_gene_regulators(regulondb, genes = regulation)
    } else if (!is(regulation, "regulondb_result")) {
        stop(
            "The parameter gene_regulators must be a regulondb_result object
            resulting from a call to get_gene_regulators or a vector of genes.",
            call. = FALSE
        )
    }

    # Conver onerow to multirow
    if (metadata(regulation)$format == "onerow") {
        regulation <-
            get_gene_regulators(regulondb,
                                genes = as.character(regulation$genes))
    }

    # Convert table to multirow
    if (metadata(regulation)$format == "table") {
        regulation <- get_gene_regulators(regulondb,
                                    genes = as.character(rownames(regulation)))
    }

    TF_counts <-
        data.frame(table(regulation$regulators), stringsAsFactors = FALSE)

    summary <- apply(TF_counts, 1, function(x) {
        #Get rows for a specific TF
        TF_data <- regulation[regulation[, "regulators"] == x[1], ]

        #Count regulated effects
        effect <- table(TF_data$effect)
        #Adds regulation +, - or +/- in case their frequency is 0.
        if (!"+" %in% names(effect)) {
            effect <- c(effect, "+" = 0)
        }
        if (!"-" %in% names(effect)) {
            effect <- c(effect, "-" = 0)
        }
        if (!"+/-" %in% names(effect)) {
            effect <- c(effect, "+/-" = 0)
        }

        #Concatenates row to include the summary information for each TF
        summary_row <- c(
            TF = x[1],
            #TF name
            regulated_number = x[2],
            #Number of genes regulated per TF
            regulated_percentage =
                (as.numeric(x[2]) / length(regulation$genes)) * 100,
            #Percent of genes in query regulated per TF
            activator = effect["+"],
            #Frequency of activation interactions
            repressor = effect["-"],
            #Frequency of repression interactions
            dual = effect["+/-"] ,
            #Frequency of dual interactions
            regulated = paste(TF_data$genes,
                            collapse = ", ")#List of genes regulated per TF
        )

        return(summary_row)

    })

    summary <- data.frame(t(summary)) #Format summary as a data.frame
    colnames(summary) <-
        c(
            "TF",
            "Regulated_genes_per_TF",
            "Percent",
            "Activator",
            "Repressor",
            "Dual",
            "Regulated_genes"
        )

    summary <- dataframe_to_dbresult(summary, regulondb, "NETWORK")
    return(summary)
}
