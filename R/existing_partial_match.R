#' @title Constructs a logical condition to query database
#' @description Given a list of filters, this function builds a logical
#' condition to query database using intervals.
#' The output is used in [existing_intervals()] and [non_existing_intervals()].
#' @author Carmina Barberena Jonás, Jesús Emiliano Sotelo Fonseca,
#' José Alquicira Hernández
#' @param filters List of filters to be used. The names should correspond to
#' the attribute and the values correspond to the condition for selection.
#' @param partialmatch name of the condition(s) with a string pattern for full
#' or partial match in the query.
#' @param operator A string indicating if all the filters (AND) or some of them
#' (OR) should be met.
#' @return A `character(1)` with the sql logical condition to query the dataset.
#' @examples
#'
#' ## Build the SQL query for existing partial matches for ara
#' existing_partial_match(
#'     filters = list(
#'         name = c("ara"),
#'         strand = c("forward"),
#'         posright = c("2000", "40000")
#'     ),
#'     partialmatch = "name",
#'     operator = "AND"
#' )
#'
#' @export
existing_partial_match <-
    function(filters, partialmatch, operator) {
        existing.partialmatch.index <- names(filters) %in% partialmatch
        existing.partialmatch <- filters[existing.partialmatch.index]
        condition.format.partialmatch <-
            mapply(paste0, filters[existing.partialmatch.index], "%'", SIMPLIFY = FALSE)
        condition.format.partialmatch <-
            mapply(
                paste,
                names(condition.format.partialmatch),
                condition.format.partialmatch,
                sep = " like  '%" ,
                SIMPLIFY = FALSE
            )
        condition.format.partialmatch <-
            lapply(condition.format.partialmatch, function(x) {
                paste0("(", paste(x, collapse = " AND "), ")")
            })
        condition.partialmatch <-
            paste(unlist(condition.format.partialmatch),
                collapse = paste0(" ", operator, " "))
        return(condition.partialmatch)
    }
