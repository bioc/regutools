#' @title Constructs a particular logical condition to query database
#' @description Given a list of filters, this function builds a logical
#' condition to query database using intervals.
#' The output is used in [build_condition()].
#' @author Carmina Barberena Jonás, Jesús Emiliano Sotelo Fonseca,
#' José Alquicira Hernández, Joselyn Chávez
#' @param filters List of filters to be used. The names should correspond to
#' the attribute and the values correspond to the condition for selection.
#' @param interval the filters with values considered as interval.
#' @param partialmatch name of the condition(s) with a string pattern for full
#' or partial match in the query.
#' @param operator A string indicading if all the filters (AND) or some of
#' them (OR) should be met.
#' @return A `character(1)` with the sql logical condition to query the dataset.
#' @examples
#'
#' ## Build the SQL query for existing interval partial matches for ara
#' existing_intervals(
#'     filters = list(
#'         name = "ara",
#'         strand = "for",
#'         posright = c("2000", "40000")
#'     ),
#'     interval = c("posright"),
#'     operator = "AND",
#'     partialmatch = c("name", "strand")
#' )
#' @export


existing_intervals <-
    function(filters, interval, operator, partialmatch) {
        existing.interv.index <- names(filters) %in% interval
        existing.interv <- filters[existing.interv.index]

        # Check that interval's value is a pair
        # If they are more than two, then drop the remaining
        existing.interv <- lapply(existing.interv, function(x) {
            if (length(x) > 2) {
                warning("Only the first two values of interval will be considered.",
                    call. = FALSE
                )
                x[seq_len(2L)]
            } else if (length(x) == 1) {
                stop("Two values in the interval filter are required. ",
                    call. = FALSE
                )
            } else {
                x
            }
        })
        filters[existing.interv.index] <- existing.interv

        condition.format.interv <-
            mapply(paste0, filters[existing.interv.index], "", SIMPLIFY = FALSE)
        condition.format.interv <-
            lapply(condition.format.interv, function(y) {
                paste(c(">=", "<="), y)
            })
        condition.format.interv <-
            mapply(
                paste,
                names(condition.format.interv),
                condition.format.interv,
                sep = " ",
                SIMPLIFY = FALSE
            )
        condition.format.interv <-
            lapply(condition.format.interv, function(x) {
                paste0("(", paste(x, collapse = " AND "), ")")
            })
        condition.interv <-
            paste(unlist(condition.format.interv),
                collapse = paste0(" ", operator, " ")
            )
        if (!is.null(partialmatch)) {
            condition.partialmatch <-
                existing_partial_match(filters, partialmatch, operator)
            condition.pmandin <-
                paste(
                    condition.partialmatch,
                    condition.interv,
                    sep = operator,
                    collapse = operator
                )
            return(condition.pmandin)
        }
        return(condition.interv)
    }
