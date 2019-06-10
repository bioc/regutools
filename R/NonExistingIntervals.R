#' @title Constructs a logical condition to query database
#' @description Given a list of filters, this function builds a logical condition to query database using intervals.
#' The output is used in \code{\link{BuildCondition}}.
#' @author Carmina Barberena Jonás, Jesús Emiliano Sotelo Fonseca, José Alquicira Hernández
#' @param filters List of filters to be used. The names should correspond to the attribute and the values correspond to the condition for selection.
#' @param interval the filters whose values will be considered as interval
#' @param operator A string indicating if all the filters (AND) or some of them (OR) should be met
#' @return A string with a sql logical condition to query the dataset
#' @examples
#' NonExistingIntervals(filters = list(name="ara",strand = "for"),
#'                      interval = NULL,
#'                      operator = "AND",
#'                      partialmatch = c("name", "strand"))
#' @export

NonExistingIntervals<-function(filters, interval, operator, partialmatch){
  if (!(length(partialmatch)+length(interval)==length(filters))){
  non.interv.index<-!((names(filters) %in% interval)|(names(filters) %in% partialmatch))
  non.interv <-filters[non.interv.index]
  condition.format.non.interv <- mapply(paste0, filters[non.interv.index], "'", SIMPLIFY = FALSE)
  condition.format.non.interv <- mapply(paste, names(condition.format.non.interv), condition.format.non.interv, sep = " = '" , SIMPLIFY = FALSE)
  condition.format.non.interv <- lapply(condition.format.non.interv, function(x){
    paste0("(", paste(x, collapse = " OR "), ")") })
  condition.non.interv <- paste(unlist(condition.format.non.interv), collapse = paste0(" ", operator, " "))
  }
  if (!is.null(partialmatch)){
    condition.partialmatch<-ExistingPartialMatch(filters, partialmatch, operator)
    if ((length(partialmatch)+length(interval)==length(filters))){
      return(condition.partialmatch)
    }else{
    condition.partialmatch<-ExistingPartialMatch(filters, partialmatch, operator)
    condition.pmandnoin<-paste(condition.partialmatch, condition.non.interv, sep = operator, collapse = operator)
    return(condition.pmandnoin)
  }}else{
  return(condition.non.interv)
  }
}
