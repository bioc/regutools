#' @title Constructs a logical condition to query database
#' @description Given a list of filters, this function builds a logical condition to query database using intervals.
#' The output is used in [build_condition()].
#' @author Carmina Barberena Jonás, Jesús Emiliano Sotelo Fonseca, José Alquicira Hernández
#' @param interval the filters whose values will be considered as interval
#' @inheritParams existing_partial_match
#' @return A string with a sql logical condition to query the dataset
#' @examples
#' non_existing_intervals(filters = list(name="ara",strand = "for"),
#'                      interval = NULL,
#'                      operator = "AND",
#'                      partialmatch = c("name", "strand"))
#' @export

non_existing_intervals<-function(filters, interval, operator, partialmatch){
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
    condition.partialmatch<- existing_partial_match(filters, partialmatch, operator)
    if ((length(partialmatch)+length(interval)==length(filters))){
      return(condition.partialmatch)
    }else{
    condition.partialmatch<- existing_partial_match(filters, partialmatch, operator)
    condition.pmandnoin<-paste(condition.partialmatch, condition.non.interv, sep = operator, collapse = operator)
    return(condition.pmandnoin)
  }}else{
  return(condition.non.interv)
  }
}
