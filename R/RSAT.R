#' @title Call RSAT.
#' @description This function access RSAT via SOAP.
#' @author José Alquicira Hernández & Joselyn Chávez
#' @param method Name of the method to be used from RSAT. See \code{\link{RetrieveSeq}}, \code{\link{RandomSeq}}
#' @param parameters List of parameters provided to method
#' @return an R object with results retrieved from RSAT
#' @examples
#'
#' # Get 50 random protein sequences of length 100 in fasta format. Use number 300 as seed.
#'
#' RSAT(method = "random_seq",
#'     parameters = list(sequence_length = 100,
#'                       repetition = 50,
#'                       type = "protein",
#'                       seed = 300,
#'                       format = "fasta"))
#' @export

RSAT <- function(method, parameters = NULL){

  if(!is.null(parameters)){
    request <- BuildXml(method = method,
                        parameters = parameters)
  }else{
    request <- BuildXml(method = method)
  }


  out <- tryCatch({
    # Communicates to RSAT
    res <- POST("http://embnet.ccg.unam.mx/rsat/web_services/RSATWS.cgi",
                body = request,
                content_type("text/xml; charset=utf-8"))
    stop_for_status(res)
    res <- content(res)
    # Extracts result from XML response
    res.format <- xmlToList(xmlParse(res))
    res.format$Body[[1]]$response$client$text
  },
  error = function(cond){
    #message(cond)
    message("\nError details:")
    # Extracts result from XML response
    res <- content(res)
    res.format <- xmlToList(xmlParse(res))
    return(res.format$Body)
  }
  )
  return(out)
}
