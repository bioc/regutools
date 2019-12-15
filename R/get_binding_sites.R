#' @title Get the binding sites for a Transcription Factor (TF)
#' @description Retrieve the binding sites and genome location for a given transcription factor.
#' @author José Alquicira Hernández, Jacques van Helden, Joselyn Chávez
#' @param regulondb A regulondb object.
#' @param TF name of the transcription factor.
#' @param seq.format Default: \code{table}. Supported: \code{fasta}, \code{wconsensus}.
#' @return A dataframe with the following columns corresponding to the TFBSs associated with the TF.
#' If seq.format = "fasta", returns a string with the TFBS sequences with breakline character.
#' If the transcription factor does not exist, returns NA.
#' \itemize{
#' \item ID (identifier)
#' \item left position
#' \item right position
#' \item strand
#' \item sequence
#' }
#' @examples
#' # Get the binding sites for AraC
#'
#' get_binding_sites(e_coli_regulondb, TF = "AraC")
#'
#' # Get the fasta binding sites sequences for AraC
#'
#' fasta.seq <- get_binding_sites(e_coli_regulondb,
#'                                TF = "AraC", seq.format = "fasta")
#' cat(fasta.seq)
#' @export

get_binding_sites <- function(regulondb, TF, seq.format = "table") {
    stopifnot(validObject(regulondb))
    tfbs.raw <- tryCatch({
      get_dataset(regulondb = regulondb,
              dataset = "TF",
              attributes = c("name", "tfbs_unique"),
              filters = list(name = TF))
    },error = function(cond) return(NULL))

  if(is.null(tfbs.raw)){
    return(NA)
  }else{
    # convert raw info into a table with 1 row per TFBS and 1 col per attribute
    tfbs.table <-  as.data.frame(strsplit(x = tfbs.raw$tfbs_unique, split = " "), col.names = "V1")
    tfbs.table <- t(as.data.frame(strsplit(x = as.character(tfbs.table$V1), split =  "\t")))

    tfbs.table <- as.data.frame(tfbs.table)
    colnames(tfbs.table) <- c("ID", "left", "right", "strand", "sequence")
    rownames(tfbs.table) <- NULL
    for (i in 1:5) {names(tfbs.table[,i]) <- NULL}

    if (seq.format == "table") {return(tfbs.table)}
    if (seq.format == "fasta") {
      header <- paste0(">", tfbs.table$ID, "_", tfbs.table$left, ":", tfbs.table$right, ":", tfbs.table$strand)
      seq.string <- paste0(header,"\n",as.vector(tfbs$sequence), "\n")
      return(seq.string)
    } else { stop("seq.format must be 'table' or 'fasta' ", call. = FALSE)}

  }
}
