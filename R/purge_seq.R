#' @title Mask repeated fragments of an input sequence
#' @description Returns your sequence with substitution of repeated fragments by 'n' .
#' @author Joselyn Chávez
#' @param sequence Sequence to purge
#' @param format Sequence format. Supported: 'IG' (Intelligenetics), 'WC' (wconsensus), 'raw', 'fasta'. Defaul is 'fasta'.
#' @param match.length Minimal match length. Default is 40.
#' @param mismatch Number of mismatches allowed. Default is 3.
#' @param both.strand Discard duplications on the direct strand only (FALSE) or on the reverse complement as well (TRUE). Default is TRUE.
#' @param delete Delete repeats instead of masking them if TRUE.
#' @param mask.short Mask (replace by N characters) sequences shorter than the specified length.
#' @return An string with results retrieved from RSAT
#' @examples
#' # Purge sequence with repeats of length 10
#' res <- purge_seq(sequence = ">seq1\nTCGCAATCTACTTACTATTGCTTACTATTGCTTACTATTGCGAGCGCAGA")

#' res <- purge_seq(sequence = "TCGCAATCTACTTACTATTGCTTACTATTGCTTACTATTGCGAGCGCAGA",
#'                    match.length = 10,
#'                    format = "raw")
#' cat(res)
#' @export

purge_seq <- function(sequence, format = "fasta", match.length = NULL, mismatch = NULL,
                        both.strand = TRUE, delete = NULL, mask.short = NULL){

  if(sequence == '' | is.na(sequence) | is.null(sequence)){
    stop("A sequence is required to purge")
  }

  delete <- ifelse(delete, 1, 0)
  both.strand <- ifelse(both.strand, 2,1)


  parameters <- list(sequence = sequence,
                     format = format,
                     match_length = match.length,
                     mismatch = mismatch,
                     str = both.strand,
                     delete = delete,
                     mask_short = mask.short)

  res <- RSAT(method = 'purge_seq', parameters = parameters)

  #if (format == "fasta") { res <- strsplit(res, split = ">")[[1]][-1] }

  if (format == "WC" || format == "IG") {
    ifelse(format == "WC", res <- strsplit(res, split = '\\\n;', fixed = T), res <- strsplit(res, split = '1\n;', fixed = T))
    res_temp <- strsplit(res[[1]][1], split = ";")[[1]][-1]
    res[[1]][1] <- paste0(res_temp, collapse = ";") }

  return(res)
}
