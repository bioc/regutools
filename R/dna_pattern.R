#' @title Searches all occurrences of a pattern within DNA sequences.
#' @description Returns position of pattern match in the input sequence .
#' @author Joselyn Chávez
#' @param sequence Input sequence
#' @param format Input sequence format. Supported: 'IG' (Intelligenetics), 'WC' (wconsensus), 'raw', 'fasta'. Defaul is 'fasta'.
#' @param substitution Number of substitutions allowed.
#' @param pattern Pattern to match.
#' @param id Pattern identifier.
#' @param origin Origin for the calculation of positions (0 for end of sequence).
#' @param noov No overlapping of oligos allowed if TRUE. Disable the detection of overlapping matches for self-overlapping patterns (e.g. TATATA, GATAGA).
#' @param both.strand Oligonucleotide occurrences found on both stands are summed (TRUE) or not (FALSE). Default is TRUE
#' @param sort.oligo sort oligomers according to overrepresentation if TRUE.
#' @param threshold Threshold on match count.
#' @return A data frame with columns
#' \itemize{
#' \item {patternID}
#' \item{strand}
#' \item{pattern}
#' \item{seqID}
#' \item{start}
#' \item{end}
#' \item{matching_seq}
#' \item{score}
#' }
#' @examples
#' # Find pattern in both strands of input sequence
#' res <- dna_pattern(sequence = "TCGCAATCTACTTACTATTGCTTACTATTGCTTACTATTGCGAGCGCAGA",
#'                   format = "raw",
#'                   pattern = "TCGCAAT")
#' res <- dna_pattern(sequence = ">seq1\nTCGCAATCTACTTACTATTGCTTACTATTGCTTACTATTGCGAGCGCAGA",
#'                   pattern = "TCGCAAT",
#'                   both.strand = F,
#'                   id = "patt_1")
#' cat(res)
#' @export

dna_pattern <- function(sequence, format = "fasta", substitution = NULL, pattern,
                     id = NULL, origin = NULL, noov = NULL, both.strand = TRUE,
                     sort.oligo = NULL, threshold = NULL){

  if(sequence == '' | is.na(sequence) | is.null(sequence)){
    stop("A sequence is required")
  }
  if(pattern == '' | is.na(pattern) | is.null(pattern)){
    stop("A pattern is required")
  }
  noov <- ifelse(noov, 1, 0)
  both.strand <- ifelse(both.strand, 2,1)
  sort.oligo <- ifelse(sort.oligo, 1, 0)

  parameters <- list(sequence = sequence,
                     format = format,
                     subst = substitution,
                     pattern = pattern,
                     id = id,
                     origin  = origin,
                     noov = noov,
                     str = both.strand,
                     sort = sort.oligo,
                     th = threshold)

  res <- RSAT(method = 'dna_pattern', parameters = parameters)

  res <- as.data.frame(strsplit(res, split = "\n"), col.names = "V1")
  res <- as.data.frame(t(as.data.frame(strsplit(as.character(res$V1), split = "\t"))))
  colnames(res) <- c("patternID", "strand", "pattern", "seqID", "start", "end", "matching_seq", "score")
  rownames(res) <- NULL
  for (i in 1:8) {names(res[,i]) <- NULL}

  return(res)
}
