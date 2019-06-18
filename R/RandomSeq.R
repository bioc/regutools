#' @title Generate random sequences
#' @description Generates random sequences with RSAT.
#' @author José Alquicira Hernández & Joselyn chávez
#' @param n Number of sequences to generate
#' @param seq.length Length of sequence to generate.
#' @param format Format of sequence(s) to generate
#' @param type Type of sequence(s) to generate (protein | DNA | other).
#' @param seed Seed for the random generator.
#' @param alphabet Alphabet. Must be followed by residue frequencies expressed precisely this way:    a:t # c:g #
#' @param expfreq Expected frequencies of oligomers in sequence(s) to generate. When this option is used, the sequences are generated according to a Markov chain.
#' @param bg.model Background model. Automatically load a pre-calibrated expected frequency file from the RSAT genome distribution. When this option is used, the options organism and \code{oligo.length} are also required, to indicate the organism and the oligonucleotide length, respectively.
#' This option is incompatible with the option \code{expfreq}.
#' Type of sequences used as background model for estimating expected oligonucleotide frequencies (supported models):
#' \itemize{
#' \item \code{equi} (equiprobable residue frequencies [default]).
#' \item \code{upstream} (all upstream sequences, allowing overlap with upstream ORFs. Requires to speciy a model organism).
#' \item \code{upstream-noorf} (all upstream sequences, preventing overlap with upstream ORFs. Requires to specify a model organism).
#' \item \code{intergenic} (intergenic frequencies. Whole set of intergenic regions, including upstream and downstream sequences. Requires to specify a model organism).
#' }
#' @param organism Name of the organism when using a background model.
#' @param oligo.length Length of oligomer when using a background model.
#' @param length.file Length file. Allows to generate random sequences with the same lengths as a set of reference sequences.The length file contains two columns : sequence ID (ignored) and sequence length.
#' @param lw Line width (0 for whole sequence on one line).
#' @return a string with results retrieved from RSAT
#' @examples
#' res <- RandomSeq(n = 3, seq.length = 100, format = "fasta",
#'                 type = "DNA", seed = 12,
#'                 organism = "Escherichia_coli_K12")
#' cat(res)
#' @export



RandomSeq <- function(n, seq.length, format = NULL, type = NULL, seed = NULL,
                      alphabet = NULL, expfreq = NULL, bg.model = NULL,
                      organism = NULL, oligo.length = NULL,
                      length.file = NULL, lw = NULL){

  !missing(n) || stop('parameter n is missing. A number of sequences to be generated is required')
  !missing(seq.length) || stop('parameter seq.length is missing. A length of sequences to be generated is required')

  if(!is.null(type)){
    (type %in% c("DNA", "protein", "other")) || stop('Only "DNA", "protein" and "other" sequence types are allowed')
  }


  parameters <- list(repetition = n,
                     sequence_length = seq.length,
                     format = format,
                     type = type,
                     seed = seed,
                     alphabet = alphabet,
                     expfreq = expfreq,
                     bg_model = bg.model,
                     organism = organism,
                     oligo_length = oligo.length,
                     length_file = length.file,
                     line_width = lw)

  res <- RSAT(method = "random_seq",
              parameters = parameters)

  if (format == "fasta") { res <- strsplit(res, split = ">")[[1]][-1] }

  if (format == "WC" || format == "IG") {
    ifelse(format == "WC", res <- strsplit(res, split = '\\\n;', fixed = T), res <- strsplit(res, split = '1\n;', fixed = T))
    res_temp <- strsplit(res[[1]][1], split = ";")[[1]][-1]
    res[[1]][1] <- paste0(res_temp, collapse = ";")}

  return(res)
}
