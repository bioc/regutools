#' Plot elements within a genome region
#'
#' @param regulondb A regulondb object.
#' @param from left position.
#' @param plot_map By default TRUE, if FALSE returns a GRanges object with the elements found.
#' @param to right position
#'
#' @author Joselyn Chavez
#' @return A plot with genomic elements found within a genome region, including genes and regulators.
#' @importFrom Gviz GeneRegionTrack plotTracks GenomeAxisTrack
#' @importFrom stringr str_split_fixed
#' @examples
#' ## Download the database if necessary
#' if(!file.exists(file.path(tempdir(), 'regulondb_sqlite3.db'))) {
#'     download_database(tempdir())
#' }
#'
#' ## Build the regulondb object
#' e_coli_regulondb <-
#'     regulondb(
#'         database_path = file.path(tempdir(), "regulondb_sqlite3.db"),
#'         organism = "chr",
#'         database_version = "1",
#'         genome_version = "1"
#'     )
#' plot_GRanges(e_coli_regulondb, from = 5000, to = 10000)
plot_GRanges <- function(regulondb, from, to, plot_map = TRUE) {
    # validate ranges
        if (!is.numeric(from) || !is.numeric(to)) {
        stop("Parameter 'from' and 'to' must be a number", call. = FALSE)
    }

    # search for dna_objects ("-10 promoter box", -35 promoter box", "gene", "promoter", "Regulatory Interaction","sRNA interaction","terminator")
    dna_objects <- get_dataset(regulondb,
                               dataset = "DNA_OBJECTS",
                               filters = list(posright = c(from, to)),
                               interval = "posright",
                               output_format = "GRanges")
    # search for TF binding sites

    # tf <- get_dataset(regulondb, dataset = "TF", attributes = c("name")) %>% sort() %>% unique()
    # tf_result <- get_binding_sites(e_coli_regulondb, transcription_factor = tf$name[1])
    #
    # for (i in 2:3) {
    #     partial_result <- get_binding_sites(e_coli_regulondb, transcription_factor = tf$name[i])
    #     tf_result <- c(tf_result, partial_result)
    # }

    # construct track

    dna_annotation <- AnnotationTrack(dna_objects,
                                      genome = "eschColi_K12",
                                      chromosome = "chr",
                                      name = "DNA_objects",
                                      options(ucscChromosomeNames=TRUE),
                                      transcriptAnnotation= "name",
                                      background.panel="#FFFEDB",
                                      background.title="brown")
    # optional: return GRanges result
    if (plot_map == FALSE) {
        return(dna_objects)
        stop()
    }

    # plot
    plotTracks(list(GenomeAxisTrack(),
                    dna_annotation)
               )
}

