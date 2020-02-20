#' Plot elements within a genome region
#'
#' @param regulondb A regulondb object.
#' @param from left position.
#' @param plot_map By default TRUE, if FALSE returns a GRanges object with the elements found.
#' @param to right position
#' @param genome A valid UCSC genome name.
#'
#' @author Joselyn Chavez
#' @return A plot with genomic elements found within a genome region, including genes and regulators.
#' @importFrom Gviz GeneRegionTrack plotTracks GenomeAxisTrack AnnotationTrack
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
#' @export
plot_GRanges <- function(regulondb, genome = "eschColi_K12", from, to, plot_map = TRUE) {
    # validate ranges
        if (!is.numeric(from) || !is.numeric(to)) {
        stop("Parameter 'from' and 'to' must be a number", call. = FALSE)
    }

    # search for dna_objects ("-10 promoter box", -35 promoter box", "gene", "promoter", "Regulatory Interaction","sRNA interaction","terminator")
    dna_objects <- regutools::get_dataset(regulondb,
                               dataset = "DNA_OBJECTS",
                               filters = list(posright = c(from, to)),
                               interval = "posright",
                               output_format = "GRanges")

    # optional: return GRanges result
    if (plot_map == FALSE) {
        return(dna_objects)
        stop()
    }

    # construct tracks
    dna_objects_type <- dna_objects$type %>% sort %>% unique
    if (grep("gene", dna_objects_type) > 0) {
        dna_objects_type <- c("gene", dna_objects_type[!dna_objects_type == "gene"])}

    list_dna_annotation <- list()
    for (i in dna_objects_type) {
        list_dna_annotation[i] <- Gviz::AnnotationTrack(dna_objects[dna_objects$type == i,],
                                                genome,
                                                chromosome = "chr",
                                                name = i,
                                                options(ucscChromosomeNames=TRUE),
                                                transcriptAnnotation= "name",
                                                background.panel="#FFFEDB",
                                                background.title="brown",
                                                fontsize=10) }

    # plot
    Gviz::plotTracks(c(Gviz::GenomeAxisTrack(),
                          list_dna_annotation) )
}

