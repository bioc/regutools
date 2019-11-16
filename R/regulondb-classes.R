#' The regulondb class
#' @description The regulondb class is an extension of the SQLiteConnection, which as the name
#' suggests, consists of an SQLite connection to a database with the table design of the
#' RegulonDb database. In addition to the slots defined in the SQLiteConnection object,
#' the regulondb class also contains additional slots to store
#' information about database versions, organism information and genome build
#' versions.
#' @slot organism A character vector with the name of the organism of the database.
#' @slot genome_version A character vector with the version of the genome build.
#' @slot database_version A character vector with the version of regulondb build.
#' @exportClass
setClass(
  "regulondb",
  contains="SQLiteConnection",
  slots=list(
    organism="character",
    genome_version="character",
    database_version="character"
    )
)

#' Constructor function of a regulondb class
#' @aliases build_regulondb
#' @description The build_regulondb function is a constructor function of a regulondb
#' class.
#' @slot organism A character vector with the name of the organism of the database.
#' @slot genome_version A character vector with the version of the genome build.
#' @slot database_version A character vector with the version of regulondb build.
#' @export
build_regulondb <- function( database_path, organism, genome_version, database_version ){
  stopifnot( file.exists( database_path ) )
  stopifnot( is(class(organism), "character") )
  stopifnot( is(class(genome_version), "character") )
  stopifnot( is(class(database_version), "character") )
  dbconn <- dbConnect( SQLite(), database_path )
  new( "regulondb", dbconn, organism=organism,
       genome_version=genome_version,
       database_version=database_version)
}

