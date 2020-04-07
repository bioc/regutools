#' Connect to the regulondb database
#' @description If the file does not exist, this function downloads
#' the RegulonDB SQLite database file into the `path` prior to loading it.
#' @param path The directory where the database will be downloaded to.
#' @param overwrite `logical(1)` specifying whether to overwrite the database
#' file.
#' @param ah An `AnnotationtHub` object
#' [AnnotationHub-class][AnnotationHub::AnnotationHub-class]. Can be `NULL`
#' if you want to force to use the backup download mechanism.
#' @return An [SQLiteConnection-class][RSQLite::SQLiteConnection-class]
#' connection to the RegulonDB database.
#' @export
#' @importFrom utils download.file
#' @import AnnotationHub
#' @importFrom AnnotationDbi dbFileConnect
#'
#' @examples
#'
#' ## Connect to the RegulonDB database if necessary
#' if(!exists('regulondb_conn')) regulondb_conn <- connect_database()
#'
#' ## Connect to the database without using AnnotationHub
#' regulondb_conn_noAH <- connect_database(ah = NULL, overwrite = TRUE)
#'
connect_database <-
    function(path = tempdir(),
            overwrite = FALSE,
            ah = AnnotationHub::AnnotationHub()) {
        if (!is.null(ah)) {
            ## Check input
            stopifnot(methods::is(ah, 'AnnotationHub'))

            ## Query AH
            q <-
                AnnotationHub::query(
                    ah,
                    pattern = c(
                        'RegulonDB SQLite database version v10.6.2_DM for the regutools Bioconductor package',
                        'regutools;RegulonDB;v10.6.2_DM'
                    )
                )
            if (length(q) == 1) {
                ## Return the connection
                return(q[[1]])
            }
        }


        ## Otherwise, use the Dropbox version
        destfile <-
            file.path(path, "regulondb_sqlite3_v10.6.2_DM.db")
        url <-
            "https://www.dropbox.com/s/eod8vdq4fthvjcr/regulondb_v10.6.2_DM_sqlite3.db?dl=1"
        e <- simpleError("Download error")
        if (overwrite) {
            tryCatch(
                download.file(url = url, destfile, mode = 'wb'),
                error = function(e)
                    e
            )
        } else {
            if (!file.exists(destfile)) {
                tryCatch(
                    download.file(url = url, destfile, mode = 'wb'),
                    error = function(e)
                        e
                )
            } else {
                message(
                    paste(
                        "regulondb_sqlite3_v10.6.2_DM.db already exists in",
                        "local directory, setoverwrite = TRUE if you want to",
                        "replace the existing file."
                    )
                )
            }
        }
        AnnotationDbi::dbFileConnect(destfile)
    }
