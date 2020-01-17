#' Download regulondb database
#' @description If the file does not exist, this function downloads
#' regulondb_sqlite3.db into the extdata directory of the regutools package.
#' @param path The directory where the database will be downloaded to.
#' @param overwrite `logical(1)` specifying whether to overwrite the database
#' file.
#' @return local file regulondb_sqlite3.db
#' @export
#' @importFrom utils download.file
#'
#' @examples
#'
#' ## Download the database if necessary
#' if(!file.exists(file.path(tempdir(), 'regulondb-sqlite3.db'))) {
#'     download_database(tempdir())
#' }
#'
download_database <- function(path = tempdir(), overwrite = FALSE) {
    destfile <-  file.path(path, "regulondb_sqlite3.db")
    url <-
        "https://www.dropbox.com/s/eod8vdq4fthvjcr/regulondb_v10.6.2_DM_sqlite3.db?dl=1"
    e <- simpleError("Download error")
    if (overwrite) {
        tryCatch(
            download.file(url = url, destfile),
            error = function(e)
                e
        )
    } else {
        if (!file.exists(destfile)) {
            tryCatch(
                download.file(url = url, destfile),
                error = function(e)
                    e
            )
        } else {
            print(
                "regulondb_sqlite3.db already exists in local directory, set overwrite = TRUE if you want to replace existing file."
            )
        }
    }
}
