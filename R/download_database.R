#' Download regulondb database
#' @description If the file does not exist, this function downloads regulondb_sqlite3.db into the extdata directory of the regutools package.
#' @return local file regulondb_sqlite3.db
#' @export
#'
#' @examples download_database()
download_database <- function(path=".") {
  destfile <-  file.path(path, "regulondb_sqlite3.db")
  if (!file.exists(destfile)) {
    e <- simpleError("Download error")
    tryCatch(download.file(url = "https://www.dropbox.com/s/3xg6b2c0era21q8/regulondb_sqlite3.db?dl=1",
                           destfile ), error = function(e) e)
  } else {print("regulondb_sqlite3.db already exists in local directory")}
}
