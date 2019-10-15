#' Download regulondb database from online file
#'
#' @return local file regulondb_sqlite3.db
#' @export
#'
#' @examples download_database()
download_database <- function() {
  destfile <-  paste0(system.file(package = "regutools"),"/extdata/regulondb_sqlite3.db")
  if (!file.exists(destfile)) {
    e <- simpleError("Download error")
    tryCatch(download.file(url = "https://www.dropbox.com/s/3xg6b2c0era21q8/regulondb_sqlite3.db?dl=0",
                           destfile, method = "wget"), error = function(e) e)
  } else {print("regulondb_sqlite3.db already exists in local directory")}
}
