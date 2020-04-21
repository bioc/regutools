## For this script
dir.create(here::here('inst', 'scripts'), recursive = TRUE)

## Download the DB file

dir.create(here::here('data-raw'))
usethis::use_data_raw()
unlink(here::here('data-raw', 'DATASET.R'))
## This function was dropped during the development of regutools
# regutools::download_database(here::here('data-raw'), overwrite = FALSE, ah = NULL)

## Rename if necessary
# file.rename(file.path(here::here('data-raw', 'regulon_sqlite3.db')),
#     file.path(here::here(
#         'data-raw', 'regulon_sqlite3_v10.6.2_DM.db'
#     )))
