context("get_regulators")
test_that("Functions to retrieve gene regulation work as expected", {
    fileToDb <- file.path(tempdir(), "regulondb_sqlite3.db")
    if (!file.exists(fileToDb))
        download_database(tempdir())
    regdb <-
        regulondb(
            fileToDb,
            organism = "prueba",
            genome_version = "prueba",
            database_version = "prueba"
        )
    rs1 <- get_gene_regulators(regdb, c("araC", "fis", "crp"))
    expect_s4_class(rs1, "regulondb_result")
    expect_true(metadata(rs1)$format == "multirow")
    rs2 <-
        get_gene_regulators(regdb, c("araC", "fis", "crp"), format = "onerow")
    expect_s4_class(rs2, "regulondb_result")
    expect_true(metadata(rs2)$format == "onerow")
    rs3 <-
        get_gene_regulators(regdb, c("araC", "fis", "crp"), format = "table")
    expect_s4_class(rs3, "regulondb_result")
    expect_true(metadata(rs3)$format == "table")
    expect_s4_class(get_regulatory_network(regdb), "regulondb_result")
    expect_s4_class(get_regulatory_summary(regdb, gene_regulators = c("araC", "modB")),
        "regulondb_result")
})
