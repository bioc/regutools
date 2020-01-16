## ----echo = TRUE--------------------------------------------------------------
library('regutools')

## Other packages used
library('Biostrings')

## to be replaced with AnnotationHub call ##
download_database(tempdir(), overwrite=FALSE)

ecoli_regulondb <- 
  regulondb(
    database_path = file.path( tempdir(), "regulondb_sqlite3.db" ),
    organism="E.coli", database_version="version 1", 
    genome_version="version 1" )

ecoli_regulondb


## ----echo = TRUE--------------------------------------------------------------
list_datasets( ecoli_regulondb )

## -----------------------------------------------------------------------------
head( list_attributes( ecoli_regulondb, "GENE" ), 8 )

## ----echo = TRUE--------------------------------------------------------------

get_dataset( regulondb=ecoli_regulondb, 
            dataset="GENE",
            attributes=c("posleft", "posright", "strand", "name"), 
            filters=list("name"=c("araC","crp","lacI")) )


## ----echo = TRUE--------------------------------------------------------------
get_dataset( ecoli_regulondb, attributes=c("posright","name"),
        filters=list( "posright"=c(1, 5000) ),
        interval = "posright",
        dataset="GENE" )

## -----------------------------------------------------------------------------
res <- get_dataset( regulondb=ecoli_regulondb, 
            dataset="GENE",
            attributes=c("posleft", "posright", "strand", "name"), 
            filters=list("name"=c("araC","crp","lacI")) )
slotNames(res)

## -----------------------------------------------------------------------------
convert_to_granges(res)

## -----------------------------------------------------------------------------
get_dataset( regulondb=ecoli_regulondb, 
            dataset="GENE",
            attributes=c("posleft", "posright", "strand", "name"), 
            filters=list("name"=c("araC","crp","lacI")), 
            output_format="GRanges" )

## -----------------------------------------------------------------------------
res_dnastring <- get_dataset( regulondb=ecoli_regulondb, 
            dataset="GENE",
            attributes=c("posleft", "posright", "strand", "name", "dna_sequence"), 
            filters=list("name"=c("araC","crp","lacI")) )
res_dnastring <- convert_to_biostrings( res_dnastring, seq_type="DNA" )
res_dnastring
mcols( res_dnastring )

## -----------------------------------------------------------------------------
res_prodstring <- get_dataset( regulondb=ecoli_regulondb, 
            dataset="GENE",
            attributes=c("posleft", "posright", "strand", "name", "product_sequence"), 
            filters=list("name"=c("araC","crp","lacI")) )
res_prodstring <- convert_to_biostrings( res_prodstring, seq_type="product" )
res_prodstring
mcols( res_prodstring )

## ----echo = TRUE--------------------------------------------------------------
get_dataset( ecoli_regulondb, 
             attributes=c("posright","name"),
             filters=list("name"="ara"),
             partialmatch ="name",
             dataset="GENE" )

## ----echo = TRUE--------------------------------------------------------------

get_dataset(
  ecoli_regulondb, 
  attributes = c("name", "strand", "posright", "product_name"), 
           dataset = "GENE",
           filters = list(posright=c("2000","4000000") ),
           interv="posright" )

## -----------------------------------------------------------------------------
nrow(get_dataset(
  ecoli_regulondb, 
  attributes = c("name", "strand", "posright", "product_name"), 
           dataset = "GENE",
           filters = list(name=c("ARA"),
                          product_name=c("Ara"),
                          strand=c("forward"),
                          posright=c("2000","4000000")
           ),
           and=FALSE,
           partialmatch = c("name", "product_name") ,
           interv="posright" ))

## -----------------------------------------------------------------------------
nrow(get_dataset(
  ecoli_regulondb, 
  attributes = c("name", "strand", "posright", "product_name"), 
           dataset = "GENE",
           filters = list(name=c("ARA"),
                          product_name=c("Ara"),
                          strand=c("forward"),
                          posright=c("2000","4000000")
           ),
           and=TRUE,
           partialmatch = c("name", "product_name") ,
           interv="posright" ))

## ----echo = TRUE--------------------------------------------------------------
get_gene_regulators( ecoli_regulondb, c("araC","fis","crp") )

## ----echo = TRUE--------------------------------------------------------------
head( get_regulatory_network( ecoli_regulondb ) )

## -----------------------------------------------------------------------------
get_regulatory_summary( ecoli_regulondb, gene_regulators = c("araC", "modB") )

