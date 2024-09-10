# activate clusterprofiler conda env first (.snakemake/conda)

library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(dplyr)
library(GO.db) 

get_go_gmt <- function(org.db, overwrite=FALSE) {
# Extract GO terms from org.db

    gmt_file <- file.path("resources", paste0("go_terms.",deparse(substitute(org.db)),".gmt"))

    if (file.exists(gmt_file) && !overwrite) {
        print("Existing gmt file not overwritten")
    }

    go_terms <- AnnotationDbi::select(org.db, 
                                    keys = keys(org.db, keytype = "GO"), 
                                    columns = c("GO", "SYMBOL"), 
                                    keytype = "GO")

    # GSEApy complains if not all upper because of Enrichr
    go_terms$SYMBOL <- toupper(go_terms$SYMBOL)

    # Add GO term descriptions
    go_terms$DESCRIPTION <- Term(go_terms$GO)


    go_gmt <- go_terms %>%
    group_by(GO, DESCRIPTION) %>%
    summarise(Genes = paste(SYMBOL, collapse = "\t"))

    # Combine GO ID and description for the first column
    go_gmt <- go_gmt %>%
    mutate(Term = paste(GO, DESCRIPTION, sep = "\t"))

    go_fmt <- go_gmt[, c("Term", "Genes")]

    # Write to GMT format
    write.table(go_gmt, gmt_file, sep="\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
}

get_go_gmt(org.Hs.eg.db)
get_go_gmt(org.Mm.eg.db)