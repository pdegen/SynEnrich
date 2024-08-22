suppressMessages(library(clusterProfiler))

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
keytype <- args[2]
organismKEGG <- args[3]
outfile <- args[4]

df <- read.csv(input_file, row.names = 1)


create_gene_converter <- function(df, OrgDb=org.Hs.eg.db) {

    if ("ENTREZID" %in% names(df)) {
        print("Entrez ID already in df...")
        write.csv(df$ENTREZID, outfile)
    }

    df[[keytype]] <- row.names(df)

    # Convert to ENTREZ ID
    # We will lose some genes here because not all IDs will be converted

    ids<-bitr(row.names(df), fromType = keytype, toType = "ENTREZID", OrgDb=OrgDb)
    write.csv(ids, outfile, row.names=FALSE)
}

if (organismKEGG == "hsa") {
    suppressMessages(library(org.Hs.eg.db))
    OrgDb <- org.Hs.eg.db
} else if (organismKEGG == "mmu") {
    suppressMessages(library(org.Mm.eg.db))
    OrgDb <- org.Mm.eg.db
} else {
    stop(paste("Organism not yet implemented:", organismKEGG))
}

create_gene_converter(df, OrgDb=OrgDb)