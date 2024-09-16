suppressMessages(library(clusterProfiler))

create_gene_converter <- function(input, outfile, keytype, OrgDb=org.Hs.eg.db) {

    if (is.data.frame(input)) {

        if ("ENTREZID" %in% names(input)) {
            print("Entrez ID already in df...")
            write.csv(input$ENTREZID, outfile)
        }

        input[[keytype]] <- row.names(input)
        ids <- bitr(row.names(input), fromType = keytype, toType = c("ENTREZID", "SYMBOL"), OrgDb = OrgDb)

    } else if (is.vector(input) || is.list(input)) {
        ids <- bitr(input, fromType = keytype, toType = c("ENTREZID", "SYMBOL"), OrgDb = OrgDb)
    } else {
        stop("Input must be either a data frame or a list/vector.")
    }
    write.csv(ids, outfile, row.names=FALSE)
}

if (!interactive()) {

    args <- commandArgs(trailingOnly = TRUE)
    input_file <- args[1]
    keytype <- args[2]
    organismKEGG <- args[3]
    outfile <- args[4]
    df <- read.csv(input_file, row.names = 1)

    if (organismKEGG == "hsa") {
        suppressMessages(library(org.Hs.eg.db))
        OrgDb <- org.Hs.eg.db
    } else if (organismKEGG == "mmu") {
        suppressMessages(library(org.Mm.eg.db))
        OrgDb <- org.Mm.eg.db
    } else {
        stop(paste("Organism not yet implemented:", organismKEGG))
    }

    create_gene_converter(df, outfile, keytype, OrgDb=OrgDb)
}
