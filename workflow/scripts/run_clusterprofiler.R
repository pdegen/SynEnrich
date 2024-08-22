suppressMessages(library(clusterProfiler))

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
keytype <- args[2]
organismKEGG <- args[3]
metric <- args[4]
outfile_go <- args[5]
outfile_kegg <- args[6]

print(paste("Args:",args))
print(paste("Reading clsuterProfiler input:", input_file))

#metric <- strsplit(basename(outfile_go), "syn.clusterProfiler\\.")[[1]][2]
#metric <- strsplit(metric, "\\.")[[1]][1]
print(paste("Metric:", metric))

df <- read.csv(input_file, row.names = 1)

# Check if the metric is in the columns
if (!(metric %in% colnames(df))) {
    if (metric == "neg_signed_logpval") {
        message(paste("Adding", metric, "to df"))
        df$neg_signed_logpval <- -sign(df$logFC) * log10(df$PValue)
    } else {
        stop(paste("Metric", metric, "not in columns!"))
    }
}

run_clusterProfiler <- function(df, outfile_go, outfile_kegg,
                                metric,
                                organism.KEGG,
                                organism.GO, 
                                overwrite=FALSE, seed=1234) 
{
  set.seed(seed)

  print(outfile_go)
  print(outfile_kegg)

  if (file.exists(outfile_go) && file.exists(outfile_kegg) && !overwrite) {
    print("Existing files not overwritte, skipping")
    return
  }

  start_time <- Sys.time()

  geneList <- df[[metric]]
  names(geneList) <- df$ENTREZID
  geneList = sort(geneList, decreasing = TRUE)

  if ((length(outfile_go) > 0 && !file.exists(outfile_go)) || overwrite) {

    print("Running GO...")
    ego3 <- gseGO(geneList     = geneList,
                  OrgDb        = organism.GO,
                  ont          = "ALL", ## CC MF BP
                  minGSSize    = 10,
                  maxGSSize    = 500,
                  pvalueCutoff = 1,
                  eps = 0,
                  seed = TRUE,
                  verbose = FALSE)
    write.csv(ego3,outfile_go)
    print(paste("Wrote GO to:", outfile_go))
  }

  if ((length(outfile_kegg) > 0 && !file.exists(outfile_kegg)) || overwrite)  {

    print("Running KEGG...")
    kegg <- gseKEGG(geneList     = geneList,
                  organism        = organism.KEGG,
                  minGSSize    = 10,
                  maxGSSize    = 500,
                  pvalueCutoff = 1,
                  eps = 0,
                  seed = TRUE,
                  verbose = FALSE)
    write.csv(kegg,outfile_kegg)
    print(paste("Wrote KEGG to:", outfile_kegg))
  }

  end_time <- Sys.time()
  print(end_time - start_time)
}

convert_df <- function(df, OrgDb=org.Hs.eg.db) {

  if ("ENTREZID" %in% names(df)) return(df)
  
  df[[keytype]] <- row.names(df)
  # Convert to ENTREZ ID
  # We will lose some genes here because not all IDs will be converted

  ids<-bitr(row.names(df), fromType = keytype, toType = "ENTREZID", OrgDb=OrgDb)
  df <- merge(df, ids, by = keytype, all.x = TRUE)
  print(paste("Before",nrow(df)))
  df <- na.omit(df)
  print(paste("After",nrow(df)))
  return(df)
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

df <- convert_df(df, OrgDb=OrgDb)
#print(head(df))
run_clusterProfiler(df, outfile_go, outfile_kegg, metric, overwrite=FALSE, organism.KEGG=organismKEGG, organism.GO = OrgDb) 