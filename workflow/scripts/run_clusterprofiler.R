library(clusterProfiler)
library(org.Hs.eg.db)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
outfile_go <- args[2]
outfile_kegg <- args[3]
metric <- args[4]

print(paste("Reading clsuterProfiler input:", input_file))

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
                                metric, overwrite=FALSE, 
                                organism.KEGG="hsa",
                                organism.GO = org.Hs.eg.db, seed=123) 
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

  if (!file.exists(outfile_go) || overwrite) {

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
  }

  if (!file.exists(outfile_kegg) || overwrite) {

    kegg <- gseKEGG(geneList     = geneList,
                  organism        =  organism.KEGG,
                  minGSSize    = 10,
                  maxGSSize    = 500,
                  pvalueCutoff = 1,
                  eps = 0,
                  seed = TRUE,
                  verbose = FALSE)
    write.csv(kegg,outfile_kegg)
  }

  end_time <- Sys.time()
  print(end_time - start_time)
}

convert_df <- function(df, OrgDb=org.Hs.eg.db) {

  if ("ENTREZID" %in% names(df)) return(df)
  
  df$ENSEMBL <- row.names(df)
  # Convert to ENTREZ ID
  # We will lose some genes here because not all IDs will be converted

  ids<-bitr(row.names(df), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=OrgDb)
  df <- merge(df, ids, by = "ENSEMBL", all.x = TRUE)
  print(paste("Before",nrow(df)))
  df <- na.omit(df)
  print(paste("After",nrow(df)))
  return(df)
}

df <- convert_df(df, OrgDb=org.Hs.eg.db)
run_clusterProfiler(df, outfile_go, outfile_kegg, metric, overwrite=FALSE, organism.KEGG="hsa", organism.GO = org.Hs.eg.db) 