run_clusterProfiler <- function(df, outfile,
                                metric,
                                library_,
                                organism.KEGG,
                                organism.GO, 
                                overwrite=FALSE, 
                                minGSSize = 10,
                                maxGSSize = 500,
                                seed=1234,
                                keytype_gmt = "SYMBOL") 
{
  set.seed(seed)

  print(paste0("ClusterProfiler target", outfile))

  if (file.exists(outfile) && !overwrite) {
    print("Existing files not overwritte, skipping")
    return
  }

  start_time <- Sys.time()

  geneList <- df[[metric]]

  if ((endsWith(library_, ".gmt") && !file.exists(outfile)) || overwrite) {
    if (!file.exists(library_)) {
      library_ = file.path("./resources/Ontologies",library_)
      if (!file.exists(library_))
        stop(paste0("gmt file not found:", library_))
    }
    print(paste0("Running with custom gmt file:", library_))
    gmt <- read.gmt(library_)
    gmt$gene <- toupper(gmt$gene)
    names(geneList) <- toupper(df[[keytype_gmt]])
    geneList = sort(geneList, decreasing = TRUE)
    ego3 <- GSEA(geneList     = geneList,
              TERM2GENE = gmt,
              minGSSize    = minGSSize,
              maxGSSize    = maxGSSize,
              pvalueCutoff = 1,
              eps = 0,
              seed = TRUE,
              verbose = FALSE)

    write.csv(ego3,outfile)
    print(paste("Wrote ClusterProfiler output to:", outfile))

  } else if ((library_=="GO" && !file.exists(outfile)) || overwrite) {

    names(geneList) <- df$ENTREZID
    geneList = sort(geneList, decreasing = TRUE)

    print("Running GO...")
    ego3 <- gseGO(geneList     = geneList,
                  OrgDb        = organism.GO,
                  ont          = "ALL", ## CC MF BP
                  minGSSize    = minGSSize,
                  maxGSSize    = maxGSSize,
                  pvalueCutoff = 1,
                  eps = 0,
                  seed = TRUE,
                  verbose = FALSE)
    write.csv(ego3,outfile)
    print(paste("Wrote GO to:", outfile))

  } else if ((library_=="KEGG" && !file.exists(outfile)) || overwrite)  {

    names(geneList) <- df$ENTREZID
    geneList = sort(geneList, decreasing = TRUE)

    print("Running KEGG...")
    kegg <- gseKEGG(geneList     = geneList,
                  organism        = organism.KEGG,
                  minGSSize    = minGSSize,
                  maxGSSize    = maxGSSize,
                  pvalueCutoff = 1,
                  eps = 0,
                  seed = TRUE,
                  verbose = FALSE)
    write.csv(kegg,outfile)
    print(paste("Wrote KEGG to:", outfile))
  } else {
    stop(paste0("Invalid library:", library_))
  }

  end_time <- Sys.time()
  print(end_time - start_time)
}

convert_df <- function(df, keytype_in, keytype_target, gene_converter_file) {

  if (keytype_target %in% names(df)) return(df)
  df[[keytype_in]] <- row.names(df)
  ids <- read.csv(gene_converter_file)
  df <- merge(df, ids, by = keytype_in, all.x = TRUE)
  print(paste("Before",nrow(df)))
  df <- na.omit(df)
  print(paste("After",nrow(df)))
  return(df)
}

if (!interactive()) {

  suppressMessages(library(clusterProfiler))

  args <- commandArgs(trailingOnly = TRUE)
  #print(paste("Args:",args))

  input_file <- args[1]
  organismKEGG <- args[2]
  gene_converter_file <- args[3]
  metric <- args[4]
  library_ <- args[5] # either "GO", "KEGG", or path to gmt file
  outfile <- args[6]
  keytype <- args[7] # keytype of input file
  keytype_gmt <- args[8] # keytype of gmt file

  print(paste("Reading clusterProfiler input:", input_file))

  #metric <- strsplit(basename(outfile), "syn.clusterProfiler\\.")[[1]][2]
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

  if (organismKEGG == "hsa") {
      suppressMessages(library(org.Hs.eg.db))
      OrgDb <- org.Hs.eg.db
  } else if (organismKEGG == "mmu") {
      suppressMessages(library(org.Mm.eg.db))
      OrgDb <- org.Mm.eg.db
  } else {
      stop(paste("Organism not yet implemented:", organismKEGG))
  }

  df <- convert_df(df, keytype_in=keytype, gene_converter_file, keytype_target=keytype_gmt)
  #print(head(df))
  run_clusterProfiler(df, outfile, metric, library_, overwrite=FALSE, organism.KEGG=organismKEGG, organism.GO = OrgDb, keytype_gmt=keytype_gmt) 

}