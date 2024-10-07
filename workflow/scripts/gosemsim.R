load_godata <- function(org, ont, cache_dir="") {
  
  # Create the cache directory if it doesn't exist
  if (cache_dir == "")
    cache_dir <- "resources/.cache"
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }
  
  # Define the cache file path
  cache_file <- file.path(cache_dir, paste0("godata_", org, "_", ont, "_cache.rds"))
  
  # Check if cache exists
  if (file.exists(cache_file)) {
    message("Loading cached godata...")
    godata_obj <- readRDS(cache_file)
  } else {
    message("Cache not found. Running godata()...")
    godata_obj <- godata(org, ont = ont)
    saveRDS(godata_obj, cache_file)
  }
  
  return(godata_obj)
}

semsim <-function(depth_df, ont, qval, org, measure="Wang", savepath = "") {
    tab <- depth_df[depth_df$ONTOLOGY == ont,]
    terms <- tab$ID#tab[tab[["Combined.FDR"]] < qval, ]$ID
    terms <- unique(terms)
    
    if (length((terms))< 1) {
        print("No terms found for semantic similarity analysis")
        sim_df <- data.frame()
    } else {

        godata_ont <- load_godata(org, ont)
        sim_matrix <- mgoSim(terms, terms, semData=godata_ont, measure=measure, combine=NULL)
        sim_df <- as.data.frame(sim_matrix)
    }

    if (savepath != "") {
        dir.create(dirname(savepath), recursive = TRUE, showWarnings = FALSE)
        write.csv(sim_df, savepath)
    }

    return(sim_df)
}

semsim_from_path <- function(depth_df_path, ont, qval, org, measure="Wang", savepath="") {
    depth_df <- read.csv(depth_df_path)
    names(depth_df)[names(depth_df) == 'X'] <- 'ID'
    semsim(depth_df, ont, qval, org, measure="Wang", savepath=savepath)
}

if (!interactive()) {

  suppressMessages(library(clusterProfiler))
  suppressMessages(library(GOSemSim))
  #suppressMessages(library(org.Mm.eg.db))

  args <- commandArgs(trailingOnly = TRUE)
  print(paste("Args:",args))

  depth_df_path <- args[1]
  organismKEGG <- args[2]
  savepath <- args[3]
  qval <- args[4]

  if (organismKEGG == "hsa") {
      suppressMessages(library(org.Hs.eg.db))
      org <- "org.Hs.eg.db"
  } else if (organismKEGG == "mmu") {
      suppressMessages(library(org.Mm.eg.db))
      org <- "org.Mm.eg.db"
  } else {
      stop(paste("Organism not yet implemented:", organismKEGG))
  }

  for (subont in c("BP","CC","MF")) {
    savepath_ont <- file.path(savepath, paste0("sim_matrix_",subont,".csv"))
    semsim_from_path(depth_df_path, subont, qval, org, measure="Wang", savepath=savepath_ont)
  }

}