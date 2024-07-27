args <- commandArgs(trailingOnly=TRUE)
input_file <- args[1]
output_file <- args[2]

library(clusterProfiler)

write.csv(as.data.frame(c(0,0,0,0)), file=output_file, row.names=FALSE)