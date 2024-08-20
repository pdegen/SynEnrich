args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

library(clusterProfiler)
library(org.Hs.eg.db)


metrics <- c("neg_signed_logpval", "logFC")

print(input_file)
df <- read.csv(input_file, row.names = 1)

for (metric in metrics) {
    # Check if the metric is in the columns
    if (!(metric %in% colnames(df))) {
        if (metric == "neg_signed_logpval") {
            message(paste("Adding", metric, "to df"))
            df$neg_signed_logpval <- -sign(df$logFC) * log10(df$PValue)
        } else {
            stop(paste("Metric", metric, "not in columns!"))
        }
    }
}

print(head(df))

write.csv(as.data.frame(c(0, 0, 0, 0)), file = output_file, row.names = FALSE)
