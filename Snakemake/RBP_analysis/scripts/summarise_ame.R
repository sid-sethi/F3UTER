args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
})
options(stringsAsFactors = FALSE)

ame_file <- args[1]
pval <- as.numeric(args[2])
outpath <- args[3]
prefix <- args[4]

x = read.table(ame_file, header=TRUE) %>%
  separate(motif_ID, c("RBP_name", "RBP_id", NA), sep=":", remove = FALSE, convert = TRUE)

data <- x %>%
  filter(adj_p.value < pval) %>%
  dplyr::select(RBP_name, RBP_id) %>%
  distinct(RBP_id,.keep_all = TRUE)

write.table(data, file = str_c(outpath,"/", prefix, "_rbpNames.txt"), quote=FALSE, row.names=FALSE, sep="\t", col.names=FALSE)

# for making the RBP table
df.table <- x %>%
  filter(adj_p.value < pval) %>%
  dplyr::select(rank, RBP_name, RBP_id, p.value, adj_p.value) %>%
  distinct(RBP_id,.keep_all = TRUE)

write.table(df.table, file = str_c(outpath, "/", prefix, "_table.txt"), quote=FALSE, row.names=FALSE, sep="\t")
