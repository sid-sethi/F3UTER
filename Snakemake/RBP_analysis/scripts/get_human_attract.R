args<-commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
})

options(stringsAsFactors = FALSE, warn=-1)

db <- args[1]
outpath <- args[2]

x = read.table(db, header=TRUE, sep="\t")

h = x %>% dplyr::filter(Organism %in% "Homo_sapiens", Len > 6, Score %in% "1.000000**") %>%
  group_by(Gene_id) %>%
  slice(which.max(Len))

write.table(h, str_c(outpath, "/attract_db_human.txt"), quote = FALSE, row.names=FALSE, sep="\t")
