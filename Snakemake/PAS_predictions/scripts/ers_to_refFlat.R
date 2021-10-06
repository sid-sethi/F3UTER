args <- commandArgs(TRUE)
suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
})

er_file <- args[1]
resPath <- args[2]
prefix <- args[3]

ers= read.table(er_file, header=TRUE, sep ="\t", stringsAsFactor = FALSE)
seqlevelsStyle(ers$seqnames) <- "NCBI"

refFlat <- data.frame(
  gene_name = str_c("NM_", ers$ER_id),
  gene_id = str_c("NM_", ers$ER_id),
  seqnames = ers$seqnames,
  strand = "+",
  trs_start = ers$start - 200,
  trs_end = ers$end,
  cds_start = ers$start - 100,
  cds_end = ers$start - 1,
  exons = 1,
  exons_start = str_c(ers$start - 200, ","),
  exons_end = str_c(ers$end, ",")
)

write.table(refFlat, file = str_c(resPath, "/", prefix, "_refFlat.txt"), row.names = FALSE, quote=FALSE, col.names=FALSE, sep="\t")
