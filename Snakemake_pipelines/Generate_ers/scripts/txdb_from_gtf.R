args <- commandArgs(TRUE)
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicFeatures)
  library(rtracklayer)
})

gtf_file <- args[1]
chrom_len <- args[2]
outpath <- args[3]

gtf =rtracklayer::import(gtf_file)
gtf = keepStandardChromosomes(gtf, species = "Homo_sapiens", pruning.mode="coarse")
newStyle <- mapSeqlevels(seqlevels(gtf), "UCSC")
gtf <- renameSeqlevels(gtf, newStyle)

chromInfo = read.table(chrom_len)

chrom.df = chromInfo %>% filter(V1 %in% seqnames(gtf)) %>%
  plyr::rename(c("V1" = "seqnames", "V2" = "length")) %>%
  dplyr::select(-V3)

genome_build = "hg38"

seqlengths(gtf) <- chrom.df$length
genome(gtf) <- genome_build

gtf_TxDb <- makeTxDbFromGRanges(gtf)
saveDb(gtf_TxDb, file = str_c(outpath, "/gtf_TxDb.Rsqlite"))
