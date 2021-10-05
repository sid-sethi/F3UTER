args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
})

options(stringsAsFactor = FALSE, warn=-1)

gtf_file <- args[1]
er_file <- args[2]
outpath <- args[3]
prefix <- args[4]
dataset_name <- args[5]


er <- readRDS(er_file)
er.gr <- er %>% makeGRangesFromDataFrame(., keep.extra.columns = TRUE)


######################################
############ Refseq processing ##########
######################################

ref <- read.table(gtf_file, header=FALSE, sep="\t")
colnames(ref) <- c("seqnames", "start", "end", "id", "tmp", "strand")
utr.gr <- ref %>% makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

###################################


hits = IRanges::findOverlaps(er.gr, utr.gr, ignore.strand = TRUE, minoverlap = 0L)
intersect = pintersect(er.gr[queryHits(hits)], utr.gr[subjectHits(hits)])
percentOverlap <- width(intersect) / width(er.gr[queryHits(hits)])
hits <- hits[percentOverlap >= 0.5]


er.hits <- er.gr[unique(queryHits(hits))] %>% as.data.frame()

  df <- data.frame(
    count = nrow(er.hits),
    total = nrow(er)
  ) %>% mutate(
    perc = (count/total)*100,
    dataset = dataset_name
  )


  write.table(er.hits, str_c(outpath, "/", prefix, "_ers_refseq.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
  write.table(df, str_c(outpath, "/", prefix, "_overlap_refseq.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
