args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
  library(plyranges)
})

options(stringsAsFactors=F)

cncr_data <- args[1]
three_prime_data <- args[2]
outpath <- args[3]

# CNCR data
load(cncr_data)

# loading the known regions
load(three_prime_data)

# converting to Granges object
utr.gr = makeGRangesFromDataFrame(three_prime, keep.extra.columns = TRUE)
newStyle <- mapSeqlevels(seqlevels(utr.gr), "UCSC")
utr.gr <- renameSeqlevels(utr.gr, newStyle)

utr.res <- plyranges::join_overlap_inner_within(CNC_gr, utr.gr) %>%
  plyranges::group_by(three_prime_utr_id) %>%
  plyranges::summarise(cnc = round(mean(CNC),3)) %>%
  as.data.frame()

write.table(utr.res, file = str_c(outpath, "/cncr_scores_known_three_utrs.txt"), row.names = FALSE, quote = FALSE, sep="\t")
