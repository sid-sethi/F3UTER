args <- commandArgs(TRUE)
library(tidyverse)
library(rtracklayer)
library(plyranges)

options(stringsAsFactors=F)

# CNCR data from Chen et al. 2020
load("/CNCR_score/CNC_gr.rda")


# loading predictions
inFile <- args[1]
outFile <- args[2]

ER <- read.table(inFile, header=TRUE, sep="\t")
er.gr <- makeGRangesFromDataFrame(ER, keep.extra.columns = TRUE)

res <- plyranges::join_overlap_inner_within(CNC_gr, er.gr) %>%
  plyranges::group_by(id) %>%
  plyranges::summarise(cnc = round(mean(CNC),3)) %>%
  as.data.frame()

file_name <- outFile
write.table(res, file = file_name , row.names = FALSE, quote = FALSE, sep="\t")
