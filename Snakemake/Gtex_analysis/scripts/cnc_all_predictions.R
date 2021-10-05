args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
  library(plyranges)
})

options(stringsAsFactors=F)

all_tissues <- args[1]
cnc_data <- args[2]
outpath <- args[3]

# CNCR data
load(cnc_data)

# loading predictions


data <- readRDS(all_tissues)

ER <- data %>%
  dplyr::filter(! tissue %in% c("brain_cerebellum", "cortex"))

er.gr <- makeGRangesFromDataFrame(ER, keep.extra.columns = TRUE)

res <- plyranges::join_overlap_inner_within(CNC_gr, er.gr) %>%
  plyranges::group_by(id) %>%
  plyranges::summarise(cnc = round(mean(CNC),3)) %>%
  as.data.frame() %>%
  dplyr::left_join(ER, by = "id")

write.table(res, file = str_c(outpath, "/cncr_all_predictions.txt") , row.names = FALSE, quote = FALSE, sep="\t")
