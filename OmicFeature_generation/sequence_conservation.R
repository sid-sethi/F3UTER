library(tidyverse)
library(stringr)
library(rtracklayer)

options(stringsAsFactors=F)

# phastCons score file from UCSC
bw_path = "/PhastCons/hg38.phastCons7way.bw"

# Functions -------------------------------------------------------------------------------------------

get_conservation_score <- function(bw_path, gr, summaryFun  = "mean"){

  BigWigFile <- BigWigFile(bw_path)

  phast_cons_score <- bw_path %>% str_replace(".*/", "") %>% str_extract("phastCons.*way")

  gr_w_scores <- summary(BigWigFile, gr, size = 1L, type = summaryFun) %>% unlist()

  stopifnot((width(gr)  == width(gr_w_scores)))

  elementMetadata(gr)[[str_c(summaryFun, "_", phast_cons_score)]] <- gr_w_scores$score

  gr_w_scores <- gr

  return(gr_w_scores)

}


# loading the UTR training regions
# this code can be modified to be used on other training regions and on ERs
load("/GS.RData")

utr.gr = makeGRangesFromDataFrame(utr, keep.extra.columns = TRUE)

# adding "chr" in front of seqnames
newStyle <- mapSeqlevels(seqlevels(utr.gr), "UCSC")
utr.gr <- renameSeqlevels(utr.gr, newStyle)


utr.export = get_conservation_score(bw_path, utr.gr) %>% 
  as.data.frame() %>%
  mutate(class = "3-UTR") %>%
  select(c(three_prime_utr_id, mean_phastCons7way, class)) %>%
  plyr::rename(c("three_prime_utr_id" = "id"))


file_name = "/Features/phastCons_feat.txt"
write.table(utr.export, file = file_name, row.names = FALSE, quote = FALSE, sep="\t")
