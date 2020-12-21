library(tidyverse)
library(scales)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38.masked)
library(rtracklayer)

options(stringsAsFactors=F)

# loading the UTR training regions
# this code can be modified to be used on other training regions and on ERs
load("/GS.RData")

# converting to Granges object
utr.gr = makeGRangesFromDataFrame(utr, keep.extra.columns = TRUE)

# adding "chr" in front of seqnames
newStyle <- mapSeqlevels(seqlevels(utr.gr), "UCSC")
utr.gr <- renameSeqlevels(utr.gr, newStyle)

# Extract sequences (using the masked version)
utr.seq = BSgenome::getSeq(Hsapiens, utr.gr)
utr.seq@ranges@NAMES <- utr.gr$three_prime_utr_id

# calculating nucleotide frequency
utr.nf = letterFrequency(utr.seq,c("A","T","G","C"), as.prob = TRUE) %>% as.data.frame()
utr.dinf = data.frame(AA = (str_count(utr.seq, "AA")/utr.seq@ranges@width),
               TT = (str_count(utr.seq, "TT")/utr.seq@ranges@width),
               GG = (str_count(utr.seq, "GG")/utr.seq@ranges@width),
               CC = (str_count(utr.seq, "CC")/utr.seq@ranges@width),
               AT = (str_count(utr.seq, "AT")/utr.seq@ranges@width),
               AG = (str_count(utr.seq, "AG")/utr.seq@ranges@width),
               AC = (str_count(utr.seq, "AC")/utr.seq@ranges@width),
               TA = (str_count(utr.seq, "TA")/utr.seq@ranges@width),
               TG = (str_count(utr.seq, "TG")/utr.seq@ranges@width),
               TC = (str_count(utr.seq, "TC")/utr.seq@ranges@width),
               GA = (str_count(utr.seq, "GA")/utr.seq@ranges@width),
               GT = (str_count(utr.seq, "GT")/utr.seq@ranges@width),
               GC = (str_count(utr.seq, "GC")/utr.seq@ranges@width),
               CA = (str_count(utr.seq, "CA")/utr.seq@ranges@width),
               CT = (str_count(utr.seq, "CT")/utr.seq@ranges@width),
               CG = (str_count(utr.seq, "CG")/utr.seq@ranges@width)
               )
utr.export = bind_cols(utr, utr.nf, utr.dinf) %>% mutate(class = "3-UTR") %>%
  dplyr::select(c(three_prime_utr_id:class)) %>%
  plyr::rename(c("three_prime_utr_id" = "id"))


write.table(utr.export, file = "/Features/nt_di_freq.txt", row.names = FALSE, quote = FALSE, sep="\t")
