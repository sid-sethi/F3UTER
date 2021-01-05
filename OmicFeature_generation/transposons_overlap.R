library(tidyverse)
library(rtracklayer)

options(stringsAsFactors=F)


# load repeats/transposons data
x = "/Repeats/hg38.repeatMasker.mod.fa.out"
d = read.table(x, header=FALSE)
colnames(d) = c("seqnames", "start", "end", "repeat_type")
repeats = d %>% dplyr::filter(!repeat_type %in% c("RNA", "rRNA", "scRNA", "snRNA", "srpRNA", "tRNA", "Unknown", "Simple_repeat", "Satellite/telo", "Low_complexity", "Satellite", "Satellite/acro", "Satellite/centr"))

rep.gr = makeGRangesFromDataFrame(repeats, keep.extra.columns = TRUE)
rep.gr = keepStandardChromosomes(rep.gr, species = "Homo_sapiens", pruning.mode="coarse")


# loading the UTR training regions
# this code can be modified to be used on other training regions and on ERs
load("/GS.RData")

# converting to Granges object
utr.gr = makeGRangesFromDataFrame(utr, keep.extra.columns = TRUE)


# adding "chr" in front of seqnames
newStyle <- mapSeqlevels(seqlevels(utr.gr), "UCSC")
utr.gr <- renameSeqlevels(utr.gr, newStyle)


hits = IRanges::findOverlaps(utr.gr, rep.gr, ignore.strand = TRUE)
ovr = pintersect(utr.gr[queryHits(hits)], rep.gr[subjectHits(hits)])
ovr$fracOverlap <- round( width(ovr) / width(utr.gr[queryHits(hits)]) , 3)

overlaps <- ovr %>% as.data.frame() %>%
  group_by(three_prime_utr_id) %>%
  summarise(fracOverlap = sum(fracOverlap))

utr.export <- dplyr::left_join(utr, overlaps, by = "three_prime_utr_id") %>%
  mutate(fracOverlap = replace_na(fracOverlap, 0)) %>%
  mutate(fracOverlap = ifelse(fracOverlap > 1, 1, fracOverlap)) %>%
  mutate(class = "3-UTR") %>%
  dplyr::select(c(three_prime_utr_id:class)) %>%
  plyr::rename(c("three_prime_utr_id" = "id"))


write.table(utr.export, file = "/Features/repeats.txt", row.names = FALSE, quote = FALSE, sep="\t")
