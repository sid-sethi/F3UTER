library(tidyverse)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38.masked)


# load GTEx ER predictions
# Files: all_tissues_mergedData.rds and tissue-specific categorised files
all <- readRDS("/PredictionApp_AllTissues/all_tissues_appData.rds") %>%
  dplyr::select(-c(strand))

all.pos <- all %>%  dplyr::filter(predicted.prob > 0.60)
all.neg <- all %>%  dplyr::filter(predicted.prob <= 0.60)

# converting to Granges object
pos.gr = makeGRangesFromDataFrame(all.pos, keep.extra.columns = TRUE)
neg.gr = makeGRangesFromDataFrame(all.neg, keep.extra.columns = TRUE)

# adding "chr" in front of seqnames
newStyle <- mapSeqlevels(seqlevels(pos.gr), "UCSC")
pos.gr <- renameSeqlevels(pos.gr, newStyle)

newStyle <- mapSeqlevels(seqlevels(neg.gr), "UCSC")
neg.gr <- renameSeqlevels(neg.gr, newStyle)

# Extract sequences (using the masked version)
pos.seq = BSgenome::getSeq(Hsapiens, pos.gr)
pos.seq@ranges@NAMES <- pos.gr$id
neg.seq = BSgenome::getSeq(Hsapiens, neg.gr)
neg.seq@ranges@NAMES <- neg.gr$id

writeXStringSet(pos.seq, "/Intergenic_3/Fasta/all.pos.fasta", format="fasta")
writeXStringSet(neg.seq, "/Intergenic_3/Fasta/all.neg.fasta", format="fasta")
