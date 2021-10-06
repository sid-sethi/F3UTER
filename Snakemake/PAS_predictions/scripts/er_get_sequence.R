args <- commandArgs(TRUE)
suppressPackageStartupMessages({
  library(tidyverse)
  library(BSgenome)
  library(BSgenome.Hsapiens.UCSC.hg38.masked)
})

er_file <- args[1]
resPath <- args[2]
prefix <- args[3]

x = read.table(er_file, header=TRUE, sep ="\t", stringsAsFactor = FALSE)

# converting to Granges object
x.gr = makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)

# adding "chr" in front of seqnames
newStyle <- mapSeqlevels(seqlevels(x.gr), "UCSC")
x.gr <- renameSeqlevels(x.gr, newStyle)

# Extract sequences (using the masked version)
x.seq = BSgenome::getSeq(Hsapiens, x.gr)
x.seq@ranges@NAMES <- x.gr$ER_id

df <- x.seq %>%
  as.data.frame() %>%
  tibble::rownames_to_column("er_id")

write.table(df, file = str_c(resPath, "/", prefix, "_seq.txt"), row.names = FALSE, quote=FALSE, sep="\t")
