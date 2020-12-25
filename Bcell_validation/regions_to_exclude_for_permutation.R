args<-commandArgs(TRUE)

library(rtracklayer)
library(tidyverse)

cat("\n***************************\n")

# ER data to add to exculsion
# 3' and 5' Intergenic ERs (need both)

er_path = args[1]
out_path = args[2]
tissue = args[3]

print(str_c(Sys.time(), " - processing ER data for masking"))

files = c(
  str_c(er_path, "/3_prime/", tissue, ".txt"),
  str_c(er_path, "/5_prime/", tissue, ".txt")
)

er_data = files %>%
  map_dfr(read.table, header=TRUE) %>%
  dplyr::select(seqnames:strand)


print(str_c(Sys.time(), " - processing GTF data for masking"))
gtf <- rtracklayer::import("/Homo_sapiens.GRCh38.94.gtf")
gtf = keepStandardChromosomes(gtf, species = "Homo_sapiens", pruning.mode="coarse")
newStyle <- mapSeqlevels(seqlevels(gtf), "UCSC")
gtf <- renameSeqlevels(gtf, newStyle)

# extracting coordinates of all genes
genes = gtf[gtf$type %in% c("gene")]
genes.m = GenomicRanges::reduce(genes, ignore.strand=FALSE) %>% as.data.frame()


print(str_c(Sys.time(), " - merging above datasets"))
data_to_exclude = dplyr::bind_rows(genes.m, er_data)

write.table(data_to_exclude, file = str_c(out_path, "/", tissue, "_regions_to_exclude.txt"), sep="\t", row.names=FALSE, col.names=FALSE, quote = FALSE)
