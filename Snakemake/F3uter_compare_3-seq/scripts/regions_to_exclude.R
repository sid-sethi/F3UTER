args<-commandArgs(TRUE)

suppressPackageStartupMessages({
  library(rtracklayer)
  library(tidyverse)
})

cat("\n***************************\n")

# ER data to add to exculsion
# requires both 3' and 5' Intergenic ERs

three_prime <- args[1]
five_prime <- args[2]
outpath = args[3]
tissue = args[4]
gtf_file = args[5]

print(str_c(Sys.time(), " - processing ER data for masking"))

files <- c(three_prime, five_prime)

er_data = files %>%
  map_dfr(read.table, header=TRUE) %>%
  dplyr::select(seqnames:strand)


print(str_c(Sys.time(), " - processing GTF data for masking"))
gtf <- rtracklayer::import(gtf_file)
gtf = keepStandardChromosomes(gtf, species = "Homo_sapiens", pruning.mode="coarse")
newStyle <- mapSeqlevels(seqlevels(gtf), "UCSC")
gtf <- renameSeqlevels(gtf, newStyle)

# extracting coordinates of all genes
genes = gtf[gtf$type %in% c("gene")]
genes.m = GenomicRanges::reduce(genes, ignore.strand=FALSE) %>% as.data.frame()


print(str_c(Sys.time(), " - merging above datasets"))
data_to_exclude = dplyr::bind_rows(genes.m, er_data)

write.table(data_to_exclude, file = str_c(outpath, "/", tissue, "_regions_to_exclude.txt"), sep="\t", row.names=FALSE, col.names=FALSE, quote = FALSE)
