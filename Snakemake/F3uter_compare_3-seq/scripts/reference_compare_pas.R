args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
})

three_prime_data <- args[1]
outpath <- args[2]
tissue <- args[3]
pas_dir <- args[4]

######################################################
######################################################

print(str_c(Sys.time(), " - ", tissue))
print(str_c(Sys.time(), " - processing PAS data"))

pas_files = list.files(pas_dir, pattern = ".bed", full.names = TRUE)

pas_gr = GRanges()
for(f in pas_files){
  print(f)
  pa = read.table(f,header=F,sep="\t")
  colnames(pa) = c("chr","start","end","id","exp_thisSample","strand", "perc_sample_support", "number_protocol_support", "avg_exp_allSamples", "cluster_annotation","pas")
  gr = makeGRangesFromDataFrame(pa, keep.extra.columns = TRUE)
  pas_gr <- c(pas_gr,gr)
}


pac <- IRanges::reduce(pas_gr)
pac = keepStandardChromosomes(pac, species = "Homo_sapiens", pruning.mode="coarse")


######### three prime data #########
print(str_c(Sys.time(), " - processing known three prime UTRs"))
load(three_prime_data)

utr.gr = makeGRangesFromDataFrame(three_prime, keep.extra.columns = TRUE)

utr.ov = IRanges::subsetByOverlaps(utr.gr, pac, ignore.strand = FALSE, maxgap= -1, type = "any") %>% length()
utr.ov.perc <- round((utr.ov/length(utr.gr))*100,2)

cat(
  "Overlap of known 3'UTRs with poly(A) site clusters in ", tissue, " = ", utr.ov.perc, "\n",
  file = str_c(outpath, "/", tissue, "_reference_compare_pas.txt")
)
