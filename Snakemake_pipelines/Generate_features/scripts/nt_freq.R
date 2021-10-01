args<-commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(BSgenome)
  library(BSgenome.Hsapiens.UCSC.hg38.masked)
  library(rtracklayer)
  library(stringr)
})

options(stringsAsFactors=F)


erFile = args[1]
outPath = args[2]
tissue = args[3]

tissues = c(tissue)

for(tissue in tissues){

  print(tissue)

  #ER data
  er = read.table(erFile, header=T, sep="\t") %>%
    dplyr::select(-c(strand))

  # converting to Granges object
  er.gr = makeGRangesFromDataFrame(er, keep.extra.columns = TRUE)

  # Extract sequences (using the masked version)
  er.seq = BSgenome::getSeq(Hsapiens, er.gr)
  er.seq@ranges@NAMES <- er.gr$ER_id

  # calculating nucleotide frequency
  er.nf = letterFrequency(er.seq,c("A","T","G","C"), as.prob = TRUE) %>% as.data.frame()
  # di-nucleotide
  dinf = data.frame(AA = (str_count(er.seq, "AA")/er.seq@ranges@width),
                 TT = (str_count(er.seq, "TT")/er.seq@ranges@width),
                 GG = (str_count(er.seq, "GG")/er.seq@ranges@width),
                 CC = (str_count(er.seq, "CC")/er.seq@ranges@width),
                 AT = (str_count(er.seq, "AT")/er.seq@ranges@width),
                 AG = (str_count(er.seq, "AG")/er.seq@ranges@width),
                 AC = (str_count(er.seq, "AC")/er.seq@ranges@width),
                 TA = (str_count(er.seq, "TA")/er.seq@ranges@width),
                 TG = (str_count(er.seq, "TG")/er.seq@ranges@width),
                 TC = (str_count(er.seq, "TC")/er.seq@ranges@width),
                 GA = (str_count(er.seq, "GA")/er.seq@ranges@width),
                 GT = (str_count(er.seq, "GT")/er.seq@ranges@width),
                 GC = (str_count(er.seq, "GC")/er.seq@ranges@width),
                 CA = (str_count(er.seq, "CA")/er.seq@ranges@width),
                 CT = (str_count(er.seq, "CT")/er.seq@ranges@width),
                 CG = (str_count(er.seq, "CG")/er.seq@ranges@width)
                 )

  er.export = bind_cols(er, er.nf, dinf) %>%
    dplyr::select(ER_id:CG) %>%
    mutate(tissue = tissue) %>%
    plyr::rename(c("ER_id" = "id"))

  file_name = str_c(outPath, "/", tissue, "_nt_freq.txt")
  write.table(er.export, file = file_name , row.names = FALSE, quote = FALSE, sep="\t")

  rm(er.nf, er.export, er.seq, er.gr, dinf)

}
