args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(BSgenome)
  library(BSgenome.Hsapiens.UCSC.hg38.masked)
  library(rtracklayer)
  library(stringr)
  library(TFBSTools)
})

options(stringsAsFactors=F)

erFile = args[1]
outPath = args[2]
tissue = args[3]

tissues = c(tissue)

# building the position frequency matrix for PAS
pfm <- PFMatrix(ID="Motif_all", name="PolyA_signals",
                bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                profileMatrix=matrix(c(9L,  9L, 0L,  11L,  9L,  12L,
                                       1L, 1L,  0L, 0L,  1L,  0L,
                                       1L,  1L,  1L,  1L, 1L,  0L,
                                       1L,  1L,  11L,  0L,  1L, 0L),
                                     byrow=TRUE, nrow=4,
                                     dimnames=list(c("A", "C", "G", "T"))
                )
)

# conversion to position weight matrix
pwm <- TFBSTools::toPWM(pfm, type="log2probratio", pseudocounts=0.8,
                        bg=c(A=0.25, C=0.25, G=0.25, T=0.25))


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


  # searching sequences with 95% match to PWM (time intensive commands)
  er.sites <- searchSeq(pwm, er.seq, min.score="95%", strand="*")


  er.relscore = relScore(er.sites)
  er.hits.len = sapply(er.relscore, length) %>% as.data.frame()
  er.hits.len$id = rownames(er.hits.len)
  er.export = er.hits.len %>%
    mutate(polyA_signal = as.factor(ifelse(. > 0, "1", .)), tissue = tissue) %>%
    select(-c(.))


  file_name = str_c(outPath, "/", tissue, "_polyA_signal.txt")
  write.table(er.export, file = file_name , row.names = FALSE, quote = FALSE, sep="\t")

  rm(er.export, er.seq, er.sites)
}
