args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(BSgenome)
  library(BSgenome.Mmusculus.UCSC.mm10)
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

  print(str_c(Sys.time(), " - scanning PAS occurences"))

  #ER data
  er = read.table(erFile, header=T, sep="\t")

  # converting to Granges object
  er.gr = makeGRangesFromDataFrame(er, keep.extra.columns = TRUE)

  # adding "chr" in front of seqnames
  newStyle <- mapSeqlevels(seqlevels(er.gr), "UCSC")
  er.gr <- renameSeqlevels(er.gr, newStyle)

  # Extract sequences
  er.seq = BSgenome::getSeq(Mmusculus, er.gr)
  er.seq@ranges@NAMES <- er.gr$id


  # searching sequences with 95% match to PWM (time intensive commands)
  er.sites <- searchSeq(pwm, er.seq, min.score="95%", strand="*")


  er.relscore = relScore(er.sites)
  er.hits.len = sapply(er.relscore, length) %>% as.data.frame()
  er.hits.len$id = rownames(er.hits.len)
  er.export = er.hits.len %>%
    mutate(polyA_signal = as.factor(ifelse(. > 0, "1", .)), class = er$class) %>%
    select(-c(.))


  file_name = str_c(outPath, "/", tissue, "_polyA_signal.txt")
  write.table(er.export, file = file_name , row.names = FALSE, quote = FALSE, sep="\t")

  rm(er.export, er.seq, er.sites)
}
