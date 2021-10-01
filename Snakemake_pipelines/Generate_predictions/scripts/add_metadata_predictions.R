args<-commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
})

options(stringsAsFactors=F)

predictions = args[1]
erFile = args[2]
outpath = args[3]
tissue = args[4]

tissues <- c(tissue)


for(tissue in tissues){

  print(tissue)

  #ER data
  er = read.table(erFile, header=T, sep="\t")

  ## prediction data
  pred = read.table(predictions, header=TRUE)

  data = left_join(pred, er, by = c("id" = "ER_id"))

  write.table(data, file = str_c(outpath, "/", tissue, "_predictions_meta.txt"), row.names=FALSE, col.names = TRUE, quote = FALSE, sep="\t")

  # saving data as a bed format
  bed = data %>%
    mutate(score = 0, bed_strand = ".", thickstart = start, thickend = end, color = ifelse(predicted.prob>0.60, "0,175,187", "231,184,0")) %>%
    dplyr::select(seqnames, start, end, id, score, bed_strand, thickstart, thickend, color)

  write.table(bed, file = str_c(outpath, "/", tissue, "_predictions.bed"), row.names=FALSE, col.names = FALSE, quote = FALSE, sep="\t")

  rm(pred, data)

}
