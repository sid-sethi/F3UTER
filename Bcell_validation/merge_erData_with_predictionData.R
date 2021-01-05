args<-commandArgs(TRUE) # 1st  = ML output folder, 2nd = erPath
library(tidyverse)
library(stringr)
library(rtracklayer)

options(stringsAsFactors=F)

path = args[1]
erPath = args[2]

if(!dir.exists(paste(path, "/PredictionMerge", sep=""))){
  system(paste("mkdir -m a=rwx ",path, "/PredictionMerge", sep=""))
}

# new output path
out.path = paste(path, "/PredictionMerge", sep="")

tissues = c(
  "GCB1", "GCB2"
)


for(tissue in tissues){

  print(tissue)

  #ER data
  er = read.table(paste(erPath, "/", tissue, ".txt", sep=""), header=T, sep="\t")

  ## prediction data
  pred = read.table(paste(path, "/Prediction/", tissue,".pred.txt", sep=""), header=TRUE)

  data = left_join(pred, er, by = c("id" = "ER_id"))

  saveRDS(data, file = paste(out.path, "/", tissue, "_mergedData.rds", sep=""))

  # saving data as a bed format
  bed = data %>%
    mutate(score = 0, bed_strand = ".", thickstart = start, thickend = end, color = ifelse(predicted.prob>0.60, "0,175,187", "231,184,0")) %>%
    dplyr::select(seqnames, start, end, id, score, bed_strand, thickstart, thickend, color)

  write.table(bed, file = paste(out.path, "/", tissue, ".bed", sep=""), row.names=FALSE, col.names = FALSE, quote = FALSE, sep="\t")

  rm(pred, data)

}
