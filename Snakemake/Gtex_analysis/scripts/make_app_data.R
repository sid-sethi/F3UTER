args<-commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(rtracklayer)
})

options(stringsAsFactors=F)

pred_file = args[1]
erPath = args[2]
outpath = args[3]
tissue = args[4]


  print(str_c(Sys.time(), " - ", tissue))

  #ER data
  er = read.table(paste(erPath, "/", tissue, ".txt", sep=""), header=T, sep="\t")

  ## prediction data
  pred = read.table(pred_file, header=TRUE)

  data = left_join(pred, er, by = c("id" = "ER_id"))

  saveRDS(data, file = str_c(outpath, "/", tissue, "_appData.rds", sep=""))
