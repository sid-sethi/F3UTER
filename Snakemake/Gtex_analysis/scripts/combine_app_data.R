args<-commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
})

options(stringsAsFactors=F)

dir <- args[1]
outpath = args[2]


files <- list.files(dir, pattern = "_appData.rds", full.names = TRUE)

res_all_tissues = data.frame()

for(f in files){

  print(str_c(Sys.time(), " - ", f))

  x <- readRDS(f)

  res_all_tissues = rbind(res_all_tissues, x)

  rm(x)

}

saveRDS(res_all_tissues, file = str_c(outpath, "/all_tissues_appData.rds"))
