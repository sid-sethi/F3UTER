args<-commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(rtracklayer)
})

options(stringsAsFactors=F)

erFile = args[1]
outPath = args[2]
tissue = args[3]
bw_path = args[4]

tissues = c(tissue)

####################################################################################################################
###################################################      Functions      ############################################
####################################################################################################################

get_conservation_score <- function(bw_path, gr, summaryFun  = "mean"){

  BigWigFile <- BigWigFile(bw_path)

  phast_cons_score <- bw_path %>% str_replace(".*/", "") %>% str_extract("phastCons.*way")

  gr_w_scores <- summary(BigWigFile, gr, size = 1L, type = summaryFun) %>% unlist()

  stopifnot((width(gr)  == width(gr_w_scores)))

  elementMetadata(gr)[[str_c(summaryFun, "_", phast_cons_score)]] <- gr_w_scores$score

  gr_w_scores <- gr

  return(gr_w_scores)

}

####################################################################################################################
####################################################################################################################
####################################################################################################################


for(tissue in tissues){

  print(tissue)

  #ER data
  er = read.table(erFile, header=T, sep="\t")

  # converting to Granges object
  er.gr = makeGRangesFromDataFrame(er, keep.extra.columns = TRUE)

  er.export = get_conservation_score(bw_path, er.gr) %>%
    as.data.frame() %>%
    mutate(tissue = tissue) %>%
    select(c(ER_id, mean_phastCons7way, tissue)) %>%
    plyr::rename(c("ER_id" = "id"))

  file_name = str_c(outPath, "/", tissue, "_phastcons.txt")
  write.table(er.export, file = file_name , row.names = FALSE, quote = FALSE, sep="\t")

  rm(er.export, er.gr)
}
