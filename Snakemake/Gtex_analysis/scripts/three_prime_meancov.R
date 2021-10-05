args<-commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(derfinder)
  library(rtracklayer)
})

options(stringsAsFactors=F, warn=-1)

data = args[1]
outPath = args[2]
tissue = args[3]
coverage_dir = args[4]


# mean coverage
meanCov <- function(start, end){
  mean(cov[start:end])
}


# data
load(data) # three_prime

er.exp <- data.frame()

for(i in c(1:22,"X","Y")){

  # per chr
  chr.er = three_prime %>% dplyr::filter(seqnames %in% i)

  if(nrow(chr.er) > 0) {

    # load exp data
    file = paste(coverage_dir, "/gtex_", tissue, "_chr", i, "_mean_cov.rda", sep="")

    if(i != 1){ rm(tissue_coverage_w_mean_normalised) }

    load(file)
    cov = tissue_coverage_w_mean_normalised$meanCoverage %>% as.numeric()

    # analyse ER regions
    result.er = chr.er %>%
      rowwise() %>%
      mutate(meanCov = meanCov(start, end)) %>%
      as.data.frame()

    er.exp <- rbind.data.frame(er.exp, result.er)

  }

}

er.export = left_join(three_prime, er.exp, by = "three_prime_utr_id") %>%
  mutate(tissue = tissue) %>%
  dplyr::select(c(three_prime_utr_id, meanCov)) %>%
  plyr::rename(c("three_prime_utr_id" = "id"))


file_name = str_c(outPath, "/", tissue, "_meancov.txt")
write.table(er.export, file = file_name , row.names = FALSE, quote = FALSE, sep="\t")
