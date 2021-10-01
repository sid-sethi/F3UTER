args<-commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(derfinder)
  library(rtracklayer)
})

options(stringsAsFactors=F)

erFile = args[1]
outPath = args[2]
tissue = args[3]
coverage_dir = args[4]

tissues = c(tissue)


####################################################################################################################
###################################################      Functions      ############################################
####################################################################################################################

# mean coverage
meanCov <- function(start, end){
  mean(cov[start:end])
}

# percentage difference between reads at boundaries
pd <- function(start, end){
  pd = (abs(cov[start] - cov[end])/mean(c(cov[start], cov[end])))*100
  return(round(pd,3))
}

# entropy efficiency (EE)
EE <- function(start, end){

  cov.vector = cov[start:end]

  # only consider bases which have coverage > 0
  cov.vector = cov.vector[cov.vector>0]

  length = length(cov.vector)
  total_coverage = sum(cov.vector)
  p.xi = cov.vector/total_coverage
  log.p.xi = log(p.xi)

  # final equation
  ee = -(sum(p.xi*log.p.xi))/log(length)
  return(ee)
}

####################################################################################################################
####################################################################################################################
####################################################################################################################


for(tissue in tissues){

  print(tissue)

  #ER data
  er = read.table(erFile, header=T, sep="\t")

  er.exp <- data.frame()

  for(i in c(1:22,"X","Y")){

    #print(i)

    # per chr
    chr.er = er %>% dplyr::filter(seqnames %in% paste("chr", i, sep=""))

    if(nrow(chr.er) > 0) {

      # load exp data
      file = paste(coverage_dir, "/", tissue, "_", i, "_mean_cov.rda", sep="")

      if(i != 1){ rm(tissue_coverage_w_mean_normalised) }

      load(file)
      cov = tissue_coverage_w_mean_normalised$meanCoverage %>% as.numeric()

      # analyse ER regions
      result.er = chr.er %>%
        rowwise() %>%
        mutate(meanCov = meanCov(start, end),
               pd = pd(start, end),
               EE = EE(start, end)) %>%
        as.data.frame()

      er.exp <- rbind.data.frame(er.exp, result.er)

    }

  }

  er.export = left_join(er, er.exp, by = "ER_id") %>%
    mutate(tissue = tissue) %>%
    dplyr::select(c(ER_id, meanCov, pd, EE, tissue)) %>%
    plyr::rename(c("ER_id" = "id"))


  file_name = str_c(outPath, "/", tissue, "_exp_feat.txt")
  write.table(er.export, file = file_name , row.names = FALSE, quote = FALSE, sep="\t")

  rm(er.export, er)

}
