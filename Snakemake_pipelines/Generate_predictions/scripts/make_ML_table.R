args<-commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
})

options(stringsAsFactors=F)

path = args[1] # path to directory containing features
outpath = args[2]
tissue = args[3]

tissues = c(tissue)


for(tissue in tissues){

  print(tissue)

  pa_signal = str_c(path, "/", tissue, "_polyA_signal.txt")
  ntf = str_c(path, "/", tissue, "_nt_freq.txt")
  exp = str_c(path, "/", tissue, "_exp_feat.txt")
  phastCons = str_c(path, "/", tissue, "_phastcons.txt")
  dsp = str_c(path, "/", tissue, "_structural_feat.txt")
  repeats = str_c(path, "/", tissue, "_repeats.txt")


  df1 = read.table(pa_signal, header=TRUE, sep="\t") %>% dplyr::select(c(id, polyA_signal)) %>% mutate(polyA_signal = as.factor(polyA_signal))
  df2 = read.table(exp, header=TRUE, sep="\t") %>% dplyr::select(-c(id, tissue, meanCov))
  df3 = read.table(phastCons, header=TRUE, sep="\t") %>% dplyr::select(-c(id, tissue))
  df4 = read.table(ntf, header=TRUE, sep="\t") %>% dplyr::select(-c(id, tissue))
  df5 = read.table(dsp, header=TRUE, sep="\t") %>% dplyr::select(-c(id, tissue))
  df6 = read.table(repeats, header=TRUE, sep="\t") %>% dplyr::select(-c(id, tissue)) %>% plyr::rename(c("fracOverlap" = "RepeatsOverlap"))

  data = cbind(df1, df2, df3, df4, df5, df6) %>%
    tidyr::drop_na() %>%
    remove_rownames() %>%
    column_to_rownames("id")

  file = str_c(outpath, "/", tissue, "_ml_table.rds")
  saveRDS(data, file= file)

  rm(data)
}
