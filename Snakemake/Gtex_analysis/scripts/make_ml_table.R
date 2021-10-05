args<-commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
})

options(stringsAsFactors=F)

path = args[1]
outpath = args[2]
tissue <- args[3]



  print(str_c(Sys.time(), "-", tissue, " - compiling ML table"))

  pa_signal = paste(path, "/PAS/", tissue, "_polyA_signal.txt", sep="")
  ntf = paste(path, "/Nt_freq/", tissue, "_nt_freq.txt", sep="")
  exp = paste(path, "/Expression/", tissue, "_exp_feat.txt", sep="")
  phastCons = paste(path, "/Phastcons/", tissue, "_phastcons.txt", sep="")
  dsp = paste(path, "/Structural_properties/", tissue, "_structural_feat.txt", sep="")
  repeats = paste(path, "/Transposons/", tissue, "_repeats.txt", sep="")


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

  file = str_c(outpath, "/table_", tissue, ".rds")
  saveRDS(data, file= file)
