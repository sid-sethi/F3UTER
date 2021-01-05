args<-commandArgs(TRUE) # 1st  =  output folder (dir ML_tables will be created)

library(tidyverse)

options(stringsAsFactors=F)

path1 = "/Features_threePrime"
path2 = "/Features_ICE"
path3 = "/Features_ncRNA/lncRNA"
path4 = "/Features_ncRNA/ncRNA"
path5 = "/Features_ncRNA/pseudoGene"
path6 = "/Features_fivePrime"
out = args[1]


if(!dir.exists(paste(out, "/ML_tables", sep=""))){
  system(paste("mkdir -m a=rwx ",out, "/ML_tables", sep=""))
}

# output path
out.path = paste(out, "/ML_tables", sep="")


pa_signal = c(
  paste(path1, "/polyA_signal.txt", sep=""),
  paste(path2, "/polyA_signal.txt", sep=""),
  paste(path3, "/polyA_signal.txt", sep=""),
  paste(path4, "/polyA_signal.txt", sep=""),
  paste(path5, "/polyA_signal.txt", sep=""),
  paste(path6, "/polyA_signal.txt", sep="")
)


ntf_di = c(
  paste(path1, "/nt_di_freq.txt", sep=""),
  paste(path2, "/nt_di_freq.txt", sep=""),
  paste(path3, "/nt_di_freq.txt", sep=""),
  paste(path4, "/nt_di_freq.txt", sep=""),
  paste(path5, "/nt_di_freq.txt", sep=""),
  paste(path6, "/nt_di_freq.txt", sep="")
)

exp = c(
  paste(path1, "/all_tissue_exp_feat_colwise.txt", sep=""),
  paste(path2, "/all_tissue_exp_feat_colwise.txt", sep=""),
  paste(path3, "/all_tissue_exp_feat_colwise.txt", sep=""),
  paste(path4, "/all_tissue_exp_feat_colwise.txt", sep=""),
  paste(path5, "/all_tissue_exp_feat_colwise.txt", sep=""),
  paste(path6, "/all_tissue_exp_feat_colwise.txt", sep="")
)

phastCons = c(
  paste(path1, "/phastCons_feat.txt", sep=""),
  paste(path2, "/phastCons_feat.txt", sep=""),
  paste(path3, "/phastCons_feat.txt", sep=""),
  paste(path4, "/phastCons_feat.txt", sep=""),
  paste(path5, "/phastCons_feat.txt", sep=""),
  paste(path6, "/phastCons_feat.txt", sep="")
)

# dna structural properties
dsp = c(
  paste(path1, "/structural_feat_mean.txt", sep=""),
  paste(path2, "/structural_feat_mean.txt", sep=""),
  paste(path3, "/structural_feat_mean.txt", sep=""),
  paste(path4, "/structural_feat_mean.txt", sep=""),
  paste(path5, "/structural_feat_mean.txt", sep=""),
  paste(path6, "/structural_feat_mean.txt", sep="")
)

# repeats
repeats = c(
  paste(path1, "/repeats.txt", sep=""),
  paste(path2, "/repeats.txt", sep=""),
  paste(path3, "/repeats.txt", sep=""),
  paste(path4, "/repeats.txt", sep=""),
  paste(path5, "/repeats.txt", sep=""),
  paste(path6, "/repeats.txt", sep="")
)



df1 = pa_signal %>%
  map(read.table, header=TRUE) %>%
  bind_rows() %>%
  distinct(id, .keep_all = TRUE) %>%
  dplyr::select(c(id, polyA_signal)) %>%
  mutate(polyA_signal = as.factor(polyA_signal))


df2 = exp %>%
  map(read.table, header=TRUE) %>%
  bind_rows() %>%
  distinct(id, .keep_all = TRUE) %>%
  #dplyr::select(-c(id), -contains("meanCov_")) %>%
  transmute(meanPd = rowMeans(dplyr::select(., contains("pd_"))),
            meanEE = rowMeans(dplyr::select(., contains("EE_")))
  )


df3 = phastCons %>%
  map(read.table, header=TRUE) %>%
  bind_rows() %>%
  distinct(id, .keep_all = TRUE) %>%
  dplyr::select(-c(id, class))


df4 = dsp %>%
  map(read.table, header=TRUE) %>%
  bind_rows() %>%
  distinct(id, .keep_all = TRUE) %>%
  dplyr::select(-c(id, class))

df5 = repeats %>%
  map(read.table, header=TRUE) %>%
  bind_rows() %>%
  distinct(id, .keep_all = TRUE) %>%
  plyr::rename(c("fracOverlap" = "RepeatsOverlap")) %>%
  dplyr::select(-c(id, class))

df6 = ntf_di %>%
  map(read.table, header=TRUE) %>%
  bind_rows() %>%
  distinct(id, .keep_all = TRUE) %>%
  dplyr::select(-c(id))


data = cbind(df1, df2, df3, df4, df5, df6) %>%
  mutate(class = as.factor(ifelse(class %in% "3-UTR", "UTR", class))) %>%
  tidyr::drop_na() %>%
  remove_rownames() %>%
  column_to_rownames("id")



file = paste(out.path, "/ml_table.rds", sep="")
saveRDS(data, file= file)
