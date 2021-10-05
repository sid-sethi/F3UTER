args<-commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
})


gtex_file <- args[1]
outpath <- args[2]

data <- read.table(gtex_file, header=TRUE, sep="\t")

data <- data %>% mutate(Name = str_split_fixed(Name, "\\.", n=2)[,1])

names <- colnames(data)
names <- names %>%
  str_replace_all("\\.+", "_") %>%
  str_replace("_$", "") %>%
  tolower()

for(i in 3:ncol(data)){

  file_name <- names[i]
  print(str_c(Sys.time(), "-", i, "-", file_name))
  df <- data[,c(1,i)] %>% dplyr::filter(.[[2]] > 0.1)

  write.table(df, str_c(outpath, "/", file_name, ".txt"), quote=FALSE, row.names=FALSE, sep="\t")

  rm(df)
}
