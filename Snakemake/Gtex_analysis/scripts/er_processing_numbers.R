args<-commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
})


stats_dir <- args[1]
outpath <- args[2]

files <- list.files(stats_dir, full.names = TRUE)

data = files %>%
  map_dfr(read.table, header=TRUE)


all_data <- data %>% dplyr::filter(type %in% "all_data")
all_data_perc <- sum(all_data$is_associated_gene_also_nearest)/sum(all_data$total)*100


all_proximity <- data %>% dplyr::filter(type %in% "all_proximity_data")
all_proximity_perc <- sum(all_proximity$is_associated_gene_also_nearest)/sum(all_proximity$total)*100


three_prime <- data %>% dplyr::filter(type %in% "3prime_proximity_data")
three_prime_perc <- sum(three_prime$is_associated_gene_also_nearest)/sum(three_prime$total)*100


cat("percentage of associated genes equal nearest gene in complete dataset", all_data_perc, "\n",
    "percentage of associated genes equal nearest gene in prximity dataset", all_proximity_perc, "\n",
    "percentage of associated genes equal nearest gene in 3 prime proximity dataset", three_prime_perc, "\n",
    file = str_c(outpath, "/numbers.txt")
    )
