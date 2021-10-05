args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
})

options(stringsAsFactors = FALSE, warn = -1)

rds = args[1]
gtf_file <- args[2]
outpath <- args[3]
prefix <- args[4]


gtf <- rtracklayer::import(gtf_file)
gtf = keepStandardChromosomes(gtf, species = "Homo_sapiens", pruning.mode="coarse")
gtf.df = as.data.frame(gtf,stringsAsFactor=F)
gtf.gene = gtf.df %>% dplyr::filter(type %in% "gene")

data <- readRDS(rds)

pos = data %>% dplyr::filter(predicted.prob > 0.60) %>%
  dplyr::left_join(gtf.gene, by = c("associated_gene" = "gene_id")) %>%
  distinct(associated_gene, .keep_all = TRUE) %>%
  dplyr::select(associated_gene, gene_name) %>%
  plyr::rename(c("gene_name" = "associated_gene_name"))


neg = data %>% dplyr::filter(predicted.prob <= 0.60) %>%
  dplyr::left_join(gtf.gene, by = c("associated_gene" = "gene_id")) %>%
  distinct(associated_gene, .keep_all = TRUE) %>%
  dplyr::select(associated_gene, gene_name) %>%
  plyr::rename(c("gene_name" = "associated_gene_name"))


write.table(pos, file = str_c(outpath, "/", prefix, "_pos_genes.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(neg, file = str_c(outpath, "/", prefix, "_neg_genes.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
