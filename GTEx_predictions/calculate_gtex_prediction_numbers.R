args <- commandArgs(TRUE) # 1st  = ML output folder

library(tidyverse)
library(RColorBrewer)
options(stringsAsFactors=F)

mldir = args[1]

path = str_c(mldir, "/PredictionMerge")

files <- list.files(path = path, pattern = ".rds")

# new output path
out.path = paste(mldir, "/PredictionMerge/", sep="")

print(str_c(Sys.time(), " - merging data from all tissues"))

res60 = data.frame()
res = data.frame()
distData = data.frame()
all_genes = data.frame()

i = 0
for(f in files){
  i= i + 1
  data <- readRDS(str_c(path, "/", f))

  pos = data %>% dplyr::filter(predicted.prob > 0.60) %>% nrow()
  neg = data %>% dplyr::filter(predicted.prob <= 0.60) %>% nrow()
  tis <- str_replace(f, "_mergedData.rds", "") %>% str_replace(., "brain_", "")

  pos.length = data %>% dplyr::filter(predicted.prob > 0.60) %>% .$width %>% sum()

  res60[i, "tissue"] <- tis
  res60[i, "UTR"] <- pos
  res60[i, "Non-3-UTR"] <- neg
  res60[i, "total"] <- nrow(data)
  res60[i, "percent_utr"] <- round((pos/nrow(data))*100,2)
  res60[i, "utr_total_length"] <- pos.length
  res60[i, "utr_gene_freq"] <- data %>% dplyr::filter(predicted.prob > 0.60) %>% distinct(.$associated_gene) %>% nrow()


  geneNames <- data %>% dplyr::filter(predicted.prob > 0.60) %>% distinct(associated_gene) %>%
    mutate(tissue = tis)
  all_genes <- rbind(all_genes, geneNames)

  # all prob
  breaks = seq(0,1,0.10)

  df = data %>%
    group_by(range=cut(predicted.prob, breaks= seq(0, 1, by = 0.10)) ) %>%
    summarise(frequency= n(), length = sum(width), Gene_freq= n_distinct(associated_gene)) %>%
    arrange(as.numeric(range)) %>%
    mutate(threshold = breaks[-1], tissue = tis) %>%
    as.data.frame()

  res = rbind(res, df)

  dist = data %>% dplyr::select(predicted.prob, associated_gene, distance_from_geneEnd, tissue)
  distData = rbind(distData, dist)

}

# % of positive predictions
mean_all_tissues = res60$percent_utr %>% mean()
print(str_c("average % of predictions across all tissues = ", mean_all_tissues))
res60$UTR %>% mean()
range(res60$UTR)
res60$utr_total_length %>% mean()/1000
res60$utr_total_length %>% sum()/1000000 ## total length of predictions in Mb
range(res60$utr_total_length)/1000
res60$utr_gene_freq %>% mean()
range(res60$utr_gene_freq)

brain_tissues <- c("amygdala", "anterior_cingulate_cortex_ba24", "caudate_basal_ganglia", "cerebellar_hemisphere", "cerebellum", "cortex", "frontal_cortex_ba9", "hippocampus", "hypothalamus", "nucleus_accumbens_basal_ganglia", "putamen_basal_ganglia", "spinal_cord_cervical_c_1", "substantia_nigra")

mean_brain_tissues <- res60 %>% dplyr::filter(tissue %in% brain_tissues) %>% .$percent_utr %>% mean()
print(str_c("average % of predictions across brain tissues = ", mean_brain_tissues))
res60.brain <- res60 %>% dplyr::filter(tissue %in% brain_tissues)
res60.nbrain <- res60 %>% dplyr::filter(!tissue %in% brain_tissues)
man.pval.number = wilcox.test(res60.brain$UTR,res60.nbrain$UTR)$p.value %>% format.pval(3)
man.pval.length = wilcox.test(res60.brain$utr_total_length,res60.nbrain$utr_total_length)$p.value %>% format.pval(3)
man.pval.gene_freq = wilcox.test(res60.brain$utr_gene_freq,res60.nbrain$utr_gene_freq)$p.value %>% format.pval(3)
median(res60.brain$UTR)
median(res60.nbrain$UTR)
median(res60.brain$utr_total_length)/1000
median(res60.nbrain$utr_total_length)/1000
median(res60.brain$utr_gene_freq)
median(res60.nbrain$utr_gene_freq)


# total distinct genes with unannotated 3'UTR
all_genes$associated_gene %>% unique() %>%  length()
