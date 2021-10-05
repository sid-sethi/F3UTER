args<-commandArgs(TRUE)

suppressPackageStartupMessages({
  library(rtracklayer)
  library(tidyverse)
})

options(stringsAsFactors=F)

merged_data <- args[1]
gtf_file <- args[2]
outpath <- args[3]

gtf <- rtracklayer::import(gtf_file)
gtf = keepStandardChromosomes(gtf, species = "Homo_sapiens", pruning.mode="coarse")
gtf.df=as.data.frame(gtf,stringsAsFactor=F)
gtf.gene <- gtf.df %>% dplyr::filter(type %in% "gene")
gtf_gene_gr = makeGRangesFromDataFrame(gtf.gene, keep.extra.columns = TRUE)


all_tissues <- readRDS(merged_data)

data <- all_tissues %>%
  mutate(predicted.response = ifelse(predicted.prob > 0.6, "Predicted 3-UTRs", "Predicted non-3-UTRs"),
         gene_association = ifelse(annotationType_split_read_annot %in% "partially annotated split read", "split-read", "proximity"))



all_genes <- data %>%
  dplyr::distinct(associated_gene)

pos_genes <- data %>% dplyr::filter(predicted.response %in% "Predicted 3-UTRs") %>% dplyr::distinct(associated_gene)

all_proximity_genes <- data %>%
  dplyr::filter(gene_association %in% "proximity") %>%
  dplyr::distinct(associated_gene)


genes.data <- dplyr::left_join(all_proximity_genes, gtf.gene, by = c("associated_gene" = "gene_id"))
genes_data_gr = makeGRangesFromDataFrame(genes.data, keep.extra.columns = TRUE)


x = IRanges::findOverlapPairs(genes_data_gr, gtf_gene_gr, ignore.strand = TRUE) %>%
  as.data.frame() %>%
  mutate(overlap_other_gene = ifelse(first.X.associated_gene == second.X.gene_id, "NO", "YES"),
         overlap_strand_same = ifelse(first.X.strand == second.X.strand, "YES", "NO"))


overlap_yes = x %>%
  dplyr::filter(overlap_other_gene == "YES") %>%
  group_by(first.X.associated_gene) %>%
  summarise(strand = str_c(unique(overlap_strand_same), collapse = ":"),
            biotype = str_c(unique(second.X.gene_biotype), collapse = ":"), .groups = "keep")



genes_affected <- overlap_yes %>% dplyr::filter(str_detect(biotype, "protein_coding"), str_detect(strand, "YES"))
data$overlap_other_gene <- ifelse(data$associated_gene %in% genes_affected$first.X.associated_gene, "YES", "NO")
ERs_affected <- data %>% dplyr::filter(gene_association %in% "proximity", overlap_other_gene %in% "YES") %>% .$associated_gene %>% length()


pos_genes_affected <- data %>% dplyr::filter(predicted.response %in% "Predicted 3-UTRs", gene_association %in% "proximity", overlap_other_gene %in% "YES") %>% dplyr::distinct(associated_gene)
pos_preds <- data %>% dplyr::filter(predicted.response %in% "Predicted 3-UTRs")
pos_ers_affected <- pos_preds %>% dplyr::filter(gene_association %in% "proximity", overlap_other_gene %in% "YES") %>% .$associated_gene %>% length()

cat(
  "Total number of ERs input to F3UTER - ", nrow(data), "\n",
  "Total number of distinct genes associated with all ER predictions - ", nrow(all_genes), "\n",
  "Total number of ERs predicted as 3'UTRs - ", nrow(pos_preds), "\n",
  "Number of distinct genes associated with predited 3'UTR predictions - ", nrow(pos_genes), "\n",
  "Overlap_other_gene = overlap with other protein-coding genes on the same strand", "\n",
  "Number of distinct proximity genes with overlap_other_gene - ", nrow(genes_affected), "\n",
  "% of total genes affected by overlap_other_gene - ", (nrow(genes_affected)/nrow(all_genes))*100, "\n",
  "% of total ER predictions affected by overlapping genes - ", (ERs_affected/nrow(data))*100, "\n",
  "Number of distinct proximity genes with overlap_other_gene in predicted 3'UTR predictions- ", nrow(pos_genes_affected), "\n",
  "% of predited 3'UTR-gene associations affected by overlapping genes - ", (pos_ers_affected/nrow(pos_preds))*100, "\n",
  file = str_c(outpath, "/stats_overlapping_genes_effect.txt")
)
