args <- commandArgs(TRUE) # 1st  = ML output folder

library(optparse)
option_list = list(
  make_option(c("-d", "--dir"), type="character", default=NULL,
              help="make predictions output directory", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

dir <- opt$dir


library(tidyverse)
library(rtracklayer)
library(ggpubr)
library(pheatmap)
library(gprofiler2)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
options(stringsAsFactors = FALSE)


######################################################################
## functions #####

# make_directory function
make_dir <- function(resPath, name){
  if(!dir.exists(paste(resPath,"/", name, sep=""))){
    system(paste("mkdir -m a=rwx ",resPath, "/", name, sep=""))
  }
}

# running g-profiler
run_go <- function(x){ # x is a vector of gene ids
  gostres <- gost(query = x,
                  organism = "hsapiens", ordered_query = FALSE,
                  multi_query = FALSE, significant = TRUE,
                  measure_underrepresentation = FALSE, evcodes = FALSE,
                  user_threshold = 0.05, correction_method = "fdr",
                  domain_scope = "custom",
                  custom_bg = go_bg_pc$gene_id,
                  sources = c("GO:MF", "GO:BP", "WP", "MIRNA", "HP")
  )

  if(!is.null(gostres)){
    go <- gostres$result %>% dplyr::select(-parents)
    return(go)
  }
}


go_table <- function(df, out_file){

  total = nrow(df)

  if(is.null(total)){
    return(NULL)
  } else if(total > 40){
    rows_to_print <- 40
  } else{
    rows_to_print <- total
  }

  pt <- publish_gosttable(df, highlight_terms = df[c(1:rows_to_print),],
                          use_colors = TRUE,
                          show_columns = c("source", "term_name", "term_size", "query_size", "intersection_size"),
                          filename = out_file)
}



# extracting gene lists
split_genes <- function(df) {
  all = df %>%
    dplyr::left_join(geneData, by = c("associated_gene" = "gene_id")) %>%
    dplyr::select(id:predicted.prob, nearest_any_gene_v92_distance:nearest_any_gene_v92_name, tissue, annotationType_split_read_annot, uniq_genes_split_read_annot, associated_gene, distance_from_geneTSS, distance_from_geneEnd, gene_name) %>%
    plyr::rename(c("gene_name" = "associated_gene_name"))


  pos = df %>% dplyr::filter(predicted.prob > 0.60) %>%
    dplyr::left_join(geneData, by = c("associated_gene" = "gene_id")) %>%
    dplyr::select(id:predicted.prob, nearest_any_gene_v92_distance:nearest_any_gene_v92_name, tissue, annotationType_split_read_annot, uniq_genes_split_read_annot, associated_gene, distance_from_geneTSS, distance_from_geneEnd, gene_name) %>%
    plyr::rename(c("gene_name" = "associated_gene_name"))


  neg = df %>% dplyr::filter(predicted.prob <= 0.60) %>%
    dplyr::left_join(geneData, by = c("associated_gene" = "gene_id")) %>%
    dplyr::select(id:predicted.prob, nearest_any_gene_v92_distance:nearest_any_gene_v92_name, tissue, annotationType_split_read_annot, uniq_genes_split_read_annot, associated_gene, distance_from_geneTSS, distance_from_geneEnd, gene_name) %>%
    plyr::rename(c("gene_name" = "associated_gene_name"))

  pos_minimal = df %>% dplyr::filter(predicted.prob > 0.60) %>%
    dplyr::select(seqnames:end, id, predicted.prob, tissue)

  neg_minimal = df %>% dplyr::filter(predicted.prob <= 0.60) %>%
    dplyr::select(seqnames:end, id, predicted.prob, tissue)

  res <- list(all = all, pos = pos, neg = neg, pos_minimal = pos_minimal, neg_minimal = neg_minimal)

  return(res)
}

# for shared ER data, different columns
split_genes2 <- function(df) {
  all = df %>%
    dplyr::left_join(geneData, by = c("associated_gene" = "gene_id")) %>%
    dplyr::select(id:distance_from_geneEnd, gene_name) %>%
    plyr::rename(c("gene_name" = "associated_gene_name"))


  pos = df %>% dplyr::filter(predicted.prob > 0.60) %>%
    dplyr::left_join(geneData, by = c("associated_gene" = "gene_id")) %>%
    dplyr::select(id:distance_from_geneEnd, gene_name) %>%
    plyr::rename(c("gene_name" = "associated_gene_name"))


  neg = df %>% dplyr::filter(predicted.prob <= 0.60) %>%
    dplyr::left_join(geneData, by = c("associated_gene" = "gene_id")) %>%
    dplyr::select(id:distance_from_geneEnd, gene_name) %>%
    plyr::rename(c("gene_name" = "associated_gene_name"))

  pos_minimal = df %>% dplyr::filter(predicted.prob > 0.60) %>%
    separate(id, c("tissue", "seqnames", "start", "end", "strand"), sep = ":", remove =FALSE) %>%
    dplyr::select(seqnames:end, id, predicted.prob, tissue)

  neg_minimal = df %>% dplyr::filter(predicted.prob <= 0.60) %>%
    separate(id, c("tissue", "seqnames", "start", "end", "strand"), sep = ":", remove =FALSE) %>%
    dplyr::select(seqnames:end, id, predicted.prob, tissue)


  res <- list(all = all, pos = pos, neg = neg, pos_minimal = pos_minimal, neg_minimal = neg_minimal)

  return(res)
}


# disease enrichment using clusterprofile
run_do <- function(x){
  entrez <- bitr(x, fromType="ENSEMBL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
  dores <- enrichDO(gene          = entrez$ENTREZID,
                ont           = "DO",
                pvalueCutoff  = 0.05,
                pAdjustMethod = "BH",
                #universe      = names(geneList),
                minGSSize     = 5,
                maxGSSize     = 500,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
  do <- dores@result %>% dplyr::filter(pvalue < 0.05)
  return(do)
}


# DisGenet enrichment using clusterprofile
run_dsgn <- function(x){
  entrez <- bitr(x, fromType="ENSEMBL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
  dsgnres <- enrichDGN(entrez$ENTREZID,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   minGSSize = 10, maxGSSize = 500,
                   qvalueCutoff = 0.2,
                   readable = TRUE)
  dsgn <- dsgnres@result %>% dplyr::filter(pvalue < 0.05)
  return(dsgn)
}

########################################################################
########################################################################
########################################################################
########################################################################



# new output path
make_dir(str_c(dir, "/PredictionMerge/Tis_spec"), "GeneList")
make_dir(str_c(dir, "/PredictionMerge/Brain_spec"), "GeneList")
make_dir(str_c(dir, "/PredictionMerge/Common"), "GeneList")
make_dir(str_c(dir, "/PredictionMerge/Ambiguous"), "GeneList")

out.path1 = str_c(dir, "/PredictionMerge/Tis_spec/GeneList")
out.path2 = str_c(dir, "/PredictionMerge/Brain_spec/GeneList")
out.path3 = str_c(dir, "/PredictionMerge/Common/GeneList")
out.path4 = str_c(dir, "/PredictionMerge/Ambiguous/GeneList")


tissues = c(
  "adipose_subcutaneous",
  "adrenal_gland",
  "artery_aorta",
  "artery_coronary",
  "artery_tibial",
  "adipose_visceral_omentum",
  "brain_amygdala",
  "brain_anterior_cingulate_cortex_ba24",
  "brain_caudate_basal_ganglia",
  "brain_cerebellar_hemisphere",
  "brain_frontal_cortex_ba9",
  "brain_hippocampus",
  "brain_hypothalamus",
  "brain_nucleus_accumbens_basal_ganglia",
  "brain_putamen_basal_ganglia",
  "brain_spinal_cord_cervical_c_1",
  "brain_substantia_nigra",
  "colon_sigmoid",
  "colon_transverse",
  "esophagus_gastroesophageal_junction",
  "esophagus_mucosa",
  "esophagus_muscularis",
  "heart_atrial_appendage",
  "heart_left_ventricle",
  "kidney_cortex",
  "liver",
  "lung",
  "minor_salivary_gland",
  "muscle_skeletal",
  "nerve_tibial",
  "pancreas",
  "pituitary",
  "skin_not_sun_exposed_suprapubic",
  "skin_sun_exposed_lower_leg",
  "small_intestine_terminal_ileum",
  "spleen",
  "stomach",
  "thyroid",
  "whole_blood"
)

# reading data containing gene name and ids from Ensembl
geneData <- readRDS("/Annotation/geneData_v92.rds")
go_bg_pc <- geneData %>% dplyr::filter(gene_biotype %in% "protein_coding")

tis_spec_list <- list()
brain_spec_list <- list()
common_list <- list()
ambi_list <- list()

for(tissue in tissues){

  print(str_c(Sys.time(), " - ", tissue))

  # tis_spec
  data <- readRDS(str_c(dir, "/PredictionMerge/Tis_spec/", tissue, "_tisSpec.rds"))

  data.genes <- split_genes(data)

  tis_spec_list$all <- rbind(tis_spec_list$all, data.genes$all)
  tis_spec_list$pos <- rbind(tis_spec_list$pos, data.genes$pos)
  tis_spec_list$neg <- rbind(tis_spec_list$neg, data.genes$neg)
  tis_spec_list$pos_minimal <- rbind(tis_spec_list$pos_minimal, data.genes$pos_minimal)
  tis_spec_list$neg_minimal <- rbind(tis_spec_list$neg_minimal, data.genes$neg_minimal)

  rm(data, data.genes)


  # brain_spec
  data <- readRDS(str_c(dir, "/PredictionMerge/Brain_spec/", tissue, "_brainSpec.rds"))

  data.genes <- split_genes(data)

  brain_spec_list$all <- rbind(brain_spec_list$all, data.genes$all)
  brain_spec_list$pos <- rbind(brain_spec_list$pos, data.genes$pos)
  brain_spec_list$neg <- rbind(brain_spec_list$neg, data.genes$neg)
  brain_spec_list$pos_minimal <- rbind(brain_spec_list$pos_minimal, data.genes$pos_minimal)
  brain_spec_list$neg_minimal <- rbind(brain_spec_list$neg_minimal, data.genes$neg_minimal)

  rm(data, data.genes)


  # shared
  data <- readRDS(str_c(dir, "/PredictionMerge/Common/", tissue, "_commonId.rds"))
  colnames(data) <- str_replace(colnames(data), "first.X.", "")

  data.genes <- split_genes2(data)

  common_list$all <- rbind(common_list$all, data.genes$all)
  common_list$pos <- rbind(common_list$pos, data.genes$pos)
  common_list$neg <- rbind(common_list$neg, data.genes$neg)
  common_list$pos_minimal <- rbind(common_list$pos_minimal, data.genes$pos_minimal)
  common_list$neg_minimal <- rbind(common_list$neg_minimal, data.genes$neg_minimal)

  rm(data, data.genes)


  # ambiguous
  data <- readRDS(str_c(dir, "/PredictionMerge/Ambiguous/", tissue, "_ambigs.rds"))

  data.genes <- split_genes(data)

  ambi_list$all <- rbind(ambi_list$all, data.genes$all)
  ambi_list$pos <- rbind(ambi_list$pos, data.genes$pos)
  ambi_list$neg <- rbind(ambi_list$neg, data.genes$neg)
  ambi_list$pos_minimal <- rbind(ambi_list$pos_minimal, data.genes$pos_minimal)
  ambi_list$neg_minimal <- rbind(ambi_list$neg_minimal, data.genes$neg_minimal)

  rm(data, data.genes)


}



# Saving gene lists

##  Tissue-specific
write.table(tis_spec_list$all, file = str_c(out.path1,"/tis_spec_all_data.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(tis_spec_list$pos, file = str_c(out.path1,"/tis_spec_pos_data.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(tis_spec_list$neg, file = str_c(out.path1,"/tis_spec_neg_data.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(tis_spec_list$pos_minimal, file = str_c(out.path1,"/tis_spec_pos_data_minimal.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(tis_spec_list$neg_minimal, file = str_c(out.path1,"/tis_spec_neg_data_minimal.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

tis_spec_uniq_genes_pos <- tis_spec_list$pos %>%
  distinct(associated_gene, .keep_all = TRUE) %>%
  dplyr::select(associated_gene, associated_gene_name)

tis_spec_uniq_genes_pos_2 <- tis_spec_list$pos %>%
  group_by(associated_gene, associated_gene_name) %>%
  summarise(tis_spec_tissues = str_c(unique(tissue), collapse = ":")) %>%
  as.data.frame()

tis_spec_uniq_genes_neg <- tis_spec_list$neg %>%
  distinct(associated_gene, .keep_all = TRUE) %>%
  dplyr::select(associated_gene, associated_gene_name)

write.table(tis_spec_uniq_genes_pos, file = str_c(out.path1,"/tis_spec_pos_genes.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(tis_spec_uniq_genes_pos_2, file = str_c(out.path1,"/tis_spec_pos_genes_2.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(tis_spec_uniq_genes_neg, file = str_c(out.path1,"/tis_spec_neg_genes.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

## Brain-specific
write.table(brain_spec_list$all, file = str_c(out.path2,"/brain_spec_all_data.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(brain_spec_list$pos, file = str_c(out.path2,"/brain_spec_pos_data.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(brain_spec_list$neg, file = str_c(out.path2,"/brain_spec_neg_data.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(brain_spec_list$pos_minimal, file = str_c(out.path2,"/brain_spec_pos_data_minimal.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(brain_spec_list$neg_minimal, file = str_c(out.path2,"/brain_spec_neg_data_minimal.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

brain_spec_uniq_genes_pos <- brain_spec_list$pos %>%
  distinct(associated_gene, .keep_all = TRUE) %>%
  dplyr::select(associated_gene, associated_gene_name)

brain_spec_uniq_genes_neg <- brain_spec_list$neg %>%
  distinct(associated_gene, .keep_all = TRUE) %>%
  dplyr::select(associated_gene, associated_gene_name)

write.table(brain_spec_uniq_genes_pos, file = str_c(out.path2,"/brain_spec_pos_genes.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(brain_spec_uniq_genes_neg, file = str_c(out.path2,"/brain_spec_neg_genes.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)


##  sHARED
write.table(common_list$all, file = str_c(out.path3,"/common_all_data.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(common_list$pos, file = str_c(out.path3,"/common_pos_data.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(common_list$neg, file = str_c(out.path3,"/common_neg_data.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(common_list$pos_minimal, file = str_c(out.path3,"/common_pos_data_minimal.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(common_list$neg_minimal, file = str_c(out.path3,"/common_neg_data_minimal.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

common_uniq_genes_pos <- common_list$pos %>%
  distinct(associated_gene, .keep_all = TRUE) %>%
  dplyr::select(associated_gene, associated_gene_name)

common_uniq_genes_neg <- common_list$neg %>%
  distinct(associated_gene, .keep_all = TRUE) %>%
  dplyr::select(associated_gene, associated_gene_name)

write.table(common_uniq_genes_pos, file = str_c(out.path3,"/common_pos_genes.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(common_uniq_genes_neg, file = str_c(out.path3,"/common_neg_genes.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)


## Ambiguous

write.table(ambi_list$all, file = str_c(out.path4,"/ambi_all_data.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(ambi_list$pos, file = str_c(out.path4,"/ambi_pos_data.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(ambi_list$neg, file = str_c(out.path4,"/ambi_neg_data.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(ambi_list$pos_minimal, file = str_c(out.path4,"/ambi_pos_data_minimal.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(ambi_list$neg_minimal, file = str_c(out.path4,"/ambi_neg_data_minimal.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

ambi_uniq_genes_pos <- ambi_list$pos %>%
  distinct(associated_gene, .keep_all = TRUE) %>%
  dplyr::select(associated_gene, associated_gene_name)

ambi_uniq_genes_neg <- ambi_list$neg %>%
  distinct(associated_gene, .keep_all = TRUE) %>%
  dplyr::select(associated_gene, associated_gene_name)

write.table(ambi_uniq_genes_pos, file = str_c(out.path4,"/ambi_pos_genes.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(ambi_uniq_genes_neg, file = str_c(out.path4,"/ambi_neg_genes.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
