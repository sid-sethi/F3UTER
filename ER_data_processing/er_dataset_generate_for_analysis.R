library(tidyverse)
library(stringr)
library(rtracklayer)
library(ggpubr)

options(stringsAsFactors=F)

tissues = c(
  "adipose_subcutaneous",
  "adrenal_gland",
  "artery_aorta",
  "artery_coronary",
  "artery_tibial",
  "adipose_visceral_omentum",
  ############## "bladder",
  "brain_amygdala",
  "brain_anterior_cingulate_cortex_ba24",
  "brain_caudate_basal_ganglia",
  "brain_cerebellar_hemisphere",
  #"brain_cerebellum",
  #"brain_cortex",
  "brain_frontal_cortex_ba9",
  "brain_hippocampus",
  "brain_hypothalamus",
  "brain_nucleus_accumbens_basal_ganglia",
  "brain_putamen_basal_ganglia",
  "brain_spinal_cord_cervical_c_1",
  "brain_substantia_nigra",
  ############## "cells_ebv_transformed_lymphocytes",
  ############## "cells_transformed_fibroblasts",
  ############## "cervix_ectocervix",
  "colon_sigmoid",
  "colon_transverse",
  "esophagus_gastroesophageal_junction",
  "esophagus_mucosa",
  "esophagus_muscularis",
  ############## "fallopian_tube",
  "heart_atrial_appendage",
  "heart_left_ventricle",
  "kidney_cortex",
  "liver",
  "lung",
  "minor_salivary_gland",
  "muscle_skeletal",
  "nerve_tibial",
  ############## "ovary",
  "pancreas",
  "pituitary",
  ############## "prostate",
  "skin_not_sun_exposed_suprapubic",
  "skin_sun_exposed_lower_leg",
  "small_intestine_terminal_ileum",
  "spleen",
  "stomach",
  ############## "testis",
  "thyroid",
  ############## "uterus",
  ############## "vagina",
  "whole_blood"
)


####################################################################################################################
###################################################      Functions      ############################################
####################################################################################################################

#### Adding 5' or 3' tag to ERs (David Zhang's function - modified by Sid)
get_5_3_inter <- function(er_input, ensembl_grch38_v92_gtf_gr){

  # get only intergenic regions connected to ONE gene through split read
  er_inter_sr <- er_input %>%
    filter(ensembl_grch38_v92_region_annot %in% c("intergenic"), !is.na(uniq_genes_split_read_annot), !str_detect(uniq_genes_split_read_annot, ",")) %>%
    mutate(associated_gene = uniq_genes_split_read_annot)

  # exonic intergenic
  er_exon_inter <- er_input %>%
      filter(ensembl_grch38_v92_region_annot %in% c("exon, intergenic"), !str_detect(overlap_any_gene_v92_name, ",")) %>%
      mutate(associated_gene = overlap_any_gene_v92_name)

   # intergenic with no split reads
  er_inter_nosr <- er_input %>%
      filter(ensembl_grch38_v92_region_annot %in% c("intergenic"), annotationType_split_read_annot != "partially annotated split read") %>%
      mutate(associated_gene = nearest_any_gene_v92_name)

  er_combined <-
      bind_rows(er_inter_sr, er_exon_inter, er_inter_nosr)

  # subset gtf to only genes we have a inter ER connected to through split read
  er_gene_gtf <-
    ensembl_grch38_v92_gtf_gr[ensembl_grch38_v92_gtf_gr$type == "gene" &
                                ensembl_grch38_v92_gtf_gr$gene_id %in% unique(er_combined$associated_gene)] %>%
    as.data.frame() %>%
    dplyr::select(seqnames, start, end, gene_id, strand) # select only necessary columns for clarity rather than necessity

  stopifnot(all(unique(er_combined$associated_gene) %in% er_gene_gtf$gene_id))

  # join the inter ER details with the gene details
  res <- er_combined %>%
    left_join(er_gene_gtf, by = c("associated_gene" = "gene_id")) %>%
    mutate(diff_ends = end.y - end.x,
           prime_5_3 = ifelse(diff_ends <= 0, "downstream", "upstream"),
           prime_5_3_corrected = ifelse(strand.y == "-",
                                        ifelse(prime_5_3 == "upstream", 3, 5),
                                        ifelse(prime_5_3 == "upstream", 5, 3)),
           prime_5_3_corrected_chr = ifelse(prime_5_3_corrected == 3, "3'", "5'"),
           distance_from_geneTSS = ifelse(strand.y == "+",
                                          end.x - start.y,
                                          start.x - end.y
                                          ),
           distance_from_geneEnd = ifelse(strand.y == "+",
                                          start.x - end.y,
                                          end.x - start.y
                                          )
           )

  return(res)

}

ensembl_grch38_v92_gtf_gr <- import("/Annotation/Homo_sapiens.GRCh38.92.gtf")


#make_directory function
make_dir <- function(resPath, name){
  if(!dir.exists(paste(resPath,"/", name, sep=""))){
    system(paste("mkdir -m a=rwx ",resPath, "/", name, sep=""))
  }
}


####################################################################################################################
####################################################################################################################
####################################################################################################################


erPath = "/ERs/ERs_by_tissue_processed"
main_out_path = "/ERs/Validation_datasets"

make_dir(main_out_path, "Intergenic")
make_dir(str_c(main_out_path,"/Intergenic"), "5_prime")
make_dir(str_c(main_out_path,"/Intergenic"), "3_prime")


for(tissue in tissues){

  print(str_c(Sys.time(), " - ", tissue))

  #ER data and processing
  er = read.table(paste(erPath, "/", tissue, "_processed.txt", sep=""), header=T, sep="\t")

  ########################################
  ############# INTERGENIC  ##############
  ########################################


  # ERs with split read
  inter.sr = er %>% filter(ensembl_grch38_v92_region_annot %in% "intergenic",
                           annotationType_split_read_annot %in% "partially annotated split read") %>%
    filter(!str_detect(uniq_genes_split_read_annot, ",")) %>%
    get_5_3_inter(ensembl_grch38_v92_gtf_gr) %>%
    filter(uniq_genes_split_read_annot_biotype %in% "protein_coding") %>%
    plyr::rename(c("seqnames.x"="seqnames", "start.x"="start", "end.x"="end", "strand.x"="strand")) %>%
    mutate(strand = uniq_genes_split_read_annot_strand) %>%
    dplyr::select(-c(seqnames.y, start.y, end.y, strand.y, diff_ends, prime_5_3, prime_5_3_corrected_chr))

  inter.sr.5 = inter.sr %>% filter(prime_5_3_corrected %in% 5,
                                  distance_from_geneTSS < 10000, distance_from_geneTSS > -10000)
  inter.sr.3 = inter.sr %>% filter(prime_5_3_corrected %in% 3,
                                  distance_from_geneEnd < 10000, distance_from_geneEnd > -10000)


  # ERs without split read
  inter.nosr = er %>% filter(ensembl_grch38_v92_region_annot %in% "intergenic",
                             annotationType_split_read_annot != "partially annotated split read") %>%
    get_5_3_inter(ensembl_grch38_v92_gtf_gr) %>%
    plyr::rename(c("seqnames.x"="seqnames", "start.x"="start", "end.x"="end", "strand.x"="strand")) %>%
    dplyr::select(-c(seqnames.y, start.y, end.y, strand.y, diff_ends, prime_5_3, prime_5_3_corrected_chr))


  inter.nosr.5 = inter.nosr %>% filter(prime_5_3_corrected %in% 5,
                                       nearest_any_gene_v92_distance < 10000,
                                       nearest_any_gene_v92_name_biotype %in% "protein_coding")

  inter.nosr.3 = inter.nosr %>% filter(prime_5_3_corrected %in% 3,
                                       nearest_any_gene_v92_distance < 10000,  #1000000
                                       nearest_any_gene_v92_name_biotype %in% "protein_coding")



  # combining both categories to get filtered ERs
  intergenic.5 = rbind(inter.sr.5, inter.nosr.5) %>%
    mutate(ER_id = paste(tissue, seqnames, start, end, strand, sep = ":")) %>%
    dplyr::filter(width >=40, width <=2000)

  write.table(intergenic.5,
              file = paste(main_out_path, "/Intergenic/5_prime/", tissue, ".txt", sep=""),
              row.names = FALSE, quote = FALSE, sep="\t")


  intergenic.3 = rbind(inter.sr.3, inter.nosr.3) %>%
    mutate(ER_id = paste(tissue, seqnames, start, end, strand, sep = ":")) %>%
    dplyr::filter(width >=40, width <=2000)

  write.table(intergenic.3,
              file = paste(main_out_path, "/Intergenic/3_prime/", tissue, ".txt", sep=""),
              row.names = FALSE, quote = FALSE, sep="\t")


}
