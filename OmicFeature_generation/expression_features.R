library(tidyverse)
library(derfinder)

options(stringsAsFactors=F)


####################################################################################################################
###################################################      Functions      ############################################
####################################################################################################################

# mean coverage
meanCov <- function(start, end){
  mean(cov[start:end])
}

# reads out / reads in (ratio)
RoRi <- function(start, end, strand){
  inn = ifelse(strand %in% "+", start, end)
  out = ifelse(strand %in% "+", end, start)
  RoRi = cov[out]/cov[inn]
  return(RoRi)
}

# entropy efficiency (EE)
EE <- function(start, end){

  cov.vector = cov[start:end]

  # only consider bases which have coverage > 0
  cov.vector = cov.vector[cov.vector>0]

  length = length(cov.vector)
  total_coverage = sum(cov.vector)
  p.xi = cov.vector/total_coverage
  log.p.xi = log(p.xi)

  # final equation
  ee = -(sum(p.xi*log.p.xi))/log(length)
  return(ee)
}


# percentage difference between reads at boundaries
pd <- function(start, end){
  pd = (abs(cov[start] - cov[end])/mean(c(cov[start], cov[end])))*100
  return(round(pd,3))
}

####################################################################################################################


# load regions for which the feature need to be calculated
# loading the utr training regions
load("/GS.RData")
# this code can be modified to be used on other training regions and on ERs


tissues = c(
  "adipose_subcutaneous",
  "adrenal_gland",
  "artery_aorta",
  "artery_coronary",
  "artery_tibial",
  "adipose_visceral_omentum",
  ####### "bladder",
  "brain_amygdala",
  "brain_anterior_cingulate_cortex_ba24",
  "brain_caudate_basal_ganglia",
  "brain_cerebellar_hemisphere",
  # "brain_cerebellum",
  "brain_cortex",
  "brain_frontal_cortex_ba9",
  "brain_hippocampus",
  "brain_hypothalamus",
  "brain_nucleus_accumbens_basal_ganglia",
  "brain_putamen_basal_ganglia",
  "brain_spinal_cord_cervical_c_1",
  "brain_substantia_nigra",
  ####### "cells_ebv_transformed_lymphocytes",
  ####### "cells_transformed_fibroblasts",
  ####### "cervix_ectocervix",
  "colon_sigmoid",
  "colon_transverse",
  "esophagus_gastroesophageal_junction",
  "esophagus_mucosa",
  "esophagus_muscularis",
  ##### "fallopian_tube",
  "heart_atrial_appendage",
  "heart_left_ventricle",
  "kidney_cortex",
  "liver",
  "lung",
  "minor_salivary_gland",
  "muscle_skeletal",
  "nerve_tibial",
  ###### "ovary",
  "pancreas",
  "pituitary",
  ##### "prostate",
  "skin_not_sun_exposed_suprapubic",
  "skin_sun_exposed_lower_leg",
  "small_intestine_terminal_ileum",
  "spleen",
  "stomach",
  ###### "testis",
  "thyroid",
  ###### "uterus",
  ###### "vagina",
  "whole_blood"
)


for(tissue in tissues){

  print(tissue)

  utr.exp <- data.frame()

  for(i in c(1:22,"X","Y")){

    # per chr
    chr.utr = utr %>% dplyr::filter(seqnames %in% i)

    # load exp data from GTEx
    file = paste("/GTEx_data/Mean_coverage/by_tissue_smfrze_use_me/",
                tissue, "/gtex_", tissue, "_chr", i, "_mean_cov.rda", sep="")

    if(i != 1){ rm(tissue_coverage_w_mean_normalised) }

    load(file)
    cov = tissue_coverage_w_mean_normalised$meanCoverage %>% as.numeric()

    # analyse 3'UTR
    result.utr = chr.utr %>%
      rowwise() %>%
      mutate(meanCov = meanCov(start, end),
           #RoRi = RoRi(start, end, strand),
           pd = pd(start, end),
           EE = EE(start, end)) %>%
      as.data.frame()

    utr.exp <- rbind.data.frame(utr.exp, result.utr)

  }

  utr.export = left_join(utr, utr.exp, by = "three_prime_utr_id") %>%
    mutate(class = "3-UTR") %>%
    select(c(three_prime_utr_id, meanCov, RoRi, pd, EE, class)) %>%
    plyr::rename(c("three_prime_utr_id" = "id"))

  file_name = paste("/Features/exp_feat_", tissue, ".txt", sep="")
  write.table(utr.export, file = file_name, row.names = FALSE, quote = FALSE, sep="\t")

}
