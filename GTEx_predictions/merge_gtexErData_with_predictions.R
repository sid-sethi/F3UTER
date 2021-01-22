args<-commandArgs(TRUE) # 1st  = ML output folder, 2nd = GTEx erPath

library(tidyverse)
library(stringr)
library(rtracklayer)

options(stringsAsFactors=F)

path = args[1]
erPath = args[2]

if(!dir.exists(paste(path, "/PredictionMerge", sep=""))){
  system(paste("mkdir -m a=rwx ",path, "/PredictionMerge", sep=""))
}

if(!dir.exists(paste(path, "/PredictionMerge_AllTissues", sep=""))){
  system(paste("mkdir -m a=rwx ",path, "/PredictionMerge_AllTissues", sep=""))
}

# new output path
out.path = paste(path, "/PredictionMerge", sep="")
out.path2 = paste(path, "/PredictionMerge_AllTissues", sep="")

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
  "brain_cerebellum",
  "brain_cortex",
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


res_all_tissues = data.frame()

for(tissue in tissues){

  print(str_c(Sys.time(), " - ", tissue))

  #ER data
  er = read.table(paste(erPath, "/", tissue, ".txt", sep=""), header=T, sep="\t")

  ## prediction data
  pred = read.table(paste(path, "/Prediction/", tissue,".pred.txt", sep=""), header=TRUE)

  data = left_join(pred, er, by = c("id" = "ER_id"))

  res_all_tissues = rbind(res_all_tissues,data)

  saveRDS(data, file = paste(out.path, "/", tissue, "_mergedData.rds", sep=""))

  rm(pred, data)

}

saveRDS(res_all_tissues, file = paste(out.path2, "/all_tissues_mergedData.rds", sep=""))
