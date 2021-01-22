args<-commandArgs(TRUE) # 1st  = ML directory , 2nd = Prediction trained model rda file

library(randomForest)
library(ROCR)
library(caret)
library(e1071)
library(PRROC)
library(doMC)
registerDoMC(cores = 10)
library(tidyverse)

path = args[1]
erFeat = args[3]

if(!dir.exists(paste(path, "/Prediction", sep=""))){
  system(paste("mkdir -m a=rwx ",path, "/Prediction", sep=""))
}

# new output path
out.path = paste(path, "/Prediction", sep="")

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


############## MODEL ######################
# loading MODEL
model = args[2]
load(model)
###########################################


for(tissue in tissues){

  print(tissue)

  ####################### RUNNING ON ER DATA #########################

  ml_table2 = readRDS(file = paste(erFeat, "/ML_tables/table_", tissue, ".rds", sep=""))

  data2 = ml_table2 %>%
    plyr::rename(c("pd" = "meanPd", "EE" = "meanEE"))

  preProcValues <- preProcess(data2, method = c("center", "scale"))
  test <- predict(preProcValues, data2)


  ################ predicting #######################
  test$predicted.response <- predict(rf_model, test)
  test$predicted.prob = predict(rf_model, test, type = "prob")[,2]

  res = test %>%
    rownames_to_column("id") %>%
    dplyr::select(c(id, predicted.response, predicted.prob))


  file1 = paste(out.path, "/", tissue,".pred.txt",sep="")
  write.table(res, file = file1, sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE)

  rm(test, res, file1)

}
