args<-commandArgs(TRUE) # 1st  = ML directory , 2nd = F3UTER trained model rda file

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
  "GCB1", "GCB2"
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


  ################ predicting test data #######################
  test$predicted.response <- predict(rf_model, test)
  test$predicted.prob = predict(rf_model, test, type = "prob")[,2]

  res = test %>%
    rownames_to_column("id") %>%
    dplyr::select(c(id, predicted.response, predicted.prob))


  file1 = paste(out.path, "/", tissue,".pred.txt",sep="")
  write.table(res, file = file1, sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE)

  rm(test, res, file1)

}
