args<-commandArgs(TRUE)

suppressPackageStartupMessages({
  library(randomForest)
  library(ROCR)
  library(caret)
  library(tidyverse)
})


feature_table = args[1]
outpath = args[2]
tissue = args[3]
prediction_model = args[4] # F3UTER saved model rda file

tissues = c(tissue)


############## MODEL ######################
# loading MODEL
load(prediction_model)
###########################################


for(tissue in tissues){

  print(tissue)

  ####################### RUNNING ON ER DATA #########################

  ml_table2 = readRDS(file = feature_table)

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


  file1 = str_c(outpath, "/", tissue, "_predictions.txt")
  write.table(res, file = file1, sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE)

  rm(test, res, file1)

}
