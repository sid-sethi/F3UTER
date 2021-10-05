args<-commandArgs(TRUE)

suppressPackageStartupMessages({
  library(randomForest)
  library(ROCR)
  library(caret)
  library(e1071)
  library(PRROC)
  library(tidyverse)
})

table = args[1]
model = args[2]
outpath = args[3]
tissue = args[4]

############## MODEL ######################
# loading MODEL
load(model)
###########################################


  print(str_c(Sys.time(), "-", tissue))

  ####################### RUNNING ON ER DATA #########################

  ml_table2 = readRDS(file = table)

  data2 = ml_table2 %>%
    plyr::rename(c("pd" = "meanPd", "EE" = "meanEE"))

  #data2 = ml_table2 %>%
  #  dplyr::select(-c(pd, EE))

  preProcValues <- preProcess(data2, method = c("center", "scale"))
  test <- predict(preProcValues, data2)


  ################ predicting test data #######################
  test$predicted.response <- predict(rf_model, test)
  test$predicted.prob = predict(rf_model, test, type = "prob")[,2]

  res = test %>%
    rownames_to_column("id") %>%
    dplyr::select(c(id, predicted.response, predicted.prob))


  file1 = str_c(outpath, "/", tissue, ".pred.txt")
  write.table(res, file = file1, sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE)
