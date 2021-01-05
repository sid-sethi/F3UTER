args<-commandArgs(TRUE)

library(randomForest)
library(ROCR)
library(caret)
library(e1071)
library(PRROC)
library(doMC)
registerDoMC(cores = 8)
library(tidyverse)


# output path
path = args[1]
# output names
out = args[2]


ml_table = readRDS(file = "/ML_tables/ml_table.rds")
ml_table$class = ifelse(ml_table$class %in% "UTR", "UTR", "Non_3_UTR") %>% as.factor()

data = ml_table

#####################################################
######## for bias-variance analysis ##################
#fraction = 0.8
#data = ml_table %>% dplyr::sample_frac(fraction, replace=FALSE)
########################################################


preProcValues <- preProcess(data, method = c("center", "scale"))
data <- predict(preProcValues, data)

# for test data
cm = data.frame()
auc.roc = data.frame()
auc.pr = data.frame()
dec.accu = data.frame()
dec.gini = data.frame()
res.pred = list()
res.labels = list()

# for training data
cm.train = data.frame()
auc.roc.train = data.frame()
auc.pr.train = data.frame()
dec.accu.train = data.frame()
dec.gini.train = data.frame()
res.pred.train = list()
res.labels.train = list()

cbind.fill <- function(...){
  nm <- list(...)
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow))
  do.call(cbind, lapply(nm, function (x)
    rbind(x, matrix(, n-nrow(x), ncol(x)))))
}


flds <- createFolds(data$class, k = 5, list = TRUE, returnTrain = FALSE)

counter = 0
for(j in 1:5){

  counter = counter + 1
  test = data[unlist(flds[j]),]
  train = data[-unlist(flds[j]),]


  track = paste("Fold number: ",j,sep="")
  cat("\n*************************************\n")
  print(track)
  cat("*************************************\n")


  cntrl = trainControl(method="none",sampling = "down", allowParallel = TRUE, classProbs = TRUE)
  rf_model<-train(class~.,data=train,method="rf",metric = "Kappa", trControl=cntrl,importance = TRUE)

  ################ importance of variables #############
  var.imp = importance(rf_model$finalModel)
  dec.accu = cbind.fill(dec.accu,var.imp[,3])
  dec.gini = cbind.fill(dec.gini,var.imp[,4])

  ################ predicting test data and producing confusion Matrix #######################
  test$predicted.response <- predict(rf_model, test)
  confuseMat.test = confusionMatrix(data=test$predicted.response,reference=test$class, positive = "UTR")
  #confuseMat.test
  results = t(data.frame(cbind(t(confuseMat.test$byClass),t(confuseMat.test$overall))))
  accuracy = round(results[12],2)
  sens = round(results[1],2)
  spec = round(results[2],2)
  table = confuseMat.test$table

  cm = cbind.fill(cm,results[,1])

  ###################### validating test data ####################

  ## For producing probabilities ####
  test.pr = predict(rf_model, test, type = "prob")[,2]
  test.pred = prediction(test.pr,test$class)

  test.AUC=performance(test.pred,"auc") #Calculate the AUC value
  AUC=test.AUC@y.values[[1]]

  # storing ROCR curve values
  auc.roc[counter,"ROCAUC"] = AUC
  res.pred[[counter]] <- test.pr
  res.labels[[counter]] <- test$class

  #PR curve
  fg = test.pr[test$class == "UTR"]
  bg = test.pr[test$class == "Non_3_UTR"]
  pr = pr.curve(scores.class0 = fg, scores.class1 = bg, curve=T)
  prauc = pr$auc.integral

  # storing PR curve values
  auc.pr[counter,"PRAUC"] = prauc


  #######################################################################
  ############## predicting on the TRAINING data itself #################

  ################ predicting test data and producing confusion Matrix #######################
  train$predicted.response <- predict(rf_model, train)
  confuseMat.train = confusionMatrix(data=train$predicted.response,reference=train$class, positive = "UTR")
  results.train = t(data.frame(cbind(t(confuseMat.train$byClass),t(confuseMat.train$overall))))
  accuracy.train = round(results.train[12],2)
  sens.train = round(results.train[1],2)
  spec.train = round(results.train[2],2)
  table.train = confuseMat.train$table

  cm.train = cbind.fill(cm.train,results.train[,1])

  ## For producing probabilities ####
  train.pr = predict(rf_model, train, type = "prob")[,2]
  train.pred = prediction(train.pr,train$class)

  train.AUC=performance(train.pred,"auc") #Calculate the AUC value
  AUC.train=train.AUC@y.values[[1]]

  # storing ROCR curve values
  auc.roc.train[counter,"ROCAUC"] = AUC.train
  res.pred.train[[counter]] <- train.pr
  res.labels.train[[counter]] <- train$class

  #PR curve
  fg.train = train.pr[train$class == "UTR"]
  bg.train = train.pr[train$class == "Non_3_UTR"]
  pr.train = pr.curve(scores.class0 = fg.train, scores.class1 = bg.train, curve=T)
  prauc.train = pr.train$auc.integral

  # storing PR curve values
  auc.pr.train[counter,"PRAUC"] = prauc.train

}


# writing the results to files
# Mean decrease accuracy
accu.df <- data.frame(unlist(dec.accu),stringsAsFactors=FALSE)
file1 = paste(path, out,".decAccuracy.txt",sep="")
write.table(accu.df,file = file1, sep="\t",quote=F, col.names=F)

# Mean decrease Gini
gini.df <- data.frame(unlist(dec.gini),stringsAsFactors=FALSE)
file2 = paste(path, out,".decGini.txt",sep="")
write.table(gini.df,file = file2, sep="\t",quote=F, col.names=F)

# stats from confusion matrix
cm.df <- data.frame(unlist(cm),stringsAsFactors=FALSE)
file3 = paste(path, out,".results.txt",sep="")
write.table(cm.df, file = file3, sep="\t",quote=F,row.names=T, col.names=F)

cm.df.train <- data.frame(unlist(cm.train),stringsAsFactors=FALSE)
file8 = paste(path, out,".results.train.txt",sep="")
write.table(cm.df.train, file = file8, sep="\t",quote=F,row.names=T, col.names=F)


# auc for ROC and PR
file4 = paste(path, out,".auc.txt",sep="")
write.table(auc.roc, file=file4, sep="\t",quote=F,row.names=F, col.names=F)
file5 = paste(path, out,".prauc.txt",sep="")
write.table(auc.pr, file=file5, sep="\t",quote=F,row.names=F, col.names=F)

file9 = paste(path, out,".auc.train.txt",sep="")
write.table(auc.roc.train, file=file9, sep="\t",quote=F,row.names=F, col.names=F)
file10 = paste(path, out,".prauc.train.txt",sep="")
write.table(auc.pr.train, file=file10, sep="\t",quote=F,row.names=F, col.names=F)


# saving results lists
file6 = paste(path, out,".resPred.rds",sep="")
saveRDS(res.pred, file= file6)

file7 = paste(path, out,".resLabels.rds",sep="")
saveRDS(res.labels, file= file7)

file11 = paste(path, out,".resPred.train.rds",sep="")
saveRDS(res.pred.train, file= file11)

file12 = paste(path, out,".resLabels.train.rds",sep="")
saveRDS(res.labels.train, file= file12)
