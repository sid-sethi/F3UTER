args<-commandArgs(TRUE) # 1 = output path ; 2 = ML table

library(randomForest)
library(caret)
library(doMC)
registerDoMC(cores = 10)
library(tidyverse)

# output path
out = args[1]

if(!dir.exists(paste(out, "/Trained_model", sep=""))){
  system(paste("mkdir -m a=rwx ",out, "/Trained_model", sep=""))
}

# output path
out.path = paste(out, "/Trained_model", sep="")


####################### TRAINING DATA = ALL DATA #########################

rds = args[2]

ml_table = readRDS(file = args[2])

ml_table$class = ifelse(ml_table$class %in% "UTR", "UTR", "Non_3_UTR") %>% as.factor()

data = ml_table

preProcValues <- preProcess(data, method = c("center", "scale"))
train <- predict(preProcValues, data)

############## MODEL ######################

# main command
cntrl = trainControl(method="repeatedcv", number=5, repeats=5,summaryFunction = twoClassSummary,sampling = "down", allowParallel = TRUE, classProbs = TRUE)
rf_model<-train(class~.,data=train,method="rf",trControl=cntrl, metric = "ROC",importance = TRUE)

save(rf_model, file = paste(out.path, "/f3uter_trained_model.rda", sep=""))
