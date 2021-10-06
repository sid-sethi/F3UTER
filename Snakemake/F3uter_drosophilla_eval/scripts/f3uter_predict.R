args<-commandArgs(TRUE)

suppressPackageStartupMessages({
  library(ROCR)
  library(caret)
  library(e1071)
  library(PRROC)
  library(tidyverse)
  library(rlist)
  library(matrixStats)
})


feature_table = args[1]
outpath = args[2]
tissue = args[3]
prediction_model = args[4] # F3UTER saved model rda file


############## MODEL ######################
# loading MODEL
load(prediction_model)
###########################################


ml_table2 = readRDS(file = feature_table)
ml_table2$class = ifelse(ml_table2$class %in% "three_prime", "UTR", "Non_3_UTR") %>% as.factor()

data2 = ml_table2 %>%
plyr::rename(c("pd" = "meanPd", "EE" = "meanEE"))

preProcValues <- preProcess(data2, method = c("center", "scale"))
test <- predict(preProcValues, data2)


################ predicting test data #######################

test$predicted.response <- predict(rf_model, test)
test$predicted.prob = predict(rf_model, test, type = "prob")[,2]

res = test %>%
rownames_to_column("id") %>%
dplyr::select(c(id, class, predicted.response, predicted.prob))

file1 = str_c(outpath, "/", tissue, "_predictions.txt")
write.table(res, file = file1, sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE)


# confusion matrix stats
confusemat = confusionMatrix(data=test$predicted.response, reference=test$class, positive = "UTR")
results = t(data.frame(cbind(t(round(confusemat$byClass,3)),t(round(confusemat$overall,3))))) %>%
as.data.frame() %>%
tibble::rownames_to_column("metric") %>%
plyr::rename(c("V1" = "F3UTER"))

tn <- confusemat$table[1]
fp <- confusemat$table[2]
fn <- confusemat$table[3]
tp <- confusemat$table[4]
confusetable <- data.frame(metric = c("tp", "fn", "fp", "tn"), value = c(tp, fn, fp, tn))

write.table(results, file = str_c(outpath, "/", tissue, "_performance.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
write.table(confusetable, file = str_c(outpath, "/", tissue, "_confusemat.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")


## For producing ROC/PR ####

auc = data.frame()
res.pred = list()
res.labels = list()
counter <- 1

test.pr = predict(rf_model, test, type = "prob")[,2]
test.pred = prediction(test.pr,test$class)

test.AUC=performance(test.pred,"auc") #Calculate the AUC value
AUC=test.AUC@y.values[[1]]

# storing ROCR curve values
auc[counter,"ROCAUC"] = AUC
res.pred[[counter]] <- test.pr
res.labels[[counter]] <- test$class

#PR curve
fg = test.pr[test$class == "UTR"]
bg = test.pr[test$class == "Non_3_UTR"]
pr = pr.curve(scores.class0 = fg, scores.class1 = bg, curve=T)
prauc = pr$auc.integral

# storing PR curve values
auc[counter,"PRAUC"] = prauc

write.table(auc, file = str_c(outpath, "/", tissue, ".auc.txt"), sep="\t",quote=F,row.names=FALSE, col.names=TRUE)


######################## ROCR curve ####################################

pred.all = prediction(res.pred,res.labels)
perf.all = performance(pred.all,"tpr","fpr")

# function to average results - modified function from ROCR package
avg.results <- function (perf)
{
  ## for infinite cutoff, assign maximal finite cutoff + mean difference
  ## between adjacent cutoff pairs
  if (length(perf@alpha.values)!=0) perf@alpha.values <-
      lapply(perf@alpha.values,
             function(x) { isfin <- is.finite(x);
             x[is.infinite(x)] <-
               (max(x[isfin]) +
                  mean(abs(x[isfin][-1] -
                             x[isfin][-length(x[isfin])])));
             x } )
  ## remove samples with x or y not finite
  for (i in 1:length(perf@x.values)) {
    ind.bool <- (is.finite(perf@x.values[[i]]) &
                   is.finite(perf@y.values[[i]]))

    if (length(perf@alpha.values)>0)
      perf@alpha.values[[i]] <- perf@alpha.values[[i]][ind.bool]

    perf@x.values[[i]] <- perf@x.values[[i]][ind.bool]
    perf@y.values[[i]] <- perf@y.values[[i]][ind.bool]
  }

  perf.sampled <- perf
  alpha.values <- rev(seq(min(unlist(perf@alpha.values)), max(unlist(perf@alpha.values)),
                          length = max(sapply(perf@alpha.values, length))))
  for (i in 1:length(perf.sampled@y.values)) {
    perf.sampled@x.values[[i]] <- approxfun(perf@alpha.values[[i]],
                                            perf@x.values[[i]], rule = 2, ties = mean)(alpha.values)
    perf.sampled@y.values[[i]] <- approxfun(perf@alpha.values[[i]],
                                            perf@y.values[[i]], rule = 2, ties = mean)(alpha.values)
  }
  perf.avg <- perf.sampled
  perf.avg.data <<- perf.avg
  perf.avg@x.values <- list(rowMeans(data.frame(perf.avg@x.values)))
  perf.avg@y.values <- list(rowMeans(data.frame(perf.avg@y.values)))
  perf.avg@alpha.values <- list(alpha.values)
  perf.rocr.avg <<- perf.avg
}



avg.results(perf.all)

mean.roc.x = unlist(perf.rocr.avg@x.values) %>% as.data.frame() %>%  `colnames<-` (c("x"))
mean.roc.y = unlist(perf.rocr.avg@y.values) %>% as.data.frame() %>%
  mutate(type1 = tissue, type2 = tissue) %>%
  plyr::rename(c("." = "y"))

roc.all = cbind(mean.roc.x, mean.roc.y)
write.table(roc.all, file = str_c(outpath, "/", tissue, "_rocValues.mean.txt"), sep="\t",quote=FALSE,row.names=FALSE)


######################################################################
######################## PR curve ####################################
perf.pr.all = performance(pred.all,"prec","rec")

avg.results(perf.pr.all)

mean.pr.x = unlist(perf.rocr.avg@x.values) %>% as.data.frame() %>%  `colnames<-` (c("x"))
mean.pr.y = unlist(perf.rocr.avg@y.values) %>% as.data.frame() %>%
  mutate(type1 = tissue, type2 = tissue) %>%
  plyr::rename(c("." = "y"))

mean.pr.all = cbind(mean.pr.x, mean.pr.y)
write.table(mean.pr.all, file = str_c(outpath, "/", tissue, "_prValues.mean.txt"), sep="\t",quote=FALSE,row.names=FALSE)
