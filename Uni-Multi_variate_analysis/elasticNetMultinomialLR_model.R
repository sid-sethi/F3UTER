args<-commandArgs(TRUE)


library(glmnet)
library(caret)
library(e1071)
library(doMC)
registerDoMC(cores = 8)
library(tidyverse)


# output path
path = args[1]
# output names
out = args[2]


data = readRDS(file = "/ML_tables/ml_table.rds")

data$class = recode(data$class, UTR = "three_prime_UTR", "5-UTR" = "five_prime_UTR")

preProcValues <- preProcess(data, method = c("center", "scale"))
data <- predict(preProcValues, data)

# for test data
overall = data.frame()
fn = data.frame()
fp = data.frame()
utr = data.frame()
ice = data.frame()
five = data.frame()
lncrna = data.frame()
ncrna = data.frame()
pseudo = data.frame()


cbind.fill <- function(...){
  nm <- list(...)
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow))
  do.call(cbind, lapply(nm, function (x)
    rbind(x, matrix(, n-nrow(x), ncol(x)))))
}


datax = data[,1:41]
datax$polyA_signal = as.numeric(datax$polyA_signal)
datay = data[,42]

flds <- createFolds(datay, k = 5, list = TRUE, returnTrain = FALSE)

counter = 0
for(j in 1:5){
  counter = counter + 1

  testx = datax[unlist(flds[j]),]
  testy = datay[unlist(flds[j])]

  trainx = datax[-unlist(flds[j]),]
  trainy = datay[-unlist(flds[j])]


  track = paste("Fold number: ",j,sep="")
  cat("\n*************************************\n")
  print(track)
  cat("*************************************\n")



  model <- glmnet(as.matrix(trainx), trainy,
                  family = "multinomial",
                  alpha = 0.5,
                  nlambda = 25, #We don´t need more
                  maxit = 10000)



  #print(model)
  #coef(model)

  lambda.min <- model$lambda %>% min()

  predicted.response <- predict(model, as.matrix(testx), s= lambda.min, type = "class")
  confuseMat.test = confusionMatrix(data=as.factor(predicted.response),reference=testy)
  #confuseMat.test

  res.overall <- confuseMat.test$overall %>% as.data.frame()

  results = t(confuseMat.test$byClass) %>%as.data.frame()
  res.utr <- results %>% dplyr::select(matches("three_prime_UTR"))
  res.five <- results %>% dplyr::select(matches("five_prime_UTR"))
  res.lncrna <- results %>% dplyr::select(matches("lncRNA"))
  res.ncrna <- results %>% dplyr::select(matches(" ncRNA"))
  res.ice <- results %>% dplyr::select(matches("ICE"))
  res.pseudo <- results %>% dplyr::select(matches("pseudoGene"))

  # calculating mis-classification
  table.fn <- confuseMat.test$table %>%
    as.data.frame() %>%
    filter(Reference %in% "three_prime_UTR")

  table.fp <- confuseMat.test$table %>%
    as.data.frame() %>%
    filter(Prediction %in% "three_prime_UTR")

  total.utr.ref <- table.fn$Freq %>% sum()
  total.utr.pred <- table.fp$Freq %>% sum()

  fn.table <- table.fn %>% mutate(fn.rate = (Freq/total.utr.ref)*100) %>%
    column_to_rownames("Prediction") %>%
    dplyr::select(fn.rate)
  fp.table <- table.fp %>% mutate(fp.rate = (Freq/total.utr.pred)*100) %>%
    column_to_rownames("Reference") %>%
    dplyr::select(fp.rate)


  overall <- cbind.fill(overall, res.overall)
  utr <- cbind.fill(utr, res.utr)
  ice <- cbind.fill(ice, res.ice)
  lncrna <- cbind.fill(lncrna, res.lncrna)
  ncrna <- cbind.fill(ncrna, res.ncrna)
  five <- cbind.fill(five, res.five)
  pseudo <- cbind.fill(pseudo, res.pseudo)
  fn <- cbind.fill(fn, fn.table)
  fp <- cbind.fill(fp, fp.table)
}


# writing the results to files
# stats from confusion matrix
overall.df <- data.frame(unlist(overall),stringsAsFactors=FALSE)
file1 = paste(path, out,".overall.txt",sep="")
write.table(overall.df, file = file1, sep="\t",quote=F,row.names=T, col.names=F)

utr.df <- data.frame(unlist(utr),stringsAsFactors=FALSE)
file2 = paste(path, out,".utr.txt",sep="")
write.table(utr.df, file = file2, sep="\t",quote=F,row.names=T, col.names=F)

ice.df <- data.frame(unlist(ice),stringsAsFactors=FALSE)
file3 = paste(path, out,".ice.txt",sep="")
write.table(ice.df, file = file3, sep="\t",quote=F,row.names=T, col.names=F)

lncrna.df <- data.frame(unlist(lncrna),stringsAsFactors=FALSE)
file4 = paste(path, out,".lncrna.txt",sep="")
write.table(lncrna.df, file = file4, sep="\t",quote=F,row.names=T, col.names=F)

ncrna.df <- data.frame(unlist(ncrna),stringsAsFactors=FALSE)
file5 = paste(path, out,".ncrna.txt",sep="")
write.table(ncrna.df, file = file5, sep="\t",quote=F,row.names=T, col.names=F)

five.df <- data.frame(unlist(five),stringsAsFactors=FALSE)
file6 = paste(path, out,".five.txt",sep="")
write.table(five.df, file = file6, sep="\t",quote=F,row.names=T, col.names=F)

pseudo.df <- data.frame(unlist(pseudo),stringsAsFactors=FALSE)
file7 = paste(path, out,".pseudo.txt",sep="")
write.table(pseudo.df, file = file7, sep="\t",quote=F,row.names=T, col.names=F)

fn.df <- data.frame(unlist(fn),stringsAsFactors=FALSE)
file8 = paste(path, out,".fn.txt",sep="")
write.table(fn.df, file = file8, sep="\t",quote=F,row.names=T, col.names=F)

fp.df <- data.frame(unlist(fp),stringsAsFactors=FALSE)
file9 = paste(path, out,".fp.txt",sep="")
write.table(fp.df, file = file9, sep="\t",quote=F,row.names=T, col.names=F)
