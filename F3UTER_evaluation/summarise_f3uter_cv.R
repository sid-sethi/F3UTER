args<-commandArgs(TRUE) # 1st argument: ML output directory

library(rlist)
library(matrixStats)
library(ROCR)
library(e1071)
library(PRROC)
library(tidyverse)

outPath = args[1]

system(paste("mkdir -m a=rwx ",outPath, "/Summary", sep=""))


auc <- list.files(path = outPath, pattern = "\\.auc.txt", recursive = TRUE, full.names = TRUE) %>%
  map(read.table, header=F) %>%
  bind_rows()
file1 = paste(outPath, "/Summary/auc.txt",sep="")
write.table(auc, file=file1, sep="\t",quote=F,row.names=F, col.names=F)

prauc <- list.files(path = outPath, pattern = "\\.prauc.txt", recursive = TRUE, full.names = TRUE) %>%
  map(read.table, header=F) %>%
  bind_rows()
file2 = paste(outPath, "/Summary/prauc.txt",sep="")
write.table(prauc, file=file2, sep="\t",quote=F,row.names=F, col.names=F)

decAcu <- list.files(path = outPath, pattern = "\\.decAccuracy.txt", recursive = TRUE, full.names = TRUE) %>%
  map(read.table, header=F, row.names = 1) %>%
  list.cbind() %>%
  rownames_to_column("features") %>%
  mutate(mean = rowMeans(.[,-1]), sd = rowSds(as.matrix(.[,-1])))
file3 = paste(outPath, "/Summary/decAccuracy.txt",sep="")
write.table(decAcu, file=file3, sep="\t",quote=F,row.names=F, col.names=T)

decGini <- list.files(path = outPath, pattern = "\\.decGini.txt", recursive = TRUE, full.names = TRUE) %>%
  map(read.table, header=F, row.names = 1) %>%
  list.cbind() %>%
  rownames_to_column("features") %>%
  mutate(mean = rowMeans(.[,-1]), sd = rowSds(as.matrix(.[,-1])))
file4 = paste(outPath, "/Summary/decGini.txt",sep="")
write.table(decGini, file=file4, sep="\t",quote=F,row.names=F, col.names=T)

stats <- list.files(path = outPath, pattern = "\\.results.txt", recursive = TRUE, full.names = TRUE) %>%
  map(read.table, header=F, row.names = 1) %>%
  list.cbind() %>%
  rownames_to_column("stats") %>%
  mutate(mean = rowMeans(.[,-1]), sd = rowSds(as.matrix(.[,-1])))
file5 = paste(outPath, "/Summary/stats.txt",sep="")
write.table(stats, file=file5, sep="\t",quote=F,row.names=F, col.names=T)


# combine the prediction lists RDS
res.pred <- list.files(path = outPath, pattern = "\\.resPred.rds", recursive = TRUE, full.names = TRUE) %>%
  sapply(.,readRDS, USE.NAMES=FALSE)
attributes(res.pred) <- NULL

res.labels <- list.files(path = outPath, pattern = "\\.resLabels.rds", recursive = TRUE, full.names = TRUE) %>%
  sapply(.,readRDS, USE.NAMES=FALSE)
attributes(res.labels) <- NULL

########################################################################
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

rocs.x = data.frame()
for(q in 1:length(perf.avg.data@x.values)){
  values = unlist(perf.avg.data@x.values[q]) %>% as.data.frame() %>%
    `colnames<-` (c("x"))
  rocs.x = rbind(rocs.x,values)
}

rocs.y = data.frame()
for(q in 1:length(perf.avg.data@y.values)){
  values = unlist(perf.avg.data@y.values[q]) %>%
    as.data.frame() %>%
    mutate(type1 = paste("fold",q, sep=""), type2 = "cv") %>%
    plyr::rename(c("." = "y"))
  rocs.y = rbind(rocs.y,values)
}


roc.all = cbind(rocs.x,rocs.y)
file6 = paste(outPath, "/Summary/rocValues.cv.txt",sep="")
write.table(roc.all, file=file6, sep="\t",quote=F,row.names=F)

rocs.y.wide = data.frame(perf.avg.data@y.values)
avg.y = rowMeans(rocs.y.wide)
rocs.y.wide = transform(rocs.y.wide, sd=apply(rocs.y.wide,1, sd))
rocs.y.wide$mean = avg.y
file7 = paste(outPath, "/Summary/rocMeanSd.txt",sep="")
write.table(rocs.y.wide, file=file7, sep="\t",quote=F,row.names=F)

mean.x = unlist(perf.rocr.avg@x.values) %>% as.data.frame() %>%  `colnames<-` (c("x"))
mean.y = unlist(perf.rocr.avg@y.values) %>% as.data.frame() %>%
  mutate(type1 = "mean", type2 = "mean") %>%
  plyr::rename(c("." = "y"))
mean.all = cbind(mean.x,mean.y)
file8 = paste(outPath, "/Summary/rocValues.mean.txt", sep="")
write.table(mean.all, file=file8, sep="\t",quote=F,row.names=F)


######################################################################
######################## PR curve ####################################
perf.pr.all = performance(pred.all,"prec","rec")
avg.results(perf.pr.all)

pr.x = data.frame()
for(q in 1:length(perf.avg.data@x.values)){
  values = unlist(perf.avg.data@x.values[q]) %>% as.data.frame() %>%
    `colnames<-` (c("x"))
	pr.x = rbind(pr.x,values)
}

pr.y = data.frame()
for(q in 1:length(perf.avg.data@y.values)){
  values = unlist(perf.avg.data@y.values[q]) %>%
    as.data.frame() %>%
    mutate(type1 = paste("fold",q, sep=""), type2 = "cv") %>%
    plyr::rename(c("." = "y"))
	pr.y = rbind(pr.y,values)
}

pr.all = cbind(pr.x,pr.y)
file9 = paste(outPath, "/Summary/prValues.cv.txt",sep="")
write.table(pr.all, file=file9, sep="\t",quote=F,row.names=F)

pr.y.wide = data.frame(perf.avg.data@y.values)
avg.pr.y = rowMeans(pr.y.wide)
pr.y.wide = transform(pr.y.wide, sd=apply(pr.y.wide,1, sd))
pr.y.wide$mean = avg.pr.y
file10 = paste(outPath,"/Summary/prMeanSd.txt",sep="")
write.table(pr.y.wide,file=file10,sep="\t",quote=F,row.names=F)

mean.pr.x = unlist(perf.rocr.avg@x.values) %>% as.data.frame() %>%  `colnames<-` (c("x"))
mean.pr.y = unlist(perf.rocr.avg@y.values) %>% as.data.frame() %>%
  mutate(type1 = "mean", type2 = "mean") %>%
  plyr::rename(c("." = "y"))
mean.pr.all = cbind(mean.pr.x,mean.pr.y)
file11 = paste(outPath,"/Summary/prValues.mean.txt",sep="")
write.table(mean.pr.all, file=file11, sep="\t",quote=F,row.names=F)


###########################################
######## Other plots ######################
###########################################

format_data_from_rocr_object <- function(perf.obj){

  x = data.frame()
  for(q in 1:length(perf.obj@x.values)){
    values = unlist(perf.obj@x.values[q]) %>% as.data.frame() %>%
      `colnames<-` (c("x"))
    x = rbind(x,values)
  }

  y = data.frame()
  for(q in 1:length(perf.obj@y.values)){
    values = unlist(perf.obj@y.values[q]) %>%
      as.data.frame() %>%
      mutate(type1 = paste("fold",q, sep=""), type2 = "cv") %>%
      plyr::rename(c("." = "y"))
    y = rbind(y,values)
  }
  all = cbind(x,y)
  return(all)
}


perf.mat = performance(pred.all,"mat")
mat.all = format_data_from_rocr_object(perf.mat)
mat.all.subset <- mat.all %>% filter(x %in% c(0.02, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
mat.mean <- mat.all %>% group_by(x) %>% summarise(mean = mean(y))
file12 = paste(outPath,"/Summary/mcc_vs_cutoff.png",sep="")
png(file12,bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
#plot(perf.mat, avg="vertical", spread.estimate="boxplot", col="blue")
ggplot(data=mat.all.subset, aes(x=factor(x), y=y)) +
  geom_line(data = mat.mean, aes(x=as.factor(x), y=mean, group = 1), color = "red") +
  geom_boxplot(lwd = 0.55, alpha=0.8, fill="red", width = 20, outlier.shape=21, outlier.size=0.8, outlier.fill = "grey") +
  theme_bw() +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.40,"cm"),
        legend.title = element_text(size=11),
        legend.text = element_text(size = 11),
        legend.position= c(0.85,0.20),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        panel.border = element_blank(),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_x_discrete(name = "Prediction cutoff", breaks = seq(0,1, by=0.1)) +
  scale_y_continuous(name = "Matthews correlation coefficient", limits = c(0,1))
dev.off()



perf.err = performance(pred.all,"err")
err.all = format_data_from_rocr_object(perf.err)
err.all.subset <- err.all %>% filter(x %in% c(0.02, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
err.mean <- err.all %>% group_by(x) %>% summarise(mean = mean(y))
file13 = paste(outPath,"/Summary/err_vs_cutoff.png",sep="")
png(file13,bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
#plot(perf.err, avg="vertical", spread.estimate="boxplot", col="blue")
ggplot(data=err.all.subset, aes(x=factor(x), y=y)) +
  geom_line(data = err.mean, aes(x=as.factor(x), y=mean, group = 1), color = "red") +
  geom_boxplot(lwd = 0.55, alpha=0.8, fill="red", width = 20, outlier.shape=21, outlier.size=0.8, outlier.fill = "grey") +
  theme_bw() +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.40,"cm"),
        legend.title = element_text(size=11),
        legend.text = element_text(size = 11),
        legend.position= c(0.85,0.20),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        panel.border = element_blank(),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_x_discrete(name = "Prediction cutoff", breaks = seq(0,1, by=0.1)) +
  scale_y_continuous(name = "Error rate", limits = c(0,1))
dev.off()



perf.f = performance(pred.all,"f")
f.all = format_data_from_rocr_object(perf.f)
f.all.subset <- f.all %>% filter(x %in% c(0.02, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
f.mean <- f.all %>% group_by(x) %>% summarise(mean = mean(y))
file14 = paste(outPath,"/Summary/f_vs_cutoff.png",sep="")
png(file14,bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
#plot(perf.f, avg="vertical", spread.estimate="boxplot", col="blue")
ggplot(data=f.all.subset, aes(x=factor(x), y=y)) +
  geom_line(data = f.mean, aes(x=as.factor(x), y=mean, group = 1), color = "red") +
  geom_boxplot(lwd = 0.55, alpha=0.8, fill="red", width = 20, outlier.shape=21, outlier.size=0.8, outlier.fill = "grey") +
  theme_bw() +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.40,"cm"),
        legend.title = element_text(size=11),
        legend.text = element_text(size = 11),
        legend.position= c(0.85,0.20),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        panel.border = element_blank(),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_x_discrete(name = "Prediction cutoff", breaks = seq(0,1, by=0.1)) +
  scale_y_continuous(name = "Precision-recall f-measure", limits = c(0,1))
dev.off()



perf.acc = performance(pred.all,"acc")
acc.all = format_data_from_rocr_object(perf.acc)
acc.all.subset <- acc.all %>% filter(x %in% c(0.02, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
acc.mean <- acc.all %>% group_by(x) %>% summarise(mean = mean(y))
file15 = paste(outPath,"/Summary/acc_vs_cutoff.png",sep="")
png(file15,bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
#plot(perf.acc, avg="vertical", spread.estimate="boxplot", col="blue")
ggplot(data=acc.all.subset, aes(x=factor(x), y=y)) +
  geom_line(data = acc.mean, aes(x=as.factor(x), y=mean, group = 1), color = "red") +
  geom_boxplot(lwd = 0.55, alpha=0.8, fill="red", width = 20, outlier.shape=21, outlier.size=0.8, outlier.fill = "grey") +
  theme_bw() +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.40,"cm"),
        legend.title = element_text(size=11),
        legend.text = element_text(size = 11),
        legend.position= c(0.85,0.20),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        panel.border = element_blank(),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_x_discrete(name = "Prediction cutoff", breaks = seq(0,1, by=0.1)) +
  scale_y_continuous(name = "Accuracy", limits = c(0,1))
dev.off()


perf.sens = performance(pred.all,"sens")
sens.all = format_data_from_rocr_object(perf.sens)
sens.all.subset <- sens.all %>% filter(x %in% c(0.02, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
sens.mean <- sens.all %>% group_by(x) %>% summarise(mean = mean(y))
file16 = paste(outPath,"/Summary/sens_vs_cutoff.png",sep="")
png(file16,bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
#plot(perf.sens, avg="vertical", spread.estimate="boxplot", col="blue")
ggplot(data=sens.all.subset, aes(x=factor(x), y=y)) +
  geom_line(data = sens.mean, aes(x=as.factor(x), y=mean, group = 1), color = "red") +
  geom_boxplot(lwd = 0.55, alpha=0.8, fill="red", width = 20, outlier.shape=21, outlier.size=0.8, outlier.fill = "grey") +
  theme_bw() +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.40,"cm"),
        legend.title = element_text(size=11),
        legend.text = element_text(size = 11),
        legend.position= c(0.85,0.20),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        panel.border = element_blank(),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_x_discrete(name = "Prediction cutoff", breaks = seq(0,1, by=0.1)) +
  scale_y_continuous(name = "Sensitivity", limits = c(0,1))
dev.off()


perf.spec = performance(pred.all,"spec")
spec.all = format_data_from_rocr_object(perf.spec)
spec.all.subset <- spec.all %>% filter(x %in% c(0.02, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
spec.mean <- spec.all %>% group_by(x) %>% summarise(mean = mean(y))
file17 = paste(outPath,"/Summary/spec_vs_cutoff.png",sep="")
png(file17,bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
#plot(perf.spec, avg="vertical", spread.estimate="boxplot", col="blue")
ggplot(data=spec.all.subset, aes(x=factor(x), y=y)) +
  geom_line(data = spec.mean, aes(x=as.factor(x), y=mean, group = 1), color = "red") +
  geom_boxplot(lwd = 0.55, alpha=0.8, fill="red", width = 20, outlier.shape=21, outlier.size=0.8, outlier.fill = "grey") +
  theme_bw() +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.40,"cm"),
        legend.title = element_text(size=11),
        legend.text = element_text(size = 11),
        legend.position= c(0.85,0.20),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        panel.border = element_blank(),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_x_discrete(name = "Prediction cutoff", breaks = seq(0,1, by=0.1)) +
  scale_y_continuous(name = "Specificity", limits = c(0,1))
dev.off()
