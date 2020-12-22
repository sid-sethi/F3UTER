#!/bin/Rscript
args<-commandArgs(TRUE)   #1st argument =  output directory path
setwd(args[1])
library(ggplot2)

##### ROCR ##############
file1 = "rocValues.cv.txt"
file2 = "rocValues.mean.txt"
file3 = "auc.txt"
file4 = "rocMeanSd.txt"

roc.cv = read.table(file1,header=T,sep="\t")
roc.mean = read.table(file2,header=T,sep="\t")
auc.roc = read.table(file3,header=F,sep="\t")
rocs.y.wide = read.table(file4,header=T,sep="\t")

roc.all = rbind(roc.cv,roc.mean)
m.rauc = round(mean(auc.roc$V1),3)
sd.rauc = round(sd(auc.roc$V1),3)

a1 = paste("mean"," (AUC=",m.rauc,"\u00B1",sd.rauc,")", sep="")

# Plotting ROCR curve
plot_1 = "roc_cv.png"
png(plot_1,bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
ggplot() +
geom_line(data = roc.all[roc.all$type2 == "cv",], aes(x=x, y=y, group = type1, color = type2), size = 0.2, alpha = 0.6) +
geom_line(data = roc.all[roc.all$type2 == "mean",], aes(x=x, y=y, color = type2), size = 0.3) +
geom_abline(intercept = 0, linetype = "dashed", color = "black", alpha = 0.3)+
theme_bw()+
scale_colour_manual(
  values= c("grey","blue"),
  breaks = c("cv","mean"),
  labels = c("cross validation",a1))+
geom_ribbon(data = rocs.y.wide, aes(roc.all[roc.all$type2 == "mean",]$x,ymax = mean + sd, ymin = mean - sd), alpha = 0.2, fill = "blue")+
theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
theme(plot.title = element_text(size=11, hjust = 0.5),
  legend.key.size = unit(0.30,"cm"),
  legend.title = element_blank(),
  legend.text = element_text(size = 13),
  legend.position = c(0.65,0.12),
  legend.background = element_rect(fill = "transparent",colour = "NA"),
  axis.text.x = element_text(size=10),
  axis.text.y = element_text(size=10),
  panel.border = element_rect(colour="BLACK",size=0.4),
  axis.title.x = element_text(size=13, vjust = 0.1),
  axis.title.y = element_text(size=13,angle = 90, vjust = 0.1),
  panel.background = element_rect(fill="transparent"),
  plot.background = element_rect(fill = "transparent",colour = NA)
) +
scale_y_continuous(name="True positive rate") +
scale_x_continuous(name = "False positive rate")
#ggtitle(title)
dev.off()


############ PR #####################
file5 = "prValues.cv.txt"
file6 = "prValues.mean.txt"
file7 = "prauc.txt"
file8 = "prMeanSd.txt"

pr.cv = read.table(file5,header=T,sep="\t")
pr.mean = read.table(file6,header=T,sep="\t")
auc.pr = read.table(file7,header=F,sep="\t")
pr.y.wide = read.table(file8,header=T,sep="\t")

pr.all = rbind(pr.cv,pr.mean)
m.prauc = round(mean(auc.pr$V1),3)
sd.prauc = round(sd(auc.pr$V1),3)

a2 = paste("mean"," (AUC=",m.prauc,"\u00B1",sd.prauc,")", sep="")

# Plotting PR curve
plot_2 = "pr_cv.png"
png(plot_2,bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
ggplot() +
geom_line(data = pr.all[pr.all$type2 == "cv",], aes(x=x, y=y, group = type1, color = type2), size = 0.2, alpha = 0.6) +
geom_line(data = pr.all[pr.all$type2 == "mean",], aes(x=x, y=y, color = type2), size = 0.3) +
theme_bw()+
scale_colour_manual(
  values= c("grey","blue"),
  breaks = c("cv","mean"),
  labels = c("cross validation",a2))+
geom_ribbon(data = pr.y.wide, aes(pr.all[pr.all$type2 == "mean",]$x,ymax = mean + sd, ymin = mean - sd), alpha = 0.2, fill = "blue")+
theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
theme(plot.title = element_text(size=11, hjust = 0.5),
  legend.key.size = unit(0.30,"cm"),
  legend.title = element_blank(),
  legend.text = element_text(size = 13),
  legend.position = c(0.35,0.12),
  legend.background = element_rect(fill = "transparent",colour = "NA"),
  axis.text.x = element_text(size=10),
  axis.text.y = element_text(size=10),
  panel.border = element_rect(colour="BLACK",size=0.4),
  axis.title.x = element_text(size=13, vjust = 0.1),
  axis.title.y = element_text(size=13,angle = 90, vjust = 0.1),
  panel.background = element_rect(fill="transparent"),
  plot.background = element_rect(fill = "transparent",colour = NA)
) +
scale_y_continuous(name="Precision") +
scale_x_continuous(name = "Recall")
#ggtitle(title)
dev.off()
