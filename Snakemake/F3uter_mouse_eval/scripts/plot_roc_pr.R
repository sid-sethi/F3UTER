args<-commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
})

roc <- args[1]
pr <- args[2]
auc <- args[3]
outpath <- args[4]
tissue <- args[5]


##### ROCR ##############

roc.mean = read.table(roc, header=TRUE, sep="\t")
auc = read.table(auc, header=TRUE, sep="\t")


a1 = str_c("AUC = ", round(auc$ROCAUC,2))

# Plotting ROCR curve
png(str_c(outpath, "/", tissue, "_roc.png"),bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
ggplot() +
geom_line(data = roc.mean, aes(x=x, y=y, color = type1), size = 0.9) +
geom_abline(intercept = 0, linetype = "dashed", color = "black", alpha = 0.3)+
theme_bw()+
scale_colour_manual(
  values= c("blue"),
  labels = c(a1)) +
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
dev.off()


############ PR #####################

pr.mean = read.table(pr, header=TRUE, sep="\t")

a2 = paste("AUC=", round(auc$PRAUC,2))

# Plotting PR curve
png(str_c(outpath, "/", tissue, "_pr.png"),bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
ggplot() +
geom_line(data = pr.mean, aes(x=x, y=y, color = type1), size = 0.9) +
theme_bw()+
scale_colour_manual(
  values= c("blue"),
  labels = c(a2)) +
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
dev.off()
