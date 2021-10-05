args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(rstatix)
  library(coin)
  library(ggpubr)
})

cnc_all <- args[1]
cnc_known_utrs <- args[2]
outpath <- args[3]
utr_tissue_specificity <- args[4]


df <- read.table(cnc_all, header=TRUE, sep="\t")
df$group <- ifelse(df$predicted.prob > 0.60, "Predicted 3-UTRs", "Predicted Non-3-UTRs")

df.pos <- df %>%
  dplyr::filter(group %in% "Predicted 3-UTRs") %>%
  dplyr::select(id, cnc, group)


utr_cnc = read.table(cnc_known_utrs, header=TRUE, sep="\t") %>%
  mutate(group = "Known 3-UTRs") %>%
  plyr::rename(c("three_prime_utr_id" = "id"))

utr_ts <- read.table(utr_tissue_specificity, header=TRUE, sep="\t")
utr_cnc <- dplyr::left_join(utr_cnc, utr_ts, by = "id")

utr_all <- utr_cnc %>%
  dplyr::filter(! outcome %in% "Not expressed") %>%
  dplyr::select(-outcome)


df2 = rbind(df.pos, utr_all)

df2.samples = table(df2$group) %>% as.data.frame()
write.table(df2.samples, file = str_c(outpath, "/nsamples_all_predictions.txt"), sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)


pos <- df2 %>% filter(group %in% "Predicted 3-UTRs") %>% .$cnc
known <- df2 %>% filter(group %in% "Known 3-UTRs") %>% .$cnc


# man-whitney test
mw<- wilcox.test(pos, known)
mw.p = format.pval(mw$p.value, 1)
label.pval = paste(ifelse(str_detect(mw.p, "<"), "p ", "p = "), mw.p, sep="")
es <- df2 %>% rstatix::wilcox_effsize(cnc ~ group) %>% .$effsize
label.es = paste("es = ", round(es, 2),sep="")
label = str_c(label.pval, label.es, sep="\n")


png(str_c(outpath, "/cnc_density_all_predictions.png"), bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
ggdensity(df2, x = "cnc", y = "..density..",
          #add = "mean",
          rug = TRUE,
          color = "group", fill = "group",
          palette = c("black", "#00AFBB")
)+
  annotate("text", x = -7, y = 0.15, label = label) +
  scale_x_continuous(name = "Mean CNC score")+
  scale_y_continuous(name = "Density")+
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.40,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        legend.position= c(0.23,0.95),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size = 9),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        #panel.border = element_blank(),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )
dev.off()
