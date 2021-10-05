args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(rstatix)
  library(ggridges)
  library(ggpubr)
})

brain_spec <- args[1]
tis_spec <- args[2]
shared <- args[3]
cnc_known_utrs <- args[4]
outpath <- args[5]
utr_tissue_specificity <- args[6]


bs <- read.table(brain_spec, header=TRUE, sep="\t") %>% dplyr::select(id, cnc, predicted.prob)  %>% mutate(group2 = "Highly brain-specific")
ts <- read.table(tis_spec, header=TRUE, sep="\t") %>% dplyr::select(id, cnc, predicted.prob) %>% mutate(group2 = "Absolute tissue-specific")
s <- read.table(shared, header=TRUE, sep="\t") %>% dplyr::select(id, cnc, predicted.prob)  %>% mutate(group2 = "Shared")

df <- rbind(bs, ts, s)

df$group <- ifelse(df$predicted.prob > 0.60, "Predicted 3-UTRs", "Predicted Non-3-UTRs")

df.pos <- df %>%
  dplyr::filter(group %in% "Predicted 3-UTRs") %>%
  dplyr::select(id, cnc, group, group2)


utr_cnc = read.table(cnc_known_utrs, header=TRUE, sep="\t") %>%
  mutate(group = "Known 3-UTRs") %>%
  plyr::rename(c("three_prime_utr_id" = "id"))

utr_ts <- read.table(utr_tissue_specificity, header=TRUE, sep="\t")
utr_cnc <- dplyr::left_join(utr_cnc, utr_ts, by = "id")


utr.dummy <- rbind(
  utr_cnc %>% dplyr::filter(outcome %in% "Absolute tissue-specific") %>% dplyr::select(-outcome) %>% mutate(group2 = "Absolute tissue-specific"),
  utr_cnc %>% dplyr::filter(outcome %in% "Highly brain-specific") %>% dplyr::select(-outcome) %>% mutate(group2 = "Highly brain-specific"),
  utr_cnc %>% dplyr::filter(outcome %in% "Shared") %>% dplyr::select(-outcome) %>% mutate(group2 = "Shared")
)

df2 = rbind(df.pos, utr.dummy)

df2.samples = table(df2$group, df2$group2)
write.table(df2.samples, file = str_c(outpath, "/nsamples_groups.txt"), sep="\t", row.names=TRUE, col.names=TRUE, quote = FALSE)


# p-values of comparisons
pred.bs <- df2 %>% filter(group2 %in% "Highly brain-specific", group %in% "Predicted 3-UTRs") %>% .$cnc
pred.ts <- df2 %>% filter(group2 %in% "Absolute tissue-specific", group %in% "Predicted 3-UTRs") %>% .$cnc
pred.s <- df2 %>% filter(group2 %in% "Shared", group %in% "Predicted 3-UTRs") %>% .$cnc
known.bs <- utr_cnc %>% dplyr::filter(outcome %in% "Highly brain-specific") %>% .$cnc
known.ts <- utr_cnc %>% dplyr::filter(outcome %in% "Absolute tissue-specific") %>% .$cnc
known.s <- utr_cnc %>% dplyr::filter(outcome %in% "Shared") %>% .$cnc


# wilcoxon and effect size
mw<- wilcox.test(pred.bs, known.bs)
mw.p = mw$p.value %>% format.pval(1)
label.pval = paste(ifelse(str_detect(mw.p, "<"), "p ", "p = "), mw.p, sep="")
es.value <- df2 %>% filter(group2 %in% "Highly brain-specific") %>% rstatix::wilcox_effsize(cnc ~ group) %>% .$effsize
label.es = paste("es = ", round(es.value, 2),sep="")
label.bs = str_c(label.pval, label.es, sep="\n")

rm(label.pval, label.es)

mw<- wilcox.test(pred.ts, known.ts)
mw.p = mw$p.value %>% format.pval(1)
label.pval = paste(ifelse(str_detect(mw.p, "<"), "p ", "p = "), mw.p, sep="")
es.value <- df2 %>% filter(group2 %in% "Absolute tissue-specific") %>% rstatix::wilcox_effsize(cnc ~ group) %>% .$effsize
label.es = paste("es = ", round(es.value, 2),sep="")
label.ts = str_c(label.pval, label.es, sep="\n")

rm(label.pval, label.es)

mw<- wilcox.test(pred.s, known.s)
mw.p = mw$p.value %>% format.pval(1)
label.pval = paste(ifelse(str_detect(mw.p, "<"), "p ", "p = "), mw.p, sep="")
es.value <- df2 %>% filter(group2 %in% "Shared") %>% rstatix::wilcox_effsize(cnc ~ group) %>% .$effsize
label.es = paste("es = ", round(es.value, 2),sep="")
label.s = str_c(label.pval, label.es, sep="\n")



#### ggridges
df2$group2 <- as.factor(df2$group2)
df2$group2 <- fct_relevel(df2$group2, "Shared", "Absolute tissue-specific", "Highly brain-specific")

png(str_c(outpath, "/cnc_density_groups.png"), bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
ggplot(df2, aes(cnc, group2, colour = group, fill=group)) +
  geom_density_ridges(alpha=0.5, scale = 1.2, lwd = 0.6)+
  scale_color_manual(values= c("black", "#00AFBB")) +
  scale_fill_manual(values= c("black", "#00AFBB")) +
  scale_y_discrete(expand = c(0,0), labels = c("Shared", "Absolute\ntissue-specific", "Highly\nbrain-specific"))+
  annotate("text",x=-7.5,y=1.2,size=3.5,label= label.s) +
  annotate("text",x=-7.5,y=2.2,size=3.5,label= label.ts) +
  annotate("text",x=-7.5,y=3.2,size=3.5,label= label.bs) +
  scale_x_continuous(expand = c(0,0), name = "Mean CNC score")+
  coord_cartesian(clip = "off") +
  theme_ridges(center_axis_labels = TRUE)  +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.30,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        legend.position= c(-0.02,0.95),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size = 11, face = "bold"),
        panel.border = element_rect(colour="BLACK",size=0.4),
        #panel.border = element_blank(),
        axis.title.x = element_text(size=11),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )
dev.off()
