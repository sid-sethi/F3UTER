args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(rstatix)
  library(coin)
  library(ggsignif)
  library(ggridges)
})


res_dir <- args[1]
outpath = args[2]
utr_tissue_specificity <- args[3]


xvalue_all <- c()
xvalue_groups <- c()
yvalue_groups <- c()

if(str_detect(res_dir, "Attract")){
  xvalue_all <- 150
  xvalue_groups <- 300
  yvalue_groups <- 65536
}else{
  xvalue_all <- 50
  xvalue_groups <- 80
  yvalue_groups <- 4096
}




all.pos <- read.table(str_c(res_dir, "/all.pos_seq_sum.txt"), header=TRUE) %>%
  mutate(group = "Predicted 3-UTRs")
all.neg <- read.table(str_c(res_dir, "/all.neg_seq_sum.txt"), header=TRUE) %>%
  mutate(group = "Predicted non-3-UTRs")

### known 3'UTRs
utr_rbp = read.table(str_c(res_dir, "/known_utrs_seq_sum.txt"), header=TRUE, sep="\t") %>%
  mutate(group = "Known 3-UTRs")
utr_ts <- read.table(utr_tissue_specificity, header=TRUE, sep="\t")
utr_rbp <- dplyr::left_join(utr_rbp, utr_ts, by = c("sequence_name" = "id"))


utr_all <- utr_rbp %>% dplyr::filter(! outcome %in% "Not expressed")



### negative control
cntrl = read.table(str_c(res_dir, "/all.pos.shuf_seq_sum.txt"), header=TRUE, sep="\t") %>%
  mutate(group = "Negative control")

df <- bind_rows(all.pos, all.neg, utr_all, cntrl)

df$group <- fct_relevel(df$group, "Negative control", "Known 3-UTRs", "Predicted non-3-UTRs", "Predicted 3-UTRs")

df.zero <- df %>% dplyr::filter(enrich_score == 0)
df <- df %>% dplyr::filter(enrich_score > 0)


no_enrich <- table(df.zero$group) %>%
  as.data.frame() %>%
  mutate(
    total = c(nrow(cntrl), nrow(utr_all), nrow(all.neg), nrow(all.pos)),
    perc = round((Freq/total)*100, 2)
  )

write.table(no_enrich, file = str_c(outpath, "/all.no_enrich.txt"), sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)

df.samples <- table(df$group) %>% as.data.frame()
write.table(df.samples, file = str_c(outpath, "/nsamples.enrich.txt"), sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)

df.zero.samples <- table(df.zero$group) %>% as.data.frame()
write.table(df.zero.samples, file = str_c(outpath, "/nsamples.no_enrich.txt"), sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)



png(str_c(outpath, "/all.noenrich.png"),bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)

ggplot(data=no_enrich, aes(x= Var1, y=perc, fill = Var1)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8, alpha=0.6) +
  theme_bw() +
  scale_fill_manual(values = c("#8B008B", "black", "#E7B800", "#00AFBB")) +
  guides(fill=FALSE) +
  theme(plot.title = element_text(size=11, hjust=0.5),
        legend.key.size = unit(0.40,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position= c(0.66,0.95),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=11, angle=40, hjust=1),
        axis.text.y = element_text(size = 10),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=11,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_y_continuous(name = "% regions with no RBP enrichment")

dev.off()



pos <- all.pos$enrich_score
neg <- all.neg$enrich_score
known <- utr_all$enrich_score
neg_cntrl <- cntrl$enrich_score


mw<- wilcox.test(pos, known)
mw.p = format.pval(mw$p.value, 1)
label.pval = paste(ifelse(str_detect(mw.p, "<"), "p ", "p = "), mw.p, sep="")
es <- df %>%  rstatix::wilcox_effsize(enrich_score ~ group, comparisons = list(c("Predicted 3-UTRs", "Known 3-UTRs")), ref.group = "Predicted 3-UTRs") %>% .$effsize
label.es = paste("es = ", round(es, 2),sep="")
label.pk = str_c("x: ", label.pval, "; ", label.es)
label.pk2 = str_c(label.pval, label.es, sep="; ")


mw<- wilcox.test(pos, neg)
mw.p = format.pval(mw$p.value, 1)
label.pval = paste(ifelse(str_detect(mw.p, "<"), "p ", "p = "), mw.p, sep="")
es <- df %>%  rstatix::wilcox_effsize(enrich_score ~ group, comparisons = list(c("Predicted 3-UTRs", "Predicted non-3-UTRs")), ref.group = "Predicted 3-UTRs") %>% .$effsize
label.es = paste("es = ", round(es, 2),sep="")
label.pn = str_c("y: ", label.pval, "; ", label.es)
label.pn2 = str_c(label.pval, label.es, sep="; ")


mw<- wilcox.test(pos, neg_cntrl)
mw.p = format.pval(mw$p.value, 1)
label.pval = paste(ifelse(str_detect(mw.p, "<"), "p ", "p = "), mw.p, sep="")
es <- df %>%  rstatix::wilcox_effsize(enrich_score ~ group, comparisons = list(c("Predicted 3-UTRs", "Negative control")), ref.group = "Predicted 3-UTRs") %>% .$effsize
label.es = paste("es = ", round(es, 2),sep="")
label.pc = str_c("z: ", label.pval, "; ", label.es)
label.pc2 = str_c(label.pval, label.es, sep="; ")


mw<- wilcox.test(neg, known)
mw.p = format.pval(mw$p.value, 1)
label.pval = paste(ifelse(str_detect(mw.p, "<"), "p ", "p = "), mw.p, sep="")
es <- df %>%  rstatix::wilcox_effsize(enrich_score ~ group, comparisons = list(c("Predicted non-3-UTRs", "Known 3-UTRs")), ref.group = "Predicted non-3-UTRs") %>% .$effsize
label.es = paste("es = ", round(es, 2),sep="")
label.nk2 = str_c(label.pval, label.es, sep="; ")


mw<- wilcox.test(neg, neg_cntrl)
mw.p = format.pval(mw$p.value, 1)
label.pval = paste(ifelse(str_detect(mw.p, "<"), "p ", "p = "), mw.p, sep="")
es <- df %>%  rstatix::wilcox_effsize(enrich_score ~ group, comparisons = list(c("Predicted non-3-UTRs", "Negative control")), ref.group = "Predicted non-3-UTRs") %>% .$effsize
label.es = paste("es = ", round(es, 2),sep="")
label.nc2 = str_c(label.pval, label.es, sep="; ")



png(str_c(outpath, "/all.density.png"),bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
ggdensity(df, x = "enrich_score", y = "..density..",
          #add = "median",
          rug = TRUE,
          color = "group",
          fill = "group",
          #palette = c("#00AFBB", "#E7B800")
          palette = c("#8B008B", "black", "#E7B800", "#00AFBB")
)+
  annotate("text",x=xvalue_all,y=0.25,label= label.pk) +
  annotate("text",x=xvalue_all,y=0.225,label= label.pn) +
  annotate("text",x=xvalue_all,y=0.20,label= label.pc) +
  scale_x_continuous(name = "Number of RBP motif binding per 100 nt", trans="log2", breaks = c(0.5, 8, 128))+
  scale_y_continuous(name = "Density")+
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.40,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        legend.position= c(0.77,0.90),
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



png(str_c(outpath, "/all.violin.png"),bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
ggplot(df, aes(x=group, y=enrich_score, fill=group)) +
  geom_violin(trim=TRUE, alpha=0.5)+
  geom_boxplot(width=0.2, lwd = 0.55, alpha = 0.6, outlier.shape=19, outlier.size=0.7, outlier.alpha=0.1) +
  theme_bw() +
  guides(fill=FALSE) +
  scale_fill_manual(values = c("#8B008B", "black", "#E7B800", "#00AFBB"))+
  scale_y_continuous(name = "Number of RBP motif binding per 100 nt", trans="log2", breaks = c(0.5, 8, 128))+
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.40,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        legend.position= c(0.77,0.95),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=10, angle=40, hjust=1),
        axis.text.y = element_text(size = 9),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        #panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=10,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  geom_signif(
    comparisons = list(c("Predicted 3-UTRs", "Predicted non-3-UTRs")),
    map_signif_level=TRUE, annotation = label.pn2,
    size=0.4, textsize=3, tip_length=0.03, margin_top =0.05,
    color = "black", fontface = 2
  ) +
  geom_signif(
    comparisons = list(c("Predicted 3-UTRs", "Known 3-UTRs")),
    map_signif_level=TRUE, annotation = label.pk2,
    size=0.4, textsize=3, tip_length=0.03, margin_top =0.15,
    color = "black", fontface = 2
  ) +
  geom_signif(
    comparisons = list(c("Predicted 3-UTRs", "Negative control")),
    map_signif_level=TRUE, annotation = label.pc2,
    size=0.4, textsize=3, tip_length=0.03, margin_top =0.25,
    color = "black", fontface = 2
  )
dev.off()




########################################
########################################

bs.p <- read.table(str_c(res_dir, "/brainSpec.pos_seq_sum.txt"), header=TRUE, sep="\t") %>% mutate(group = "Predicted 3-UTRs", group2 = "Highly brain-specific")
bs.n <- read.table(str_c(res_dir, "/brainSpec.neg_seq_sum.txt"), header=TRUE, sep="\t") %>% mutate(group = "Predicted non-3-UTRs", group2 = "Highly brain-specific")
ts.p <- read.table(str_c(res_dir, "/tisSpec.pos_seq_sum.txt"), header=TRUE, sep="\t") %>% mutate(group = "Predicted 3-UTRs", group2 = "Absolute tissue-specific")
ts.n <- read.table(str_c(res_dir, "/tisSpec.neg_seq_sum.txt"), header=TRUE, sep="\t") %>% mutate(group = "Predicted non-3-UTRs", group2 = "Absolute tissue-specific")
s.p <- read.table(str_c(res_dir, "/shared.pos_seq_sum.txt"), header=TRUE, sep="\t") %>% mutate(group = "Predicted 3-UTRs", group2 = "Shared")
s.n <- read.table(str_c(res_dir, "/shared.neg_seq_sum.txt"), header=TRUE, sep="\t") %>% mutate(group = "Predicted non-3-UTRs", group2 = "Shared")


bs.c <- read.table(str_c(res_dir, "/brainSpec.pos.shuf_seq_sum.txt"), header=TRUE, sep="\t") %>% mutate(group = "Negative control", group2 = "Highly brain-specific")
ts.c <- read.table(str_c(res_dir, "/tisSpec.pos.shuf_seq_sum.txt"), header=TRUE, sep="\t") %>% mutate(group = "Negative control", group2 = "Absolute tissue-specific")
s.c <- read.table(str_c(res_dir, "/shared.pos.shuf_seq_sum.txt"), header=TRUE, sep="\t") %>% mutate(group = "Negative control", group2 = "Shared")


utr.dummy <- rbind(
  utr_rbp %>% dplyr::filter(outcome %in% "Absolute tissue-specific") %>% mutate(group2 = "Absolute tissue-specific"),
  utr_rbp %>% dplyr::filter(outcome %in% "Highly brain-specific") %>% mutate(group2 = "Highly brain-specific"),
  utr_rbp %>% dplyr::filter(outcome %in% "Shared") %>% mutate(group2 = "Shared")
)

df2 <- bind_rows(ts.p, ts.n, bs.p, bs.n, s.p, s.n, utr.dummy, bs.c, ts.c, s.c)


df2$group2 <- as.factor(df2$group2)
df2$group2 <- fct_relevel(df2$group2, "Shared", "Absolute tissue-specific", "Highly brain-specific")
df2$group <- fct_relevel(df2$group, "Negative control", "Known 3-UTRs", "Predicted non-3-UTRs", "Predicted 3-UTRs")

df2.zero <- df2 %>% dplyr::filter(enrich_score == 0)
df2 <- df2 %>% dplyr::filter(enrich_score > 0)


no_enrich2 <- df2.zero %>%
  group_by(group2, group) %>%
  summarise(Freq = n(), .groups = "keep") %>%
  as.data.frame() %>%
  mutate(
    total = c(nrow(s.c), nrow(utr_rbp %>% dplyr::filter(outcome %in% "Shared")), nrow(s.n), nrow(s.p),
              nrow(ts.c), nrow(utr_rbp %>% dplyr::filter(outcome %in% "Absolute tissue-specific")), nrow(ts.n), nrow(ts.p),
              nrow(bs.c), nrow(utr_rbp %>% dplyr::filter(outcome %in% "Highly brain-specific")), nrow(bs.n), nrow(bs.p)
              ),
    perc = round((Freq/total)*100, 2)
  )


write.table(no_enrich2, file = str_c(outpath, "/groups.no_enrich.txt"), sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)

df2.samples = table(df2$group, df2$group2)
write.table(df2.samples, file = str_c(outpath, "/nsamples_groups.enrich.txt"), sep="\t", row.names=TRUE, col.names=TRUE, quote = FALSE)

df2.zero.samples = table(df2.zero$group, df2.zero$group2)
write.table(df2.zero.samples, file = str_c(outpath, "/nsamples_groups.no_enrich.txt"), sep="\t", row.names=TRUE, col.names=TRUE, quote = FALSE)

########################################################################
#######################################################################
############### gg ridges #############################################
#######################################################################



# p-values of comparisons
pred.bs.pos <- bs.p$enrich_score
pred.ts.pos <- ts.p$enrich_score
pred.s.pos <- s.p$enrich_score
pred.bs.neg <- bs.n$enrich_score
pred.ts.neg <- ts.n$enrich_score
pred.s.neg <- s.n$enrich_score
#known <- utr$enrich_score
known.bs <- utr_rbp %>% dplyr::filter(outcome %in% "Highly brain-specific") %>% .$enrich_score
known.ts <- utr_rbp %>% dplyr::filter(outcome %in% "Absolute tissue-specific") %>% .$enrich_score
known.s <- utr_rbp %>% dplyr::filter(outcome %in% "Shared") %>% .$enrich_score
cntrl.bs <- bs.c$enrich_score
cntrl.ts <- ts.c$enrich_score
cntrl.s <- s.c$enrich_score


# wilcoxon and effect size

mw<- wilcox.test(pred.bs.pos, pred.bs.neg)
mw.p = mw$p.value %>% format.pval(1)
label.pval = paste(ifelse(str_detect(mw.p, "<"), "p ", "p = "), mw.p, sep="")
es.value <- df2 %>% filter(group2 %in% "Highly brain-specific") %>% rstatix::wilcox_effsize(enrich_score ~ group, comparisons = list(c("Predicted 3-UTRs", "Predicted non-3-UTRs")), ref.group = "Predicted 3-UTRs") %>% .$effsize
label.es = paste("es = ", round(es.value, 2),sep="")
label.bs.pn = str_c("x: ", label.pval, "; ", label.es)
label.bs.pn2 = str_c(label.pval, label.es, sep="; ")

mw<- wilcox.test(pred.bs.pos, known.bs)
mw.p = mw$p.value %>% format.pval(1)
label.pval = paste(ifelse(str_detect(mw.p, "<"), "p ", "p = "), mw.p, sep="")
es.value <- df2 %>% filter(group2 %in% "Highly brain-specific") %>% rstatix::wilcox_effsize(enrich_score ~ group, comparisons = list(c("Predicted 3-UTRs", "Known 3-UTRs")), ref.group = "Predicted 3-UTRs") %>% .$effsize
label.es = paste("es = ", round(es.value, 2),sep="")
label.bs.pk = str_c("y: ", label.pval, "; ", label.es)
label.bs.pk2 = str_c(label.pval, label.es, sep="; ")

mw<- wilcox.test(pred.bs.pos, cntrl.bs)
mw.p = mw$p.value %>% format.pval(1)
label.pval = paste(ifelse(str_detect(mw.p, "<"), "p ", "p = "), mw.p, sep="")
es.value <- df2 %>% filter(group2 %in% "Highly brain-specific") %>% rstatix::wilcox_effsize(enrich_score ~ group, comparisons = list(c("Predicted 3-UTRs", "Negative control")), ref.group = "Predicted 3-UTRs") %>% .$effsize
label.es = paste("es = ", round(es.value, 2),sep="")
label.bs.pc = str_c("z: ", label.pval, "; ", label.es)
label.bs.pc2 = str_c(label.pval, label.es, sep="; ")



mw<- wilcox.test(pred.ts.pos, pred.ts.neg)
mw.p = mw$p.value %>% format.pval(1)
label.pval = paste(ifelse(str_detect(mw.p, "<"), "p ", "p = "), mw.p, sep="")
es.value <- df2 %>% filter(group2 %in% "Absolute tissue-specific") %>% rstatix::wilcox_effsize(enrich_score ~ group, comparisons = list(c("Predicted 3-UTRs", "Predicted non-3-UTRs")), ref.group = "Predicted 3-UTRs") %>% .$effsize
label.es = paste("es = ", round(es.value, 2),sep="")
label.ts.pn = str_c("x: ", label.pval, "; ", label.es)
label.ts.pn2 = str_c(label.pval, label.es, sep="; ")

mw<- wilcox.test(pred.ts.pos, known.ts)
mw.p = mw$p.value %>% format.pval(1)
label.pval = paste(ifelse(str_detect(mw.p, "<"), "p ", "p = "), mw.p, sep="")
es.value <- df2 %>% filter(group2 %in% "Absolute tissue-specific") %>% rstatix::wilcox_effsize(enrich_score ~ group, comparisons = list(c("Predicted 3-UTRs", "Known 3-UTRs")), ref.group = "Predicted 3-UTRs") %>% .$effsize
label.es = paste("es = ", round(es.value, 2),sep="")
label.ts.pk = str_c("y: ", label.pval, "; ", label.es)
label.ts.pk2 = str_c(label.pval, label.es, sep="; ")

mw<- wilcox.test(pred.ts.pos, cntrl.ts)
mw.p = mw$p.value %>% format.pval(1)
label.pval = paste(ifelse(str_detect(mw.p, "<"), "p ", "p = "), mw.p, sep="")
es.value <- df2 %>% filter(group2 %in% "Absolute tissue-specific") %>% rstatix::wilcox_effsize(enrich_score ~ group, comparisons = list(c("Predicted 3-UTRs", "Negative control")), ref.group = "Predicted 3-UTRs") %>% .$effsize
label.es = paste("es = ", round(es.value, 2),sep="")
label.ts.pc = str_c("z: ", label.pval, "; ", label.es)
label.ts.pc2 = str_c(label.pval, label.es, sep="; ")



mw<- wilcox.test(pred.s.pos, pred.s.neg)
mw.p = mw$p.value %>% format.pval(1)
label.pval = paste(ifelse(str_detect(mw.p, "<"), "p ", "p = "), mw.p, sep="")
es.value <- df2 %>% filter(group2 %in% "Shared") %>% rstatix::wilcox_effsize(enrich_score ~ group, comparisons = list(c("Predicted 3-UTRs", "Predicted non-3-UTRs")), ref.group = "Predicted 3-UTRs") %>% .$effsize
label.es = paste("es = ", round(es.value, 2),sep="")
label.s.pn = str_c("x: ", label.pval, "; ", label.es)
label.s.pn2 = str_c(label.pval, label.es, sep="; ")

mw<- wilcox.test(pred.s.pos, known.s)
mw.p = mw$p.value %>% format.pval(1)
label.pval = paste(ifelse(str_detect(mw.p, "<"), "p ", "p = "), mw.p, sep="")
es.value <- df2 %>% filter(group2 %in% "Shared") %>% rstatix::wilcox_effsize(enrich_score ~ group, comparisons = list(c("Predicted 3-UTRs", "Known 3-UTRs")), ref.group = "Predicted 3-UTRs") %>% .$effsize
label.es = paste("es = ", round(es.value, 2),sep="")
label.s.pk = str_c("y: ", label.pval, "; ", label.es)
label.s.pk2 = str_c(label.pval, label.es, sep="; ")

mw<- wilcox.test(pred.s.pos, cntrl.s)
mw.p = mw$p.value %>% format.pval(1)
label.pval = paste(ifelse(str_detect(mw.p, "<"), "p ", "p = "), mw.p, sep="")
es.value <- df2 %>% filter(group2 %in% "Shared") %>% rstatix::wilcox_effsize(enrich_score ~ group, comparisons = list(c("Predicted 3-UTRs", "Negative control")), ref.group = "Predicted 3-UTRs") %>% .$effsize
label.es = paste("es = ", round(es.value, 2),sep="")
label.s.pc = str_c("z: ", label.pval, "; ", label.es)
label.s.pc2 = str_c(label.pval, label.es, sep="; ")




png(str_c(outpath, "/groups_density.png"),bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
ggplot(df2, aes(enrich_score, group2, colour = group, fill=group)) +
  geom_density_ridges(alpha=0.5, scale = 1.2, lwd = 0.6)+
  scale_x_continuous(expand = c(0,0), trans="log2", name = "Number of RBP motif binding per 100 nt", breaks = c(0.5, 8, 128)) +
  scale_color_manual(values= c("#8B008B", "black", "#E7B800", "#00AFBB")) +
  scale_fill_manual(values= c("#8B008B", "black", "#E7B800", "#00AFBB")) +
  scale_y_discrete(expand = c(0,0), labels = c("Shared", "Absolute\ntissue-specific", "Highly\nbrain-specific"))+
  annotate("text",x=xvalue_groups,y=1.9,size=3.1,label= label.s.pn) +
  annotate("text",x=xvalue_groups,y=1.75,size=3.1,label= label.s.pk) +
  annotate("text",x=xvalue_groups,y=1.60,size=3.1,label= label.s.pc) +
  annotate("text",x=xvalue_groups,y=2.9,size=3.1,label= label.ts.pn) +
  annotate("text",x=xvalue_groups,y=2.75,size=3.1,label= label.ts.pk) +
  annotate("text",x=xvalue_groups,y=2.60,size=3.1,label= label.ts.pc) +
  annotate("text",x=xvalue_groups,y=4.1,size=3.1,label= label.bs.pn) +
  annotate("text",x=xvalue_groups,y=3.95,size=3.1,label= label.bs.pk) +
  annotate("text",x=xvalue_groups,y=3.80,size=3.1,label= label.bs.pc) +
  coord_cartesian(clip = "off") +
  theme_ridges(center_axis_labels = TRUE)  +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.30,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        legend.position= c(-0.4,0.93),
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



df2$group2 <- fct_relevel(df2$group2, "Highly brain-specific", "Absolute tissue-specific", "Shared")


myplot <- ggplot(df2, aes(x=group, y=enrich_score, fill=group)) +
  geom_violin(trim=TRUE, alpha=0.5)+
  geom_boxplot(width=0.2, lwd = 0.55, alpha = 0.6, outlier.shape=19, outlier.size=0.7, outlier.alpha=0.1) +
  theme_bw() +
  facet_wrap(~group2, nrow=1) +
  guides(fill=FALSE) +
  scale_fill_manual(values = c("#8B008B", "black", "#E7B800", "#00AFBB"))+
  scale_y_continuous(name = "Number of RBP motif binding per 100 nt", trans = "log2", breaks = c(0.5, 8, 128))+
  coord_cartesian(ylim = c(NA, yvalue_groups)) +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.40,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        legend.position= c(0.77,0.95),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=12, angle=40, hjust=1),
        axis.text.y = element_text(size = 10),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        #panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA),
        strip.text.x = element_text(size=12, face="bold"),
        strip.background = element_rect(fill="white")
  ) +
  geom_signif(
    comparisons = list(c("Predicted 3-UTRs", "Predicted non-3-UTRs"), c("Predicted 3-UTRs", "Known 3-UTRs"), c("Predicted 3-UTRs", "Negative control")),
    annotation = c("foo"),
    size=0.4, textsize=2.5, tip_length=0.03, margin_top =0.05, step_increase = 0.15,
    color = "black"
    #fontface = 2
  )



myplot2 <- ggplot_build(myplot)
myplot2$data[[3]]$annotation <- c(
              rep(label.bs.pn2, 3), rep(label.bs.pk2, 3), rep(label.bs.pc2, 3),
              rep(label.ts.pn2, 3), rep(label.ts.pk2, 3), rep(label.ts.pc2, 3),
              rep(label.s.pn2, 3), rep(label.s.pk2, 3), rep(label.s.pc2, 3)
            )
myplot3 <- ggplot_gtable(myplot2)

png(str_c(outpath, "/groups_violin.png"),bg="transparent",units="in",height = 4.25, width= 8.25 ,res=600)
plot(myplot3)
dev.off()






### no enrichment


png(str_c(outpath, "/groups.noenrich.png"),bg="transparent",units="in",width = 8.25, height= 3.75 ,res=600)

ggplot(data=no_enrich2, aes(x= group, y=perc, fill = group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8, alpha=0.6) +
  theme_bw() +
  scale_fill_manual(values = c("#8B008B", "black", "#E7B800", "#00AFBB")) +
  facet_wrap(~group2, nrow=1) +
  guides(fill=FALSE) +
  theme(plot.title = element_text(size=11, hjust=0.5),
        legend.key.size = unit(0.40,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position= c(0.66,0.95),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=11, angle=40, hjust=1),
        axis.text.y = element_text(size = 10),
        panel.border = element_rect(colour="BLACK",size=0.4),
        #panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=11,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA),
        strip.text.x = element_text(size=11, face="bold"),
        strip.background = element_rect(fill="white")
  )+
  #theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_y_continuous(name = "% regions with no RBP enrichment")

dev.off()
