args <- commandArgs(TRUE) # 1ST = RBP scanning data path
library(tidyverse)
library(ggpubr)
library(rstatix)

print(args[1])
setwd(args[1])

all.pos <- read.table("all.pos_seq_sum.txt", header=TRUE) %>%
  mutate(group = "Predicted 3-UTRs")
all.neg <- read.table("all.neg_seq_sum.txt", header=TRUE) %>%
  mutate(group = "Predicted non-3-UTRs")


### known 3'UTRs
utr = read.table("/RBP_scan_3/10e-4/utr_seq_sum.txt", header=TRUE, sep="\t") %>%
  mutate(group = "Known 3-UTRs")

df <- bind_rows(all.pos, all.neg, utr)

ks <- ks.test(all.pos$enrich_sum, all.neg$enrich_sum, alternative = "l")
p = ks$p.value %>% format.pval(3)
print(p)

ks <- ks.test(all.pos$enrich_sum, utr$enrich_sum, alternative = "l")
p = ks$p.value %>% format.pval(3)
print(p)

pos <- all.pos$enrich_sum
neg <- all.neg$enrich_sum
known <- utr$enrich_sum

# man-whitney test
mw<- wilcox.test(pos, known)
mw.p = mw$p.value
label.pval = paste("p = ",format.pval(mw.p,3),sep="")
es <- df %>%  rstatix::wilcox_effsize(enrich_sum ~ group, comparisons = list(c("Predicted 3-UTRs", "Known 3-UTRs")), ref.group = "Predicted 3-UTRs") %>% .$effsize
label.es = paste("es = ", round(es, 3),sep="")


mw<- wilcox.test(pos, neg)
mw.p = mw$p.value
label.pval = paste("p = ",format.pval(mw.p,3),sep="")
es <- df %>%  rstatix::wilcox_effsize(enrich_sum ~ group, comparisons = list(c("Predicted 3-UTRs", "Predicted non-3-UTRs")), ref.group = "Predicted 3-UTRs") %>% .$effsize
label.es = paste("es = ", round(es, 3),sep="")


png("all.density.png",bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
ggdensity(df, x = "enrich_sum", y = "..density..",
          #add = "median",
          rug = TRUE,
          color = "group",
          fill = "group",
          #palette = c("#00AFBB", "#E7B800")
          palette = c("black", "#00AFBB", "#E7B800")
)+
  annotate("text",x=150,y=0.25,label=c("x: p < 2e-16; es = 0.28")) +
  annotate("text",x=150,y=0.225,label=c("y: p < 2e-16; es = 0.172")) +
  scale_x_continuous(name = "Number of RBP motif binding per 100 nt", trans="log2")+
  scale_y_continuous(name = "Density")+
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.40,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        legend.position= c(0.77,0.95),
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


########################################
##### ER categories ####################
########################################

# load FIMO result files

bs.p <- read.table("/brain_spec_pos_data_minimal.txt", header=TRUE, sep="\t") # brain-specific positive
bs.n <- read.table("/brain_spec_neg_data_minimal.txt", header=TRUE, sep="\t") # brain-specific negative
ts.p <- read.table("/tis_spec_pos_data_minimal.txt", header=TRUE, sep="\t") # tissue-specific positive
ts.n <- read.table("/tis_spec_neg_data_minimal.txt", header=TRUE, sep="\t") # tissue-specific negative
com.p <- read.table("/common_pos_data_minimal.txt", header=TRUE, sep="\t") # shared positive
com.n <- read.table("/common_neg_data_minimal.txt", header=TRUE, sep="\t") # shared negative

bs.pos <- semi_join(df, bs.p, by = c("sequence_name" = "id")) %>% mutate(group2 = "Highly brain-specific") %>% plyr::rename(c("group" = "group1"))
bs.neg <- semi_join(df, bs.n, by = c("sequence_name" = "id")) %>% mutate(group2 = "Highly brain-specific") %>% plyr::rename(c("group" = "group1"))
ts.pos <- semi_join(df, ts.p, by = c("sequence_name" = "id")) %>% mutate(group2 = "Absolute tissue-specific") %>% plyr::rename(c("group" = "group1"))
ts.neg <- semi_join(df, ts.n, by = c("sequence_name" = "id")) %>% mutate(group2 = "Absolute tissue-specific") %>% plyr::rename(c("group" = "group1"))
com.pos <- semi_join(df, com.p, by = c("sequence_name" = "id")) %>% mutate(group2 = "Shared") %>% plyr::rename(c("group" = "group1"))
com.neg <- semi_join(df, com.n, by = c("sequence_name" = "id")) %>% mutate(group2 = "Shared") %>% plyr::rename(c("group" = "group1"))


df2 <- bind_rows(ts.pos, ts.neg, bs.pos, bs.neg, com.pos, com.neg)

utr = read.table("/RBP_scan_3/10e-4/utr_seq_sum.txt", header=TRUE, sep="\t") %>%
  mutate(group1 = "Known 3-UTRs")
utr.dummy <- rbind(
  utr %>% mutate(group2 = "Absolute tissue-specific"),
  utr %>% mutate(group2 = "Highly brain-specific"),
  utr %>% mutate(group2 = "Shared")
)

df3 <- rbind(df2, utr.dummy)
df3$group2 <- as.factor(df3$group2)
df3$group2 <- fct_relevel(df3$group2, "Shared", "Absolute tissue-specific", "Highly brain-specific")



# p-values of comparisons
pred.bs.pos <- df3 %>% filter(group2 %in% "Highly brain-specific", group1 %in% "Predicted 3-UTRs") %>% .$enrich_sum
pred.ts.pos <- df3 %>% filter(group2 %in% "Absolute tissue-specific", group1 %in% "Predicted 3-UTRs") %>% .$enrich_sum
pred.s.pos <- df3 %>% filter(group2 %in% "Shared", group1 %in% "Predicted 3-UTRs") %>% .$enrich_sum
pred.bs.neg <- df3 %>% filter(group2 %in% "Highly brain-specific", group1 %in% "Predicted non-3-UTRs") %>% .$enrich_sum
pred.ts.neg <- df3 %>% filter(group2 %in% "Absolute tissue-specific", group1 %in% "Predicted non-3-UTRs") %>% .$enrich_sum
pred.s.neg <- df3 %>% filter(group2 %in% "Shared", group1 %in% "Predicted non-3-UTRs") %>% .$enrich_sum
known <- utr$enrich_sum


# ks test
ks <- ks.test(pred.bs.pos, known, alternative = "l")
p = ks$p.value %>% format.pval(3)
print(p)

ks <- ks.test(pred.bs.pos, pred.bs.neg, alternative = "l")
p = ks$p.value %>% format.pval(3)
print(p)

ks <- ks.test(pred.ts.pos, known, alternative = "l")
p = ks$p.value %>% format.pval(3)
print(p)

ks <- ks.test(pred.ts.pos, pred.ts.neg, alternative = "l")
p = ks$p.value %>% format.pval(3)
print(p)

ks <- ks.test(pred.s.pos, known, alternative = "l")
p = ks$p.value %>% format.pval(3)
print(p)

ks <- ks.test(pred.s.pos, pred.s.neg, alternative = "l")
p = ks$p.value %>% format.pval(3)
print(p)


# wilcoxon and effect size
mw<- wilcox.test(pred.bs.pos, known)
mw.p = mw$p.value %>% format.pval(3)
print(mw.p)
es.value <- df3 %>% filter(group2 %in% "Highly brain-specific") %>% rstatix::wilcox_effsize(enrich_sum ~ group1, comparisons = list(c("Predicted 3-UTRs", "Known 3-UTRs")), ref.group = "Predicted 3-UTRs") %>% .$effsize
print(es.value)

mw<- wilcox.test(pred.bs.pos, pred.bs.neg)
mw.p = mw$p.value %>% format.pval(3)
print(mw.p)
es.value <- df3 %>% filter(group2 %in% "Highly brain-specific") %>% rstatix::wilcox_effsize(enrich_sum ~ group1, comparisons = list(c("Predicted 3-UTRs", "Predicted non-3-UTRs")), ref.group = "Predicted 3-UTRs") %>% .$effsize
print(es.value)

mw<- wilcox.test(pred.ts.pos, known)
mw.p = mw$p.value %>% format.pval(3)
print(mw.p)
es.value <- df3 %>% filter(group2 %in% "Absolute tissue-specific") %>% rstatix::wilcox_effsize(enrich_sum ~ group1, comparisons = list(c("Predicted 3-UTRs", "Known 3-UTRs")), ref.group = "Predicted 3-UTRs") %>% .$effsize
print(es.value)

mw<- wilcox.test(pred.ts.pos, pred.ts.neg)
mw.p = mw$p.value %>% format.pval(3)
print(mw.p)
es.value <- df3 %>% filter(group2 %in% "Absolute tissue-specific") %>% rstatix::wilcox_effsize(enrich_sum ~ group1, comparisons = list(c("Predicted 3-UTRs", "Predicted non-3-UTRs")), ref.group = "Predicted 3-UTRs") %>% .$effsize
print(es.value)

mw<- wilcox.test(pred.s.pos, known)
mw.p = mw$p.value %>% format.pval(3)
print(mw.p)
es.value <- df3 %>% filter(group2 %in% "Shared") %>% rstatix::wilcox_effsize(enrich_sum ~ group1, comparisons = list(c("Predicted 3-UTRs", "Known 3-UTRs")), ref.group = "Predicted 3-UTRs") %>% .$effsize
print(es.value)


mw<- wilcox.test(pred.s.pos, pred.s.neg)
mw.p = mw$p.value %>% format.pval(3)
print(mw.p)
es.value <- df3 %>% filter(group2 %in% "Shared") %>% rstatix::wilcox_effsize(enrich_sum ~ group1, comparisons = list(c("Predicted 3-UTRs", "Predicted non-3-UTRs")), ref.group = "Predicted 3-UTRs") %>% .$effsize
print(es.value)


png("ridges_groups.png",bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
ggplot(df3, aes(enrich_sum, group2, colour = group1, fill=group1)) +
  geom_density_ridges(alpha=0.5, scale = 1.2, lwd = 0.6)+
  scale_x_continuous(expand = c(0,0), trans="log2", name = "Number of RBP motif binding per 100 nt", breaks = c(0.25, 4, 64, 1024)) +
  scale_color_manual(values= c("black", "#00AFBB", "#E7B800")) +
  scale_fill_manual(values= c("black", "#00AFBB", "#E7B800")) +
  scale_y_discrete(expand = c(0,0), labels = c("Shared", "Absolute\ntissue-specific", "Highly\nbrain-specific"))+
  annotate("text",x=300,y=1.9,size=3.1,label=c("x: p < 2e-16; es = 0.151")) +
  annotate("text",x=300,y=1.7,size=3.1,label=c("y: p < 2.2e-10; es = 0.081")) +
  annotate("text",x=300,y=2.9,size=3.1,label=c("x: p = 3.9e-07; es = 0.039")) +
  annotate("text",x=300,y=2.7,size=3.1,label=c("y: p = 2.5e-05; es = 0.067")) +
  annotate("text",x=300,y=4.1,size=3.1,label=c("x: p < 2e-16; es = 0.18")) +
  annotate("text",x=300,y=3.9,size=3.1,label=c("y: p < 2e-16; es = 0.174")) +
  coord_cartesian(clip = "off") +
  theme_ridges(center_axis_labels = TRUE)  +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.30,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        legend.position= c(-0.4,0.95),
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
