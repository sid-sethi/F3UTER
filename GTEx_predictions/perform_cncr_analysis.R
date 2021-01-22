args <- commandArgs(TRUE) # 1ST = data path /CNCR_scores/
library(tidyverse)
library(rstatix)
library(ggpubr)
library(ggsignif)
library(ggridges)

print(args[1])
setwd(args[1])

ambi.pos <- read.table("ambi_pos_data.txt", header=TRUE) %>%
  mutate(group1 = "3-UTR", group2 = "Ambiguous")
ambi.neg <- read.table("ambi_neg_data.txt", header=TRUE) %>%
  mutate(group1 = "Non-3-UTR", group2 = "Ambiguous")
ts.pos <- read.table("tis_spec_pos_data.txt", header=TRUE) %>%
  mutate(group1 = "3-UTR", group2 = "Absolute tissue-specific")
ts.neg <- read.table("tis_spec_neg_data.txt", header=TRUE) %>%
  mutate(group1 = "Non-3-UTR", group2 = "Absolute tissue-specific")
bs.pos <- read.table("brain_spec_pos_data.txt", header=TRUE) %>%
  mutate(group1 = "3-UTR", group2 = "Highly brain-specific")
bs.neg <- read.table("brain_spec_neg_data.txt", header=TRUE) %>%
  mutate(group1 = "Non-3-UTR", group2 = "Highly brain-specific")
com.pos <- read.table("common_pos_data.txt", header=TRUE) %>%
  mutate(group1 = "3-UTR", group2 = "Shared")
com.neg <- read.table("common_neg_data.txt", header=TRUE) %>%
  mutate(group1 = "Non-3-UTR", group2 = "Shared")


df <- bind_rows(ambi.pos, ambi.neg, ts.pos, ts.neg, bs.pos, bs.neg, com.pos, com.neg)
df$group1 <- ifelse(df$group1 %in% "3-UTR", "Predicted 3-UTRs", df$group1)
df$group1 <- ifelse(df$group1 %in% "Non-3-UTR", "Predicted Non-3-UTRs", df$group1)

# loading CNCR scores of known 3'UTRs
utr = read.table("/cncr_scores_all_three_UTRs.txt", header=TRUE, sep="\t") %>%
  mutate(group1 = "Known 3-UTRs") %>%
  plyr::rename(c("three_prime_utr_id" = "id"))

df2 = df %>% dplyr::select(id, cnc, group1) %>%
  rbind(utr)


# removing negative predictions from this analysis
df3 <- df2 %>% dplyr::filter(!group1 %in% "Predicted Non-3-UTRs")

pred <- df3 %>% filter(group1 %in% "Predicted 3-UTRs") %>% .$cnc
known <- df3 %>% filter(group1 %in% "Known 3-UTRs") %>% .$cnc


# man-whitney test
mw<- wilcox.test(pred, known)
mw.p = mw$p.value
label.pval = paste("p = ",format.pval(mw.p,3),sep="")
es <- df3 %>% rstatix::wilcox_effsize(cnc ~ group1) %>% .$effsize
label.es = paste("es = ", round(es, 3),sep="")
label = str_c(label.pval, label.es, sep="\n")


png("density_all.png",bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
ggdensity(df3, x = "cnc", y = "..density..",
          #add = "mean",
          rug = TRUE,
          color = "group1", fill = "group1",
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




# split by categories
utr.dummy <- rbind(
  utr %>% mutate(group2 = "Absolute tissue-specific"),
  utr %>% mutate(group2 = "Highly brain-specific"),
  utr %>% mutate(group2 = "Shared")
)
df.dummy <- df %>% dplyr::filter(!group1 %in% "Predicted Non-3-UTRs", !group2 %in% "Ambiguous")

df4 <- rbind(df.dummy, utr.dummy)


#### ggridges
df4$group2 <- as.factor(df4$group2)
df4$group2 <- fct_relevel(df4$group2, "Shared", "Absolute tissue-specific", "Highly brain-specific")

# p-values of comparisons
pred.bs <- df4 %>% filter(group2 %in% "Highly brain-specific", group1 %in% "Predicted 3-UTRs") %>% .$cnc
pred.ts <- df4 %>% filter(group2 %in% "Absolute tissue-specific", group1 %in% "Predicted 3-UTRs") %>% .$cnc
pred.s <- df4 %>% filter(group2 %in% "Shared", group1 %in% "Predicted 3-UTRs") %>% .$cnc


# wilcoxon and effect size
mw<- wilcox.test(pred.bs, known)
mw.p = mw$p.value %>% format.pval(3)
print(mw.p)
es.value <- df4 %>% filter(group2 %in% "Highly brain-specific") %>% rstatix::wilcox_effsize(cnc ~ group1) %>% .$effsize
print(es.value)

mw<- wilcox.test(pred.ts, known)
mw.p = mw$p.value %>% format.pval(3)
print(mw.p)
es.value <- df4 %>% filter(group2 %in% "Absolute tissue-specific") %>% rstatix::wilcox_effsize(cnc ~ group1) %>% .$effsize
print(es.value)

mw<- wilcox.test(pred.s, known)
mw.p = mw$p.value %>% format.pval(3)
print(mw.p)
es.value <- df4 %>% filter(group2 %in% "Shared") %>% rstatix::wilcox_effsize(cnc ~ group1) %>% .$effsize
print(es.value)


png("ridges_groups.png",bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
ggplot(df4, aes(cnc, group2, colour = group1, fill=group1)) +
  geom_density_ridges(alpha=0.5, scale = 1.2, lwd = 0.6)+
  scale_color_manual(values= c("black", "#00AFBB")) +
  scale_fill_manual(values= c("black", "#00AFBB")) +
  scale_y_discrete(expand = c(0,0), labels = c("Shared", "Absolute\ntissue-specific", "Highly\nbrain-specific"))+
  annotate("text",x=-7.5,y=1.2,size=3.5,label=c("p = 0.182\nes = 0.01")) +
  annotate("text",x=-7.5,y=2.2,size=3.5,label=c("p = 0.632\nes = 0.004")) +
  annotate("text",x=-7.5,y=3.2,size=3.5,label=c("p = 0.121\nes = 0.011")) +
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
