library(tidyverse)
library(ggsignif)

out.path = ""

# random forest
rf <- read.table("/Summary/overall_stats.txt", header=TRUE, sep="\t")


rf.df <- pivot_longer(rf, -stats, values_to = "values") %>%
  filter(!name %in% c("mean", "sd"), stats %in% c("Accuracy", "Kappa")) %>%
  mutate(method = as.factor("random forest"))

# logistic regression
lr <- read.table("/Summary/overall_stats.txt", header=TRUE, sep="\t")


lr.df <- pivot_longer(lr, -stats, values_to = "values") %>%
  filter(!name %in% c("mean", "sd"), stats %in% c("Accuracy", "Kappa")) %>%
  mutate(method = as.factor("logistic regression"))


df <- rbind(rf.df, lr.df)


kappa.df <- df %>% filter(stats %in% "Kappa")
x <- kappa.df %>%  filter(method %in% "random forest") %>% .$values
y <- kappa.df %>%  filter(method %in% "logistic regression") %>% .$values

man.pval = wilcox.test(x,y)$p.value
label.pval = paste("p ",format.pval(man.pval,3),sep="")

png(str_c(out.path, "rf_vs_lr_kappa.png"),bg="transparent",units="in",width = 2.25, height= 3.75 ,res=600)
ggplot(data=kappa.df, aes(x=factor(method), y=values)) +
  geom_boxplot(fill="#E0E0E0", lwd = 0.55, outlier.shape=NA, alpha=0.8) +
  geom_jitter(position=position_jitter(0.28), shape=21, size=1, alpha=0.4, fill = "black") +
  theme_bw() +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.40,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        legend.position= c(0.85,0.90),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=11,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_y_continuous(name = "Kappa statistic", limits = c(0.46, 0.59)) +
  scale_x_discrete(labels = c("Random\nforest", "Logistic\nregression")) +
  geom_signif(comparisons = list(c("logistic regression", "random forest")),map_signif_level=TRUE, annotation = label.pval, size=0.4,textsize=3,tip_length=0.02, margin_top = 0.05, color = "black", fontface = 2)
dev.off()




accu.df <- df %>% filter(stats %in% "Accuracy")
x <- accu.df %>%  filter(method %in% "random forest") %>% .$values
y <- accu.df %>%  filter(method %in% "logistic regression") %>% .$values

man.pval = wilcox.test(x,y)$p.value
label.pval = paste("p ",format.pval(man.pval,3),sep="")

png(str_c(out.path, "rf_vs_lr_accu.png"),bg="transparent",units="in",width = 2.25, height= 3.75 ,res=600)
ggplot(data=accu.df, aes(x=factor(method), y=values)) +
  geom_boxplot(fill="#E0E0E0", lwd = 0.55, outlier.shape=NA, alpha=0.8) +
  geom_jitter(position=position_jitter(0.28), shape=21, size=1, alpha=0.4, fill = "black") +
  theme_bw() +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.40,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        legend.position= c(0.85,0.90),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=11,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_y_continuous(name = "Accuracy", limits = c(0.70, 0.79)) +
  scale_x_discrete(labels = c("Random\nforest", "Logistic\nregression")) +
  geom_signif(comparisons = list(c("logistic regression", "random forest")),map_signif_level=TRUE, annotation = label.pval, size=0.4,textsize=3,tip_length=0.02, margin_top = 0.05, color = "black", fontface = 2)
dev.off()
