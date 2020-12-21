library(tidyverse)

ml_table = readRDS(file = "/ML_tables/ml_table.rds")

ml_table$class <- recode(ml_table$class, UTR = "3-UTR") %>% fct_relevel("3-UTR")


##########################
### for KW test ##########
pval_df <- data.frame()
i = 0
##########################
##########################



## poly(A) signal

utr.hits = ml_table %>% dplyr::filter(class %in% "3-UTR", polyA_signal == 1)
ice.hits = ml_table %>% dplyr::filter(class %in% "ICE", polyA_signal == 1)
lncrna.hits = ml_table %>% dplyr::filter(class %in% "lncRNA", polyA_signal == 1)
ncrna.hits = ml_table %>% dplyr::filter(class %in% "ncRNA", polyA_signal == 1)
pseudo.hits = ml_table %>% dplyr::filter(class %in% "pseudoGene", polyA_signal == 1)
five.hits = ml_table %>% dplyr::filter(class %in% "5-UTR", polyA_signal == 1)


df = data.frame(number = c(nrow(utr.hits), nrow(five.hits), nrow(ice.hits), nrow(lncrna.hits), nrow(ncrna.hits), nrow(pseudo.hits)),
                group = c("3-UTR", "5-UTR", "ICE", "lncRNA", "ncRNA", "pseudoGene")
)

df$total <- table(ml_table$class) %>%  as.data.frame() %>% .$Freq
df$percentage <- (df$number/df$total)*100


p <- prop.test(x = df$number,
          n = df$total,
          alternative = "two.sided",
          correct = FALSE)


png("/Univariate_feature_plots/polya.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=df, aes(x=group, y=percentage, fill = group, color = group)) +
  geom_bar(width =0.7,stat="identity", alpha=0.8) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  geom_text(label=df$number,vjust=-0.3, size= 1.5, color = "black") +
  theme_bw() +
  guides(fill=FALSE, color = FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()


png("/Univariate_feature_plots/legend.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=df, aes(x=group, y=percentage, fill = group, color = group)) +
  geom_bar(width =0.7,stat="identity", alpha=0.8) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  geom_text(label=df$number,vjust=-0.3, size= 1.5, color = "black") +
  theme_bw() +
  guides(color = FALSE)+
  theme(
    legend.key.size = unit(0.35,"cm"),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    #legend.position= c(0.83,0.91),
    legend.background = element_rect(fill = "transparent",colour = NA),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()



k <- kruskal.test(mean_phastCons7way ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "mean_phastCons7way"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/phastCons.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=mean_phastCons7way)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()



k <- kruskal.test(aphilicity ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "aphilicity"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/aphilicity.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=aphilicity)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()


k <- kruskal.test(basestacking ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "basestacking"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/basestacking.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=basestacking)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()


k <- kruskal.test(bdnatwist ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "bdnatwist"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/bdnatwist.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=bdnatwist)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()


k <- kruskal.test(bendability ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "bendability"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/bendability.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=bendability)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()


k <- kruskal.test(bendingstiffness ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "bendingstiffness"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/bendingstiffness.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=bendingstiffness)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()


k <- kruskal.test(cpgislands ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "cpgislands"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/cpgislands.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=cpgislands)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()


k <- kruskal.test(cpnpgcpgislands ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "cpnpgcpgislands"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/cpnpgcpgislands.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=cpnpgcpgislands)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()


k <- kruskal.test(cpnpgislands ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "cpnpgislands"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/cpnpgislands.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=cpnpgislands)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()


k <- kruskal.test(dnadenaturation ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "dnadenaturation"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/dnadenaturation.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=dnadenaturation)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()


k <- kruskal.test(duplexstabilitydisruptenergy ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "duplexstabilitydisruptenergy"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/duplexstabilitydisruptenergy.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=duplexstabilitydisruptenergy)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()


k <- kruskal.test(duplexstabilityfreeenergy ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "duplexstabilityfreeenergy"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/duplexstabilityfreeenergy.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=duplexstabilityfreeenergy)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()


k <- kruskal.test(nucleosomeposition ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "nucleosomeposition"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/nucleosomeposition.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=nucleosomeposition)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()


k <- kruskal.test(propellortwist ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "propellortwist"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/propellortwist.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=propellortwist)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()


k <- kruskal.test(proteindeformation ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "proteindeformation"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/proteindeformation.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=proteindeformation)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()



k <- kruskal.test(proteindnatwist ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "proteindnatwist"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/proteindnatwist.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=proteindnatwist)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()



k <- kruskal.test(zdna ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "zdna"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/zdna.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=zdna)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()



k <- kruskal.test(RepeatsOverlap ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "RepeatsOverlap"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/RepeatsOverlap.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=RepeatsOverlap)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()



k <- kruskal.test(A ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "A"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/A.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=A)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()



k <- kruskal.test(T ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "T"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/T.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=T)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()



k <- kruskal.test(G ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "G"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/G.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=G)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()



k <- kruskal.test(C ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "C"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/C.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=C)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()



k <- kruskal.test(AA ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "AA"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/AA.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=AA)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()



k <- kruskal.test(TT ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "TT"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/TT.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=TT)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()



k <- kruskal.test(GG ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "GG"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/GG.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=GG)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()



k <- kruskal.test(CC ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "CC"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/CC.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=CC)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()



k <- kruskal.test(AT ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "AT"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/AT.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=AT)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()



k <- kruskal.test(AG ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "AG"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/AG.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=AG)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()



k <- kruskal.test(TG ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "TG"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/TG.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=TG)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()



k <- kruskal.test(GT ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "GT"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/GT.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=GT)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()



k <- kruskal.test(CT ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "CT"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/CT.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=CT)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()



k <- kruskal.test(AC ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "AC"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/AC.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=AC)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()



k <- kruskal.test(TC ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "TC"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/TC.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=TC)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()



k <- kruskal.test(GC ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "GC"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/GC.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=GC)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()



k <- kruskal.test(CG ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "CG"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/CG.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=CG)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()



k <- kruskal.test(TA ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "TA"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/TA.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=TA)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()



k <- kruskal.test(GA ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "GA"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/GA.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=GA)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()



k <- kruskal.test(CA ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "CA"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/CA.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=CA)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()



k <- kruskal.test(meanPd ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "meanPd"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/meanPd.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=meanPd)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()



k <- kruskal.test(meanEE ~ class, data = ml_table) %>% .$p.value %>% format.pval()
i = i + 1
pval_df[i, "feature"] <- "meanEE"
pval_df[i, "pval"] <- k

png("/Univariate_feature_plots/meanEE.png", height = 1.25, width = 1.75, res = 600, units = "in", bg = "transparent")
ggplot(data=ml_table, aes(x=class, y=meanEE)) +
  geom_boxplot(aes(fill=class), lwd = 0.55, outlier.shape=16, alpha=0.8, outlier.size = 1.1, outlier.alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999" )) +
  scale_fill_manual(values = c("#00AFBB", "#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999")) +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))+
  coord_cartesian(ylim=c(0.90,1))
dev.off()
