options(stringsAsFactors=F)
library(tidyverse)
library(scales)
library(caret)
library(ggpubr)


ml_table = readRDS(file = "/ML_tables/ml_table.rds")
data <- ml_table

preProcValues <- preProcess(data, method = c("center", "scale"))
data <- predict(preProcValues, data)
df = data[,2:40]


# UMAP

library(umap)

res.umap = umap(df)

#saveRDS(res.umap, file = "/Umap/res_umap_all.rds")
#res.umap <- readRDS("/Umap/res_umap_all.rds")

layout <- data.frame(res.umap$layout) %>% mutate(class=data$class)
layout$class <- recode(layout$class, UTR = "3-UTR") %>% fct_relevel("pseudoGene")

png("/Umap/res_umap_all.png", height = 3.75, width = 3, res = 600, units = "in", bg = "transparent")
ggplot(layout, aes(x=X1, y=X2, color=class, fill=class)) +
  scale_color_manual(values = c("#FF9999", "#404040", "#D1A700", "#CC00CC", "#009900", "#00AFBB")) +
  geom_density_2d() +
  theme_bw() +
  scale_y_continuous(name = "UMAP2", limits = c(-10,4))+
  scale_x_continuous(name = "UMAP1", limits = c(-9,8)) +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.35,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.position= c(0.2,0.15),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size = 9),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        panel.border = element_blank(),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()



# dot/scatter plot
png("/Umap/res_umap_all_scatter.png", height = 5.75, width = 6.25, res = 600, units = "in", bg = "transparent")
ggplot(layout, aes(x=X1, y=X2, color=class, fill=class)) +
  scale_color_manual(values = c("#FF9999", "#404040", "#D1A700", "#CC00CC", "#009900", "#00AFBB")) +
  geom_point(shape=1, size=0.5, alpha=0.5) +
  theme_bw() +
  scale_y_continuous(name = "UMAP2", limits = c(-10,4))+
  scale_x_continuous(name = "UMAP1", limits = c(-9,8)) +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.65,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position= c(0.12,0.18),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size = 11),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        panel.border = element_blank(),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4))
dev.off()
