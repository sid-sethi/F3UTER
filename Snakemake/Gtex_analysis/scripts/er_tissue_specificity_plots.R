args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
})
options(stringsAsFactors = FALSE, warn=-1)

numbers_rds <- args[1]
cor_data <- args[2]
formatting_csv <- args[3]
outpath <- args[4]


formatting <- read_delim(formatting_csv, delim = ",") %>% mutate(tissue_color_hex = str_c("#", tissue_color_hex)) %>%
  drop_na() %>%
  dplyr::select(gtex_tissues_name_formatted_2, gtex_tissue_group, gtex_tissues_name_to_plot, tissue_color_hex, OMIM_gtex_name)


cor.data <- readRDS(cor_data)

cor.data <- cor.data %>% column_to_rownames("query_tis") %>%  as.matrix()
png(str_c(outpath, "/all_corData_heatmap.png"), height = 5.75, width = 6.75, res = 600, units = "in")
pheatmap(cor.data,
         fontsize_col = 7,
         fontsize_row = 7,
         treeheight_col = 10,
         treeheight_row = 10)
dev.off()




numbers <- readRDS(numbers_rds)

### plots numbers
# tis_spec
res <- numbers %>% dplyr::filter(group %in% "tissue-specific") %>%
      dplyr::left_join(formatting, by = c("tissue" = "gtex_tissues_name_formatted_2"))


col.df <- data.frame(tissue = levels(reorder(res$gtex_tissues_name_to_plot, -res$frequency))) %>%
  dplyr::left_join(formatting, by = c("tissue" = "gtex_tissues_name_to_plot"))


png(str_c(outpath, "/number_tis_spec_allProb.png"), height = 4.20, width = 6.75, res = 600, units = "in")
ggplot(data = res, aes(x=reorder(gtex_tissues_name_to_plot, -frequency), y=frequency, fill = factor(range))) +
  geom_bar(width =0.8, stat = "identity", color = "black", lwd = 0.3) +
  theme_bw()  +
  scale_fill_brewer(palette = c("RdYlBu")) +
  guides(fill = guide_legend(title = "prediction probability", nrow = 5, title.position = "top")) +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.30,"cm"),
        legend.title = element_text(size=7),
        legend.text = element_text(size = 7),
        legend.direction = "horizontal",
        legend.position= c(0.90,0.82),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=8, angle=90, vjust = 0.3, hjust=1, face = ifelse(col.df$gtex_tissue_group == "brain", "bold", "plain")),
        axis.text.y = element_text(size = 8),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=9,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_y_continuous(name = "# total predictions") +
  geom_point(data = col.df, # object that has the order of the tissues to be plotted
             aes(x = tissue, y = -30), # set
             colour = "black", fill = col.df$tissue_color_hex, shape = 22, size = 2)+
  coord_cartesian(ylim = c(0, 1200), clip = "off")
dev.off()


# brain_spec
res <- numbers %>% dplyr::filter(group %in% "highly brain-specific") %>%
  dplyr::left_join(formatting, by = c("tissue" = "gtex_tissues_name_formatted_2"))

col.df <- data.frame(tissue = levels(reorder(res$gtex_tissues_name_to_plot, -res$frequency))) %>%
  dplyr::left_join(formatting, by = c("tissue" = "gtex_tissues_name_to_plot"))


png(str_c(outpath, "/number_brain_spec_allProb.png"), height = 4.20, width = 6.75, res = 600, units = "in")
ggplot(data = res, aes(x=reorder(gtex_tissues_name_to_plot, -frequency), y=frequency, fill = factor(range))) +
  geom_bar(width =0.8, stat = "identity", color = "black", lwd = 0.3) +
  theme_bw()  +
  scale_fill_brewer(palette = c("RdYlBu")) +
  guides(fill = guide_legend(title = "prediction probability", nrow = 5, title.position = "top")) +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.30,"cm"),
        legend.title = element_text(size=7),
        legend.text = element_text(size = 7),
        legend.direction = "horizontal",
        legend.position= c(0.90,0.82),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=8, angle=90, vjust = 0.3, hjust=1, face = ifelse(col.df$gtex_tissue_group == "brain", "bold", "plain")),
        axis.text.y = element_text(size = 8),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=9,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_y_continuous(name = "# total predictions") +
  geom_point(data = col.df, # object that has the order of the tissues to be plotted
             aes(x = tissue, y = -50), # set
             colour = "black", fill = col.df$tissue_color_hex, shape = 22, size = 2)+
  coord_cartesian(ylim = c(0, 2000), clip = "off")
dev.off()


# common
res <- numbers %>% dplyr::filter(group %in% "shared") %>%
  dplyr::left_join(formatting, by = c("tissue" = "gtex_tissues_name_formatted_2"))

col.df <- data.frame(tissue = levels(reorder(res$gtex_tissues_name_to_plot, -res$frequency))) %>%
  dplyr::left_join(formatting, by = c("tissue" = "gtex_tissues_name_to_plot"))


png(str_c(outpath, "/number_shared_allProb.png"), height = 4.20, width = 6.75, res = 600, units = "in")
ggplot(data = res, aes(x=reorder(gtex_tissues_name_to_plot, -frequency), y=frequency, fill = factor(range))) +
  geom_bar(width =0.8, stat = "identity", color = "black", lwd = 0.3) +
  theme_bw()  +
  scale_fill_brewer(palette = c("RdYlBu")) +
  guides(fill = guide_legend(title = "prediction probability", nrow = 5, title.position = "top")) +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.30,"cm"),
        legend.title = element_text(size=7),
        legend.text = element_text(size = 7),
        legend.direction = "horizontal",
        legend.position= c(0.90,0.82),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=8, angle=90, vjust = 0.3, hjust=1, face = ifelse(col.df$gtex_tissue_group == "brain", "bold", "plain")),
        axis.text.y = element_text(size = 8),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=9,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_y_continuous(name = "# total predictions") +
  geom_point(data = col.df, # object that has the order of the tissues to be plotted
             aes(x = tissue, y = -8), # set
             colour = "black", fill = col.df$tissue_color_hex, shape = 22, size = 2)+
  coord_cartesian(ylim = c(0, 320), clip = "off")
dev.off()


# ambiguous
res <- numbers %>% dplyr::filter(group %in% "ambiguous") %>%
  dplyr::left_join(formatting, by = c("tissue" = "gtex_tissues_name_formatted_2"))

col.df <- data.frame(tissue = levels(reorder(res$gtex_tissues_name_to_plot, -res$frequency))) %>%
  dplyr::left_join(formatting, by = c("tissue" = "gtex_tissues_name_to_plot"))


png(str_c(outpath, "/number_ambigs_allProb.png"), height = 4.20, width = 6.75, res = 600, units = "in")
ggplot(data = res, aes(x=reorder(gtex_tissues_name_to_plot, -frequency), y=frequency, fill = factor(range))) +
  geom_bar(width =0.8, stat = "identity", color = "black", lwd = 0.3) +
  theme_bw()  +
  scale_fill_brewer(palette = c("RdYlBu")) +
  guides(fill = guide_legend(title = "prediction probability", nrow = 5, title.position = "top")) +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.30,"cm"),
        legend.title = element_text(size=7),
        legend.text = element_text(size = 7),
        legend.direction = "horizontal",
        legend.position= c(0.90,0.82),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=8, angle=90, vjust = 0.3, hjust=1, face = ifelse(col.df$gtex_tissue_group == "brain", "bold", "plain")),
        axis.text.y = element_text(size = 8),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=9,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_y_continuous(name = "# total predictions") +
  geom_point(data = col.df, # object that has the order of the tissues to be plotted
             aes(x = tissue, y = -50), # set
             colour = "black", fill = col.df$tissue_color_hex, shape = 22, size = 2)+
  coord_cartesian(ylim = c(0, 2000), clip = "off")
dev.off()



# plot high probability numbers

range_select = c("(0.6,0.7]", "(0.7,0.8]", "(0.8,0.9]", "(0.9,1]")

# extracting the colours used
blue_colors <- brewer.pal(n=10, name = "RdYlBu")[7:10]


## tis_spec
res <- numbers %>% dplyr::filter(group %in% "tissue-specific", range %in% range_select) %>%
  dplyr::left_join(formatting, by = c("tissue" = "gtex_tissues_name_formatted_2"))

col.df <- data.frame(tissue = levels(reorder(res$gtex_tissues_name_to_plot, -res$frequency))) %>%
  dplyr::left_join(formatting, by = c("tissue" = "gtex_tissues_name_to_plot"))

png(str_c(outpath, "/numbers_tis_spec_highProb.png"), height = 4.20, width = 6.75, res = 600, units = "in")
ggplot(data = res, aes(x=reorder(gtex_tissues_name_to_plot, -frequency), y=frequency, fill = factor(range))) +
  geom_bar(width =0.8, stat = "identity", color = "black", lwd = 0.3) +
  theme_bw()  +
  scale_fill_manual(values = blue_colors) +
  guides(fill = guide_legend(title = "prediction probability", nrow = 5, title.position = "top")) +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.30,"cm"),
        legend.title = element_text(size=7),
        legend.text = element_text(size = 7),
        legend.direction = "horizontal",
        legend.position= c(0.90,0.82),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=8, angle=90, vjust = 0.3, hjust=1, face = ifelse(col.df$gtex_tissue_group == "brain", "bold", "plain")),
        axis.text.y = element_text(size = 8),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=9,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_y_continuous(name = "# of 3'UTR predictions") +
  geom_point(data = col.df, # object that has the order of the tissues to be plotted
             aes(x = tissue, y = -1), # set
             colour = "black", fill = col.df$tissue_color_hex, shape = 22, size = 2)+
  coord_cartesian(ylim = c(0, 40), clip = "off")
dev.off()



# brain_spec
res <- numbers %>% dplyr::filter(group %in% "highly brain-specific", range %in% range_select) %>%
  dplyr::left_join(formatting, by = c("tissue" = "gtex_tissues_name_formatted_2"))

col.df <- data.frame(tissue = levels(reorder(res$gtex_tissues_name_to_plot, -res$frequency))) %>%
  dplyr::left_join(formatting, by = c("tissue" = "gtex_tissues_name_to_plot"))


png(str_c(outpath, "/numbers_brain_spec_highProb.png"), height = 4.20, width = 6.75, res = 600, units = "in")
ggplot(data = res, aes(x=reorder(gtex_tissues_name_to_plot, -frequency), y=frequency, fill = factor(range))) +
  geom_bar(width =0.8, stat = "identity", color = "black", lwd = 0.3) +
  theme_bw()  +
  #scale_fill_brewer(palette = c("RdYlBu")) +
  scale_fill_manual(values = blue_colors) +
  guides(fill = guide_legend(title = "prediction probability", nrow = 5, title.position = "top")) +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.30,"cm"),
        legend.title = element_text(size=7),
        legend.text = element_text(size = 7),
        legend.direction = "horizontal",
        legend.position= c(0.90,0.82),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=8, angle=90, vjust = 0.3, hjust=1, face = ifelse(col.df$gtex_tissue_group == "brain", "bold", "plain")),
        axis.text.y = element_text(size = 8),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=9,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_y_continuous(name = "# of 3'UTR predictions") +
  geom_point(data = col.df, # object that has the order of the tissues to be plotted
             aes(x = tissue, y = -4), # set
             colour = "black", fill = col.df$tissue_color_hex, shape = 22, size = 2)+
  coord_cartesian(ylim = c(0, 175), clip = "off")
dev.off()


# common
res <- numbers %>% dplyr::filter(group %in% "shared", range %in% range_select) %>%
  dplyr::left_join(formatting, by = c("tissue" = "gtex_tissues_name_formatted_2"))

col.df <- data.frame(tissue = levels(reorder(res$gtex_tissues_name_to_plot, -res$frequency))) %>%
  dplyr::left_join(formatting, by = c("tissue" = "gtex_tissues_name_to_plot"))


png(str_c(outpath, "/numbers_shared_highProb.png"), height = 4.20, width = 6.75, res = 600, units = "in")
ggplot(data = res, aes(x=reorder(gtex_tissues_name_to_plot, -frequency), y=frequency, fill = factor(range))) +
  geom_bar(width =0.8, stat = "identity", color = "black", lwd = 0.3) +
  theme_bw()  +
  #scale_fill_brewer(palette = c("RdYlBu")) +
  scale_fill_manual(values = blue_colors) +
  guides(fill = guide_legend(title = "prediction probability", nrow = 5, title.position = "top")) +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.30,"cm"),
        legend.title = element_text(size=7),
        legend.text = element_text(size = 7),
        legend.direction = "horizontal",
        legend.position= c(0.90,0.82),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=8, angle=90, vjust = 0.3, hjust=1, face = ifelse(col.df$gtex_tissue_group == "brain", "bold", "plain")),
        axis.text.y = element_text(size = 8),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=9,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_y_continuous(name = "# of 3'UTR predictions") +
  geom_point(data = col.df, # object that has the order of the tissues to be plotted
             aes(x = tissue, y = -0.8), # set
             colour = "black", fill = col.df$tissue_color_hex, shape = 22, size = 2)+
  coord_cartesian(ylim = c(0, 35), clip = "off")
dev.off()


# ambiguous
res <- numbers %>% dplyr::filter(group %in% "ambiguous", range %in% range_select) %>%
  dplyr::left_join(formatting, by = c("tissue" = "gtex_tissues_name_formatted_2"))

col.df <- data.frame(tissue = levels(reorder(res$gtex_tissues_name_to_plot, -res$frequency))) %>%
  dplyr::left_join(formatting, by = c("tissue" = "gtex_tissues_name_to_plot"))


png(str_c(outpath, "/numbers_ambigs_highProb.png"), height = 4.20, width = 6.75, res = 600, units = "in")
ggplot(data = res, aes(x=reorder(gtex_tissues_name_to_plot, -frequency), y=frequency, fill = factor(range))) +
  geom_bar(width =0.8, stat = "identity", color = "black", lwd = 0.3) +
  theme_bw()  +
  #scale_fill_brewer(palette = c("RdYlBu")) +
  scale_fill_manual(values = blue_colors) +
  guides(fill = guide_legend(title = "prediction probability", nrow = 5, title.position = "top")) +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.30,"cm"),
        legend.title = element_text(size=7),
        legend.text = element_text(size = 7),
        legend.direction = "horizontal",
        legend.position= c(0.90,0.82),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=8, angle=90, vjust = 0.3, hjust=1, face = ifelse(col.df$gtex_tissue_group == "brain", "bold", "plain")),
        axis.text.y = element_text(size = 8),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=9,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_y_continuous(name = "# of 3'UTR predictions") +
  geom_point(data = col.df, # object that has the order of the tissues to be plotted
             aes(x = tissue, y = -4), # set
             colour = "black", fill = col.df$tissue_color_hex, shape = 22, size = 2)+
  coord_cartesian(ylim = c(0, 150), clip = "off")
dev.off()




# plot percentages , only prob>0.60

## tis_spec
res <- numbers %>% dplyr::filter(group %in% "tissue-specific", range %in% range_select) %>%
  dplyr::left_join(formatting, by = c("tissue" = "gtex_tissues_name_formatted_2"))

col.df <- data.frame(tissue = levels(reorder(res$gtex_tissues_name_to_plot, -res$perc_tissue))) %>%
  dplyr::left_join(formatting, by = c("tissue" = "gtex_tissues_name_to_plot"))


png(str_c(outpath, "/percentage_tis_spec_highProb.png"), height = 4.20, width = 6.75, res = 600, units = "in")
ggplot(data = res, aes(x=reorder(gtex_tissues_name_to_plot, -perc_tissue), y=perc_tissue, fill = factor(range))) +
  geom_bar(width =0.8, stat = "identity", color = "black", lwd = 0.3) +
  theme_bw()  +
  scale_fill_manual(values = blue_colors) +
  guides(fill = guide_legend(title = "prediction probability", nrow = 5, title.position = "top")) +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.30,"cm"),
        legend.title = element_text(size=7),
        legend.text = element_text(size = 7),
        legend.direction = "horizontal",
        legend.position= c(0.90,0.82),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=8, angle=90, vjust = 0.3, hjust=1, face = ifelse(col.df$gtex_tissue_group == "brain", "bold", "plain")),
        axis.text.y = element_text(size = 8),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=9,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_y_continuous(name = "% of total ER predictions") +
  geom_point(data = col.df, # object that has the order of the tissues to be plotted
             aes(x = tissue, y = -0.05), # set
             colour = "black", fill = col.df$tissue_color_hex, shape = 22, size = 2)+
  coord_cartesian(ylim = c(0, 2), clip = "off")
dev.off()



# brain_spec
res <- numbers %>% dplyr::filter(group %in% "highly brain-specific", range %in% range_select) %>%
  dplyr::left_join(formatting, by = c("tissue" = "gtex_tissues_name_formatted_2"))

col.df <- data.frame(tissue = levels(reorder(res$gtex_tissues_name_to_plot, -res$perc_tissue))) %>%
  dplyr::left_join(formatting, by = c("tissue" = "gtex_tissues_name_to_plot"))


png(str_c(outpath, "/percentage_brain_spec_highProb.png"), height = 4.20, width = 6.75, res = 600, units = "in")
ggplot(data = res, aes(x=reorder(gtex_tissues_name_to_plot, -perc_tissue), y=perc_tissue, fill = factor(range))) +
  geom_bar(width =0.8, stat = "identity", color = "black", lwd = 0.3) +
  theme_bw()  +
  scale_fill_manual(values = blue_colors) +
  guides(fill = guide_legend(title = "prediction probability", nrow = 5, title.position = "top")) +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.30,"cm"),
        legend.title = element_text(size=7),
        legend.text = element_text(size = 7),
        legend.direction = "horizontal",
        legend.position= c(0.90,0.82),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=8, angle=90, vjust = 0.3, hjust=1, face = ifelse(col.df$gtex_tissue_group == "brain", "bold", "plain")),
        axis.text.y = element_text(size = 8),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=9,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_y_continuous(name = "% of total ER predictions") +
  geom_point(data = col.df, # object that has the order of the tissues to be plotted
             aes(x = tissue, y = -0.12), # set
             colour = "black", fill = col.df$tissue_color_hex, shape = 22, size = 2)+
  coord_cartesian(ylim = c(0, 5), clip = "off")
dev.off()


# common
res <- numbers %>% dplyr::filter(group %in% "shared", range %in% range_select) %>%
  dplyr::left_join(formatting, by = c("tissue" = "gtex_tissues_name_formatted_2"))

col.df <- data.frame(tissue = levels(reorder(res$gtex_tissues_name_to_plot, -res$perc_tissue))) %>%
  dplyr::left_join(formatting, by = c("tissue" = "gtex_tissues_name_to_plot"))


png(str_c(outpath, "/percentage_shared_highProb.png"), height = 4.20, width = 6.75, res = 600, units = "in")
ggplot(data = res, aes(x=reorder(gtex_tissues_name_to_plot, -perc_tissue), y=perc_tissue, fill = factor(range))) +
  geom_bar(width =0.8, stat = "identity", color = "black", lwd = 0.3) +
  theme_bw()  +
  scale_fill_manual(values = blue_colors) +
  guides(fill = guide_legend(title = "prediction probability", nrow = 5, title.position = "top")) +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.30,"cm"),
        legend.title = element_text(size=7),
        legend.text = element_text(size = 7),
        legend.direction = "horizontal",
        legend.position= c(0.90,0.82),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=8, angle=90, vjust = 0.3, hjust=1, face = ifelse(col.df$gtex_tissue_group == "brain", "bold", "plain")),
        axis.text.y = element_text(size = 8),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=9,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_y_continuous(name = "% of total ER predictions") +
  geom_point(data = col.df, # object that has the order of the tissues to be plotted
             aes(x = tissue, y = -0.17), # set
             colour = "black", fill = col.df$tissue_color_hex, shape = 22, size = 2)+
  coord_cartesian(ylim = c(0, 6), clip = "off")
dev.off()


# ambiguous
res <- numbers %>% dplyr::filter(group %in% "ambiguous", range %in% range_select) %>%
  dplyr::left_join(formatting, by = c("tissue" = "gtex_tissues_name_formatted_2"))

col.df <- data.frame(tissue = levels(reorder(res$gtex_tissues_name_to_plot, -res$perc_tissue))) %>%
  dplyr::left_join(formatting, by = c("tissue" = "gtex_tissues_name_to_plot"))


png(str_c(outpath, "/percenatge_ambigs_highProb.png"), height = 4.20, width = 6.75, res = 600, units = "in")
ggplot(data = res, aes(x=reorder(gtex_tissues_name_to_plot, -perc_tissue), y=perc_tissue, fill = factor(range))) +
  geom_bar(width =0.8, stat = "identity", color = "black", lwd = 0.3) +
  theme_bw()  +
  scale_fill_manual(values = blue_colors) +
  guides(fill = guide_legend(title = "prediction probability", nrow = 5, title.position = "top")) +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.30,"cm"),
        legend.title = element_text(size=7),
        legend.text = element_text(size = 7),
        legend.direction = "horizontal",
        legend.position= c(0.90,0.82),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=8, angle=90, vjust = 0.3, hjust=1, face = ifelse(col.df$gtex_tissue_group == "brain", "bold", "plain")),
        axis.text.y = element_text(size = 8),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=9,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_y_continuous(name = "% of total ER predictions") +
  geom_point(data = col.df, # object that has the order of the tissues to be plotted
             aes(x = tissue, y = -0.13), # set
             colour = "black", fill = col.df$tissue_color_hex, shape = 22, size = 2)+
  coord_cartesian(ylim = c(0, 5), clip = "off")
dev.off()
