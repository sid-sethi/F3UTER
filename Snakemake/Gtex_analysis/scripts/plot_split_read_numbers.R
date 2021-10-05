args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(RColorBrewer)
  library(ggridges)
})

options(stringsAsFactors=F, warn=-1)

merged_data = args[1]
app_dir = args[2]
formatting_csv <- args[3]
outpath <- args[4]


all_tissues <- readRDS(merged_data)

data <- all_tissues %>%
  dplyr::filter(! tissue %in% c("brain_cerebellum", "cortex")) %>%
  mutate(predicted.response = ifelse(predicted.prob > 0.6, "Predicted 3-UTRs", "Predicted non-3-UTRs"),
         gene_association = ifelse(annotationType_split_read_annot %in% "partially annotated split read", "split-read", "proximity"))


######### distance ER-genes ###########

dist.df <- data %>%
  dplyr::select(predicted.response, distance_from_geneEnd, gene_association) %>%
  mutate(distance_from_geneEnd = abs(distance_from_geneEnd))


dist.all.med <- dist.df %>% .$distance_from_geneEnd %>% median()/1000
dist.prox.med <- dist.df %>% dplyr::filter(gene_association %in% "proximity") %>% .$distance_from_geneEnd %>% median()/1000
dist.split.med <- dist.df %>% dplyr::filter(gene_association %in% "split-read") %>% .$distance_from_geneEnd %>% median()/1000
dist.majority <- dist.df %>% dplyr::filter(distance_from_geneEnd < 3000) %>% nrow()/nrow(dist.df)*100

pos.dist.all.med <- dist.df %>% dplyr::filter(predicted.response %in% "Predicted 3-UTRs") %>% .$distance_from_geneEnd %>% median()/1000
pos.dist.prox.med <- dist.df %>% dplyr::filter(predicted.response %in% "Predicted 3-UTRs", gene_association %in% "proximity") %>% .$distance_from_geneEnd %>% median()/1000
pos.dist.split.med <- dist.df %>% dplyr::filter(predicted.response %in% "Predicted 3-UTRs", gene_association %in% "split-read") %>% .$distance_from_geneEnd %>% median()/1000
pos.dist.majority <- dist.df %>% dplyr::filter(predicted.response %in% "Predicted 3-UTRs", distance_from_geneEnd < 3000) %>% nrow()/nrow(dist.df %>% dplyr::filter(predicted.response %in% "Predicted 3-UTRs"))*100

is_nearest_all_data <- data %>% dplyr::filter(is_associated_equal_nearest %in% "yes") %>% nrow()/nrow(data)*100
proximity_data <- data %>% dplyr::filter(gene_association %in% "proximity")
is_nearest_proximity_data <- proximity_data %>% dplyr::filter(is_associated_equal_nearest %in% "yes") %>% nrow()/nrow(proximity_data)*100

cat("median distance of ER-gene associations - all 3' ERs - ", dist.all.med, " kb\n",
    "median distance of ER-gene associations - all 3' ERs without split-reads - ", dist.prox.med, " kb\n",
    "median distance of ER-gene associations - all 3' ERs with split-reads - ", dist.split.med, " kb\n",
    "% of ER-gene associations (all 3') within 3kb - ", dist.majority, "\n",
    "median distance of ER-gene associations - 3'UTR predicted ERs - ", pos.dist.all.med, " kb\n",
    "median distance of ER-gene associations - 3'UTR predicted ERs without split-reads - ", pos.dist.prox.med, " kb\n",
    "median distance of ER-gene associations - 3'UTR predicted ERs with split-reads - ", pos.dist.split.med, " kb\n",
    "% of ER-gene associations (3'UTR predicted ERs) within 3kb - ", pos.dist.majority, "\n",
    "% of associated genes also nearest gene in all 3'ERs - ", is_nearest_all_data, "\n",
    "% of associated genes also nearest gene in proximity data (ERs without split-reads) - ", is_nearest_proximity_data, "\n",
    file = str_c(outpath, "/ER_distance_and_nearest_gene.txt"))



png(str_c(outpath, "/distance_density_splitread.png"), height = 3.75, width = 7.25, res = 600, units = "in")

ggplot(dist.df, aes(distance_from_geneEnd/1000, predicted.response, fill=gene_association)) +
  stat_density_ridges(alpha=0.7, scale = 1.2, lwd = 0.6, quantile_lines = TRUE, quantiles = 2)+
  scale_fill_brewer(palette = c("RdYlBu"))+
  guides(fill = guide_legend(title = "Gene association", nrow = 2, title.position = "top")) +
  scale_y_discrete(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0), name = "Distance of ER from gene end (kb)", breaks = seq(0, 10, 1))+
  coord_cartesian(clip = "off") +
  theme_ridges(center_axis_labels = TRUE)  +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.35,"cm"),
        legend.title = element_text(size=12),
        legend.text = element_text(size = 12),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size = 12, face = "bold"),
        panel.border = element_rect(colour="BLACK",size=0.4),
        #panel.border = element_blank(),
        axis.title.x = element_text(size=12),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )
dev.off()

########################################################


# overall split read percentage
all.df <- data %>%
 group_by(gene_association) %>%
 summarise(count = n(), .groups = "keep") %>%
 mutate(percentage = (count/nrow(data))*100)


# grouped by predictions
res.df <- data %>%
  group_by(predicted.response, gene_association) %>%
  summarise(count = n(), .groups = "keep")

total_pos <- data %>% dplyr::filter(predicted.response %in% "Predicted 3-UTRs") %>% nrow()

total_neg <- data %>% dplyr::filter(predicted.response %in% "Predicted non-3-UTRs") %>% nrow()

res.df <- res.df %>%
  mutate(percentage = ifelse(predicted.response %in% "Predicted 3-UTRs",
        (count/total_pos)*100,
        (count/total_neg)*100
    )
  )


write.table(all.df, str_c(outpath, "/overall_splitread_percentage.txt"), quote=FALSE, row.names=FALSE, sep="\t")
write.table(res.df, str_c(outpath, "/all_tissues_splitread_percentage.txt"), quote=FALSE, row.names=FALSE, sep="\t")

png(str_c(outpath, "/all_tissues_splitread_percentage.png"), height = 3.75, width = 4.25, res = 600, units = "in")

ggplot(data = res.df, aes(x=predicted.response, y=percentage, fill = factor(gene_association))) +
  geom_bar(width =0.8, stat = "identity", color = "black", lwd = 0.3) +
  theme_bw()  +
  scale_fill_brewer(palette = c("RdYlBu")) +
  guides(fill = guide_legend(title = "Gene association", nrow = 2, title.position = "top")) +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.40,"cm"),
        legend.title = element_text(size=10),
        legend.text = element_text(size = 10),
        legend.direction = "horizontal",
        #legend.position= c(0.90,0.82),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=12, angle=45, vjust = 1, hjust=1, face = "bold"),
        axis.text.y = element_text(size = 10),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_y_continuous(name = "% ERs")

dev.off()




files <- list.files(path = app_dir, pattern = "_appData.rds")
files <- files[!str_detect(files, c("brain_cerebellum", "brain_cortex"))]

geneData <- data.frame()

for(f in files){

  print(str_c(Sys.time(), " - ", f))

  df <- readRDS(str_c(app_dir, "/", f))

  pos.df <- df %>% dplyr::filter(predicted.prob > 0.60) %>%
    mutate(predicted.response = "Predicted 3-UTRs" ,
         gene_association = ifelse(annotationType_split_read_annot %in% "partially annotated split read", "split-read", "proximity"))

  tis <- str_replace(f, "_appData.rds", "") %>% str_replace(., "brain_", "")

  total_genes <- pos.df %>%
    distinct(.$associated_gene) %>% nrow()

  tis.df <- pos.df %>%
    group_by(associated_gene) %>%
    summarise(split_anno_merge = str_c(unique(gene_association), collapse = ":"),
              split_anno_final = ifelse(str_detect(split_anno_merge, "split-read"), "split-read", "proximity"), .groups="keep") %>%
    group_by(split_anno_final) %>%
    summarise(count = n(), .groups = "keep") %>%
    mutate(tissue_name = tis, total_count = total_genes) %>%
    as.data.frame()

  geneData <- rbind(geneData, tis.df)

}



formatting <- read_delim(formatting_csv, delim = ",") %>% mutate(tissue_color_hex = str_c("#", tissue_color_hex)) %>%
  drop_na() %>%
  dplyr::select(gtex_tissues_name_formatted_2, gtex_tissue_group, gtex_tissues_name_to_plot, tissue_color_hex)


# number of genes
gene_count = dplyr::left_join(geneData, formatting, by = c("tissue_name" = "gtex_tissues_name_formatted_2"))

col.df <- data.frame(tissue = levels(reorder(gene_count$gtex_tissues_name_to_plot, -gene_count$total_count))) %>%
  dplyr::left_join(formatting, by = c("tissue" = "gtex_tissues_name_to_plot"))

png(str_c(outpath, "/gene_splitread_pertissue.png"), height = 4.20, width = 6.75, res = 600, units = "in")
ggplot(data = gene_count, aes(x=reorder(gtex_tissues_name_to_plot, -total_count), y=count, fill = split_anno_final)) +
  geom_bar(width =0.8, stat = "identity", color = "black", lwd = 0.3, alpha=0.8) +
  theme_bw()  +
  scale_fill_brewer(palette = c("RdYlBu")) +
  guides(fill = guide_legend(title = "Gene association", nrow = 2, title.position = "top")) +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.35,"cm"),
        legend.title = element_text(size=8),
        legend.text = element_text(size = 8),
        legend.direction = "horizontal",
        legend.position= c(0.90,0.82),
        legend.background = element_rect(fill = "transparent",colour = NA),
        #axis.text.x = element_text(size=8, angle=90, vjust = 0.3, hjust=1, color = col.df$color, face = ifelse(col.df$color == "black", "plain", "bold")),
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
  scale_y_continuous(name = "Number of genes") +
  geom_point(data = col.df, # object that has the order of the tissues to be plotted
             aes(x = tissue, y = -8), # set
             colour = "black", fill = col.df$tissue_color_hex, shape = 22, size = 2)+
  coord_cartesian(ylim = c(0, 320), clip = "off")
dev.off()



# number of genes with at least one unannotated 3'UTR associated via split-read
x <- data %>%
  dplyr::filter(predicted.response %in% "Predicted 3-UTRs") %>%
  group_by(associated_gene) %>%
  summarise(split_anno_merge = str_c(unique(gene_association), collapse = ":"),
            split_anno_final = ifelse(str_detect(split_anno_merge, "split-read"), "split-read", "proximity"), .groups="keep") %>%
  group_by(split_anno_final) %>%
  summarise(count = n(), .groups = "keep") %>%
  mutate(group = "Predicted 3-UTRs")


y <- data %>%
  group_by(associated_gene) %>%
  summarise(split_anno_merge = str_c(unique(gene_association), collapse = ":"),
            split_anno_final = ifelse(str_detect(split_anno_merge, "split-read"), "split-read", "proximity"), .groups="keep") %>%
  group_by(split_anno_final) %>%
  summarise(count = n(), .groups = "keep") %>%
  mutate(group = "All ERs")

xy <- bind_rows(x, y)

write.table(xy, str_c(outpath, "/genes_with_atleast_one_splitread.txt"), quote=FALSE, row.names=FALSE, sep="\t")
