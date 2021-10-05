args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(RColorBrewer)
})

options(stringsAsFactors=F)

app_dir = args[1]
formatting_csv <- args[2]
outpath <- args[3]

files <- list.files(path = app_dir, pattern = "_appData.rds")
files <- files[!str_detect(files, c("brain_cerebellum", "brain_cortex"))]


res60 = data.frame()
res = data.frame()
distData = data.frame()
all_genes = data.frame()

i = 0
for(f in files){
  i= i + 1
  data <- readRDS(str_c(app_dir, "/", f))

  pos = data %>% dplyr::filter(predicted.prob > 0.60) %>% nrow()
  neg = data %>% dplyr::filter(predicted.prob <= 0.60) %>% nrow()
  tis <- str_replace(f, "_appData.rds", "") %>% str_replace(., "brain_", "")

  pos.length = data %>% dplyr::filter(predicted.prob > 0.60) %>% .$width %>% sum()

  res60[i, "tissue"] <- tis
  res60[i, "UTR"] <- pos
  res60[i, "Non-3-UTR"] <- neg
  res60[i, "total"] <- nrow(data)
  res60[i, "percent_utr"] <- round((pos/nrow(data))*100,2)
  res60[i, "utr_total_length"] <- pos.length
  res60[i, "utr_gene_freq"] <- data %>% dplyr::filter(predicted.prob > 0.60) %>% distinct(.$associated_gene) %>% nrow()


  geneNames <- data %>% dplyr::filter(predicted.prob > 0.60) %>% distinct(associated_gene) %>%
    mutate(tissue = tis)
  all_genes <- rbind(all_genes, geneNames)

  # all prob
  breaks = seq(0,1,0.10)

  df = data %>%
    group_by(range=cut(predicted.prob, breaks= seq(0, 1, by = 0.10)) ) %>%
    summarise(frequency= n(), length = sum(width), Gene_freq= n_distinct(associated_gene)) %>%
    arrange(as.numeric(range)) %>%
    mutate(threshold = breaks[-1], tissue = tis) %>%
    as.data.frame()


  res = rbind(res, df)

  dist = data %>% dplyr::select(predicted.prob, associated_gene, distance_from_geneEnd, tissue)
  distData = rbind(distData, dist)

}

string_out <- data.frame()

string_out[1, "metric"] <- "average % of pos predictions across all tissues"
string_out[1, "value"] <- res60$percent_utr %>% mean()

string_out[2, "metric"] <- "average number of pos predictions across all tissues"
string_out[2, "value"] <- res60$UTR %>% mean()

r <- range(res60$UTR)
string_out[3, "metric"] <- "range of pos predictions across all tissues"
string_out[3, "value"] <- str_c(r[1], ":", r[2])

string_out[4, "metric"] <- "average length of pos predictions in kb"
string_out[4, "value"] <- res60$utr_total_length %>% mean()/1000

string_out[5, "metric"] <- "total length of pos predictions in Mb"
string_out[5, "value"] <- res60$utr_total_length %>% sum()/1000000

r <- range(res60$utr_total_length)/1000
string_out[6, "metric"] <- "range of length of pos predictions in kb"
string_out[6, "value"] <- str_c(r[1], ":", r[2])

string_out[7, "metric"] <- "average number of genes"
string_out[7, "value"] <- res60$utr_gene_freq %>% mean()

r <- range(res60$utr_gene_freq)
string_out[8, "metric"] <- "range of number of genes"
string_out[8, "value"] <- str_c(r[1], ":", r[2])


brain_tissues <- c("amygdala", "anterior_cingulate_cortex_ba24", "caudate_basal_ganglia", "cerebellar_hemisphere", "cerebellum", "cortex", "frontal_cortex_ba9", "hippocampus", "hypothalamus", "nucleus_accumbens_basal_ganglia", "putamen_basal_ganglia", "spinal_cord_cervical_c_1", "substantia_nigra")

mean_brain_tissues <- res60 %>% dplyr::filter(tissue %in% brain_tissues) %>% .$percent_utr %>% mean()

string_out[9, "metric"] <- "average % of predictions across brain tissues"
string_out[9, "value"] <- mean_brain_tissues

res60.brain <- res60 %>% dplyr::filter(tissue %in% brain_tissues)
res60.nbrain <- res60 %>% dplyr::filter(!tissue %in% brain_tissues)
man.pval.number = wilcox.test(res60.brain$UTR,res60.nbrain$UTR)$p.value %>% format.pval(3)
man.pval.length = wilcox.test(res60.brain$utr_total_length,res60.nbrain$utr_total_length)$p.value %>% format.pval(3)
man.pval.gene_freq = wilcox.test(res60.brain$utr_gene_freq,res60.nbrain$utr_gene_freq)$p.value %>% format.pval(3)

string_out[10, "metric"] <- "median number of brain utrs"
string_out[10, "value"] <- median(res60.brain$UTR)

string_out[11, "metric"] <- "median number of non brain utrs"
string_out[11, "value"] <- median(res60.nbrain$UTR)

string_out[12, "metric"] <- "median length of brain utrs in kb"
string_out[12, "value"] <- median(res60.brain$utr_total_length)/1000

string_out[13, "metric"] <- "median length of non brain utrs in kb"
string_out[13, "value"] <- median(res60.nbrain$utr_total_length)/1000

string_out[14, "metric"] <- "median number of brain utrs genes"
string_out[14, "value"] <- median(res60.brain$utr_gene_freq)

string_out[15, "metric"] <- "median number of non brain utrs genes"
string_out[15, "value"] <- median(res60.nbrain$utr_gene_freq)

string_out[16, "metric"] <- "total distinct genes with unannotated 3'UTR"
string_out[16, "value"] <- all_genes$associated_gene %>% unique() %>% length()

string_out[17, "metric"] <- "pval number of utrs in brain vs non brain"
string_out[17, "value"] <- man.pval.number

string_out[18, "metric"] <- "pval length of utrs in brain vs non brain"
string_out[18, "value"] <- man.pval.length

string_out[19, "metric"] <- "pval gene-freq of utrs in brain vs non brain"
string_out[19, "value"] <- man.pval.gene_freq

write.table(string_out, file = str_c(outpath, "/numbers.txt"), quote=FALSE, row.names=FALSE, sep="\t")



print(str_c(Sys.time(), " - plotting ....."))

formatting <- read_delim(formatting_csv, delim = ",") %>% mutate(tissue_color_hex = str_c("#", tissue_color_hex)) %>%
  drop_na() %>%
  dplyr::select(gtex_tissues_name_formatted_2, gtex_tissue_group, gtex_tissues_name_to_plot, tissue_color_hex)


res = dplyr::left_join(res, formatting, by = c("tissue" = "gtex_tissues_name_formatted_2"))


col.df <- data.frame(tissue = levels(reorder(res$gtex_tissues_name_to_plot, -res$frequency))) %>%
  dplyr::left_join(formatting, by = c("tissue" = "gtex_tissues_name_to_plot"))


png(str_c(outpath, "/number_predictions_allProb.png"), height = 4.20, width = 6.75, res = 600, units = "in")
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
             aes(x = tissue, y = -100), # set
             colour = "black", fill = col.df$tissue_color_hex, shape = 22, size = 2)+
  coord_cartesian(ylim = c(0, 4000), clip = "off")
dev.off()


# all prob - length of predictions

col.df <- data.frame(tissue = levels(reorder(res$gtex_tissues_name_to_plot, -res$length))) %>%
  dplyr::left_join(formatting, by = c("tissue" = "gtex_tissues_name_to_plot"))

png(str_c(outpath, "/length_predictions_allProb.png"), height = 4.20, width = 6.75, res = 600, units = "in")
ggplot(data = res, aes(x=reorder(gtex_tissues_name_to_plot, -length), y=length/1000, fill = factor(range))) +
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
  scale_y_continuous(name = "Total length of predictions (kb)") +
  geom_point(data = col.df, # object that has the order of the tissues to be plotted
             aes(x = tissue, y = -35), # set
             colour = "black", fill = col.df$tissue_color_hex, shape = 22, size = 2)+
  coord_cartesian(ylim = c(0, 1500), clip = "off")
dev.off()



# plot high probability numbers

range_select = c("(0.6,0.7]", "(0.7,0.8]", "(0.8,0.9]", "(0.9,1]")

# extracting the colours used
blue_colors <- brewer.pal(n=10, name = "RdYlBu")[7:10]


##
res_high <- res %>% dplyr::filter(range %in% range_select)


col.df <- data.frame(tissue = levels(reorder(res_high$gtex_tissues_name_to_plot, -res_high$frequency))) %>%
  dplyr::left_join(formatting, by = c("tissue" = "gtex_tissues_name_to_plot"))


png(str_c(outpath, "/number_predictions_highProb.png"), height = 4.20, width = 6.75, res = 600, units = "in")
ggplot(data = res_high, aes(x=reorder(gtex_tissues_name_to_plot, -frequency), y=frequency, fill = factor(range))) +
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
             aes(x = tissue, y = -8), # set
             colour = "black", fill = col.df$tissue_color_hex, shape = 22, size = 2)+
  coord_cartesian(ylim = c(0, 350), clip = "off")
dev.off()


# high prob - length of predictions

col.df <- data.frame(tissue = levels(reorder(res_high$gtex_tissues_name_to_plot, -res_high$length))) %>%
  dplyr::left_join(formatting, by = c("tissue" = "gtex_tissues_name_to_plot"))

png(str_c(outpath, "/length_predictions_highProb.png"), height = 4.20, width = 6.75, res = 600, units = "in")
ggplot(data = res_high, aes(x=reorder(gtex_tissues_name_to_plot, -length), y=length/1000, fill = factor(range))) +
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
  scale_y_continuous(name = "Total length of 3'UTR predictions (kb)") +
  geom_point(data = col.df, # object that has the order of the tissues to be plotted
             aes(x = tissue, y = -7), # set
             colour = "black", fill = col.df$tissue_color_hex, shape = 22, size = 2)+
  coord_cartesian(ylim = c(0, 270), clip = "off")
dev.off()


# high prob - number of genes
# this plot has discrepency, as unique genes are calculated for each range, and not overall, so numbers are a bit inflated

col.df <- data.frame(tissue = levels(reorder(res_high$gtex_tissues_name_to_plot, -res_high$Gene_freq))) %>%
  dplyr::left_join(formatting, by = c("tissue" = "gtex_tissues_name_to_plot"))

png(str_c(outpath, "/gene_predictions_highProb.png"), height = 4.20, width = 6.75, res = 600, units = "in")
ggplot(data = res_high, aes(x=reorder(gtex_tissues_name_to_plot, -Gene_freq), y=Gene_freq, fill = factor(range))) +
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


# number of genes
res60.df = dplyr::left_join(res60, formatting, by = c("tissue" = "gtex_tissues_name_formatted_2"))

col.df <- data.frame(tissue = levels(reorder(res60.df$gtex_tissues_name_to_plot, -res60.df$utr_gene_freq))) %>%
  dplyr::left_join(formatting, by = c("tissue" = "gtex_tissues_name_to_plot"))

png(str_c(outpath, "/gene_predictions_numbers.png"), height = 4.20, width = 6.75, res = 600, units = "in")
ggplot(data = res60.df, aes(x=reorder(gtex_tissues_name_to_plot, -utr_gene_freq), y=utr_gene_freq)) +
  geom_bar(width =0.8, stat = "identity", color = "black", lwd = 0.3, fill= "#F5F3F3") +
  theme_bw()  +
  guides(fill = FALSE) +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.30,"cm"),
        legend.title = element_text(size=7),
        legend.text = element_text(size = 7),
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
