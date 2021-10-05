args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
})

options(stringsAsFactors = FALSE)

res_dir <- args[1]
outpath = args[2]
formatting_csv <- args[3]

all_pos_files <- list.files(res_dir, pattern="^all.pos_ers", recursive = TRUE, full.names = TRUE)

all_pos_files <- all_pos_files[-2]


data_all = data.frame()
for(f in all_pos_files){
  x = read.table(f, header=TRUE, sep="\t")
  x = x %>% dplyr::select(id, width, associated_gene, tissue)
  data_all = bind_rows(data_all, x)
}
data_all = data_all %>% distinct(id, .keep_all = TRUE)


data_ref = data.frame()
for(f in all_pos_files[-c(3, 4, 5)]){
  x = read.table(f, header=TRUE, sep="\t")
  x = x %>% dplyr::select(id, width, associated_gene)
  data_ref = bind_rows(data_ref, x)
}
data_ref = data_ref %>% distinct(id, .keep_all = TRUE)


data_miura = data.frame()
for(f in all_pos_files[c(3:5)]){
  x = read.table(f, header=TRUE, sep="\t")
  x = x %>% dplyr::select(id, width, associated_gene)
  data_miura = bind_rows(data_miura, x)
}
data_miura = data_miura %>% distinct(id, .keep_all = TRUE)


cat(
  "Ref: number of distinct overlapping ERs = ", nrow(data_ref), "\n",
    "Ref: total genomic space of overlapping ERs (kb) = ", data_ref$width %>% sum()/1000, "\n",
    "Ref: number of distinct associated genes = ", data_ref$associated_gene %>% unique() %>% length(), "\n",
  "Miura: number of distinct overlapping ERs = ", nrow(data_miura), "\n",
  "Miura: total genomic space of overlapping ERs (kb) = ", data_miura$width %>% sum()/1000, "\n",
  "Miura: number of distinct associated genes = ", data_miura$associated_gene %>% unique() %>% length(), "\n",
  "All: number of distinct overlapping ERs = ", nrow(data_all), "\n",
  "All: total genomic space of overlapping ERs (kb) = ", data_all$width %>% sum()/1000, "\n",
  "All: number of distinct associated genes = ", data_all$associated_gene %>% unique() %>% length(),
  file = str_c(outpath, "/overall_numbers.txt")
    )


### numbers all together per tissue #######
tis <- data_all %>%
  group_by(tissue) %>%
  summarise(count = n(),
      total_width = sum(width),
      .groups = "keep"
  ) %>%
  as.data.frame()


  cat(
    "Average number of 3'UTRs per tissue previously detected = ", mean(tis$count), "\n",
    "Average length of 3'UTRs per tissue previously detected (kb) = ", mean(tis$total_width)/1000, "\n",
    "Minimum length of 3'UTRs per tissue previously detected (kb) = ", min(tis$total_width)/1000, "\n",
    "Maximum length of 3'UTRs per tissue previously detected (kb) = ", max(tis$total_width)/1000, "\n",
    file = str_c(outpath, "/per_tissue_numbers.txt")
  )


  print(str_c(Sys.time(), " - plotting ....."))

  formatting <- read_delim(formatting_csv, delim = ",") %>% mutate(tissue_color_hex = str_c("#", tissue_color_hex)) %>%
    drop_na() %>%
    dplyr::select(OMIM_gtex_name, gtex_tissue_group, gtex_tissues_name_to_plot, tissue_color_hex)


  res = dplyr::left_join(tis, formatting, by = c("tissue" = "OMIM_gtex_name"))


  col.df <- data.frame(tissue = levels(reorder(res$gtex_tissues_name_to_plot, -res$count))) %>%
    dplyr::left_join(formatting, by = c("tissue" = "gtex_tissues_name_to_plot"))


  png(str_c(outpath, "/number_annotated_predictions_perTissue.png"), height = 4.20, width = 6.75, res = 600, units = "in")
  ggplot(data = res, aes(x=reorder(gtex_tissues_name_to_plot, -count), y=count)) +
    geom_bar(width =0.8, stat = "identity", color = "black", lwd = 0.5, fill = "#99CCFF") +
    theme_bw()  +
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
    scale_y_continuous(name = "# of 3'UTR predictions detected in reference\ngene annotations and Miura et al. data") +
    geom_point(data = col.df, # object that has the order of the tissues to be plotted
               aes(x = tissue, y = -1.5), # set
               colour = "black", fill = col.df$tissue_color_hex, shape = 22, size = 2)+
    coord_cartesian(ylim = c(0, 60), clip = "off")
  dev.off()


  # length of predictions

  col.df <- data.frame(tissue = levels(reorder(res$gtex_tissues_name_to_plot, -res$total_width))) %>%
    dplyr::left_join(formatting, by = c("tissue" = "gtex_tissues_name_to_plot"))

  png(str_c(outpath, "/length_annotated_predictions_perTissue.png"), height = 4.20, width = 6.75, res = 600, units = "in")
  ggplot(data = res, aes(x=reorder(gtex_tissues_name_to_plot, -total_width), y=total_width/1000)) +
    geom_bar(width =0.8, stat = "identity", color = "black", lwd = 0.5, fill = "#99CCFF") +
    theme_bw()  +
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
    scale_y_continuous(name = "Total length (in kb) of 3'UTR predictions detected\nin reference gene annotations and Miura et al.") +
    geom_point(data = col.df, # object that has the order of the tissues to be plotted
               aes(x = tissue, y = -1.2), # set
               colour = "black", fill = col.df$tissue_color_hex, shape = 22, size = 2)+
    coord_cartesian(ylim = c(0, 42), clip = "off")
  dev.off()





all_pos_files <- list.files(res_dir, pattern="^all.pos_overlap", recursive = TRUE, full.names = TRUE)

all_pos_files <- all_pos_files[-c(4:6)]

all_pos <- all_pos_files %>%
  map_dfr(read.table, header=TRUE) %>%
  mutate(group = "Predicted 3-UTRs")



df.overlap <- bind_rows(all_pos)
df.overlap$dataset <- fct_relevel(df.overlap$dataset, "Ensembl_v92", "Ensembl_v104", "Gencode_v38", "Refseq_v109")

png(str_c(outpath, "/all_ref_overlap.png"),bg="transparent",units="in", width = 4.25, height= 3.75 ,res=600)
ggplot(data = df.overlap, aes(x=dataset, y=perc)) +
  geom_bar(width =0.6, stat = "identity", color = "black", lwd = 0.3, position = "dodge", fill = "#00AFBB") +
  theme_bw()  +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.40,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position= c(0.25,0.95),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=11, angle=45, hjust=1),
        axis.text.y = element_text(size = 10),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=11, angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_y_continuous(name = "% predicted 3'UTRs overlapping\nreference annotation")
dev.off()

write.table(df.overlap, str_c(outpath, "/all_ref_overlap.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")


all_pos_files <- list.files(res_dir, pattern="^all.pos_biotype", recursive = TRUE, full.names = TRUE)

all_pos_files <- all_pos_files[-c(2)]

all_pos <- all_pos_files %>%
  map_dfr(read.table, header=TRUE) %>%
  mutate(group = "Predicted 3-UTRs")


df.biotype <- bind_rows(all_pos)

png(str_c(outpath, "/all_ref_biotype.png"),bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
ggplot(data = df.biotype, aes(x=dataset, y=Freq, fill = Var1)) +
  geom_bar(width = 0.6, stat = "identity", color = "black", lwd = 0.3) +
  theme_bw()  +
  scale_fill_brewer(palette = c("RdYlBu")) +
  guides(fill = guide_legend(title = "Transcript biotype", nrow = 1, title.position = "top")) +
  #facet_wrap(~dataset, nrow=1) +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.35,"cm"),
        legend.title = element_text(size=10),
        legend.text = element_text(size = 10),
        legend.position= "top",
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=11, angle=45, hjust=1),
        axis.text.y = element_text(size = 10),
        panel.border = element_rect(colour="BLACK",size=0.4),
        #panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=11,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA),
        strip.text.x = element_text(size=10, face="bold"),
        strip.background = element_rect(fill="white")
  )+
  scale_y_continuous(name = "% predicted 3'UTRs overlapping\nreference annotation")
dev.off()

write.table(df.biotype, str_c(outpath, "/all_ref_biotype.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")



all_pos_files <- list.files(res_dir, pattern="^all.pos_tsl", recursive = TRUE, full.names = TRUE)

all_pos_files <- all_pos_files[-c(2)]

all_pos <- all_pos_files %>%
  map_dfr(read.table, header=TRUE) %>%
  mutate(group = "Predicted 3-UTRs")


df.tsl <- bind_rows(all_pos)
df.tsl$Var1 <- as.factor(df.tsl$Var1)


png(str_c(outpath, "/all_ref_tsl.png"), bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
ggplot(data = df.tsl, aes(x=dataset, y=Freq, fill = Var1)) +
  geom_bar(width =0.5, stat = "identity", color = "black", lwd = 0.3) +
  theme_bw()  +
  scale_fill_brewer(palette = c("RdYlBu")) +
  #guides(fill = guide_legend(title = "Transcript support level", nrow = 1, title.position = "left")) +
  guides(fill = guide_legend(title = "Transcript support level", nrow = 6, title.position = "top")) +
  #facet_wrap(~dataset, nrow=1) +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.35,"cm"),
        legend.title = element_text(size=10),
        legend.text = element_text(size = 10),
        legend.position= "right",
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=11, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        panel.border = element_rect(colour="BLACK",size=0.4),
        #panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=11,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA),
        strip.text.x = element_text(size=10, face="bold"),
        strip.background = element_rect(fill="white")
  )+
  scale_y_continuous(name = "% predicted 3'UTRs overlapping\nreference annotation")
dev.off()

write.table(df.tsl, str_c(outpath, "/all_ref_tsl.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")



all_pos_files <- list.files(res_dir, pattern="^all.pos_overlap", recursive = TRUE, full.names = TRUE) %>% str_subset("miura")

all_pos <- all_pos_files %>%
  map_dfr(read.table, header=TRUE) %>%
  mutate(group = "Predicted 3-UTRs")


df.overlap <- bind_rows(all_pos)
df.overlap$dataset <- str_to_title(df.overlap$dataset)


png(str_c(outpath, "/all_miura_overlap.png"),bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
ggplot(data = df.overlap, aes(x=dataset, y=perc)) +
  geom_bar(width =0.6, stat = "identity", color = "black", lwd = 0.3, position = "dodge", fill= "#00AFBB") +
  theme_bw()  +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  #coord_cartesian(ylim = c(0,1.5)) +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.35,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        legend.position= c(0.82,0.95),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size = 9),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=11,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_y_continuous(name = "% predicted 3'UTRs overlapping\nMiura et al. annotation")
dev.off()

write.table(df.overlap, str_c(outpath, "/all_miura_overlap.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
