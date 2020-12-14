args <- commandArgs(TRUE) # 1st  = ML output folder

library(tidyverse)
library(RColorBrewer)
library(ggsignif)
library(ggrepel)
options(stringsAsFactors=F)

#mldir = args[1]
mldir = "/ERs/Validation_datasets/Intergenic"

# Plotting the number of GTEx predictions

path1 = str_c(mldir, "/3_prime")
path2 = str_c(mldir, "/5_prime")

files <- list.files(path = path1, pattern = ".txt")
files <- files[!str_detect(files, c("brain_cerebellum", "brain_cortex"))]

if(!dir.exists(paste(mldir, "/Plots", sep=""))){
  system(paste("mkdir -m a=rwx ",mldir, "/Plots", sep=""))
}

# new output path
out.path = paste(mldir, "/Plots", sep="")

print(str_c(Sys.time(), " - merging data from all tissues"))

res3 = data.frame()
res5 = data.frame()
distData = data.frame()
splitread3 = data.frame()
splitread5 = data.frame()

i = 0
for(f in files){
  i= i + 1

  tis <- str_replace(f, ".txt", "") %>% str_replace(., "brain_", "")

  data3 <- read.table(str_c(path1, "/", f), header=TRUE, sep="\t")
  data5 <- read.table(str_c(path2, "/", f), header=TRUE, sep="\t")

  n3 = data3 %>% nrow()
  n5 = data5 %>% nrow()
  total <- n3 + n5

  total.len3 = data3 %>% .$width %>% sum()
  total.len5 = data5 %>% .$width %>% sum()

  res3[i, "tissue"] <- tis
  res3[i, "n"] <- n3
  res3[i, "total_length"] <- total.len3
  res3[i, "total"] <- total
  res3[i, "percent_n"] <- round((n3/total)*100,2)
  res3[i, "group"] <- "3_prime"

  res5[i, "tissue"] <- tis
  res5[i, "n"] <- n5
  res5[i, "total_length"] <- total.len5
  res5[i, "total"] <- total
  res5[i, "percent_n"] <- round((n5/total)*100,2)
  res5[i, "group"] <- "5_prime"


  dist3 = data3 %>% dplyr::select(width, associated_gene, distance_from_geneEnd) %>%
    plyr::rename(c("distance_from_geneEnd" = "distance_from_associated_gene")) %>%
    mutate(tissue = tis, group = "3_prime")

  dist5 = data5 %>% dplyr::select(width, associated_gene, distance_from_geneTSS) %>%
    plyr::rename(c("distance_from_geneTSS" = "distance_from_associated_gene")) %>%
    mutate(tissue = tis, group = "5_prime")

  dist = rbind(dist3, dist5)
  distData = rbind(distData, dist)


  sr3 <- data3 %>% filter(annotationType_split_read_annot %in% "partially annotated split read") %>% nrow()
  no.sr3 <- n3 - sr3
  sr5 <- data5 %>% filter(annotationType_split_read_annot %in% "partially annotated split read") %>% nrow()
  no.sr5 <- n5 - sr5
  splitread3[i, "tissue"] <- tis
  splitread3[i, "SR"] <- sr3
  splitread3[i, "NOSR"] <- no.sr3
  splitread3[i, "group"] <- "3_prime"
  splitread5[i, "tissue"] <- tis
  splitread5[i, "SR"] <- sr5
  splitread5[i, "NOSR"] <- no.sr5
  splitread5[i, "group"] <- "5_prime"

}


splitread <- rbind(splitread3, splitread5)
# % of ERs with split read
total_ers <- sum(splitread$SR) + sum(splitread$NOSR)
sr_perc <- (sum(splitread$SR)/total_ers)*100


res <- rbind(res3, res5)


print(str_c(Sys.time(), " - plotting ....."))

formatting <- read_delim("OMIM_gtex_tissue_name_formatting.csv", delim = ",") %>% mutate(tissue_color_hex = str_c("#", tissue_color_hex)) %>%
  drop_na() %>%
  dplyr::select(gtex_tissues_name_formatted_2, gtex_tissue_group, gtex_tissues_name_to_plot, tissue_color_hex, tissue_type)


res = dplyr::left_join(res, formatting, by = c("tissue" = "gtex_tissues_name_formatted_2"))


##### 3' Vs 5'
man.pval = wilcox.test(res3$percent_n,res5$percent_n)$p.value
label.pval = paste("p = ",format.pval(man.pval,3),sep="")

png(str_c(out.path, "/3_prime_vs_5_prime.png"), height = 3.75, width = 2.25, res = 600, units = "in")
ggplot(data=res, aes(x=group, y=percent_n)) +
  geom_boxplot(fill="#E0E0E0", lwd = 0.55, outlier.shape=NA, alpha=0.8) +
  #geom_jitter(position=position_jitter(0.28), shape=21, size=1, alpha=0.7, fill="red") +
  geom_jitter(aes(fill=tissue_type), position=position_jitter(0.28), shape=21, size=1, alpha=0.7) +
  scale_fill_manual(values = c("#EEEE00", "#CCFFFF")) +
  theme_bw()+
  guides(fill=FALSE, color=FALSE)+
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.40,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        legend.position= c(0.83,0.91),
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
  theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_x_discrete(name = "Intergenic ERs\nper tissue")+
  scale_y_continuous(name = "% intergenic ERs", limits = c(0,80)) +
  geom_signif(comparisons = list(c("3_prime", "5_prime")),map_signif_level=TRUE, annotation = label.pval, size=0.4,textsize=3,tip_length=0.02, margin_top = 0.1, color = "black", fontface = 2)
dev.off()

### total length - brain vs non-brain

res2 <- res %>%
  group_by(tissue, tissue_type, gtex_tissues_name_to_plot) %>%
  summarise(n = sum(n), total_length = sum(total_length))

brain <- res2 %>% filter(tissue_type %in% "Brain")
nbrain <- res2 %>% filter(tissue_type %in% "Non-brain")

man.pval = wilcox.test(brain$total_length,nbrain$total_length)$p.value
label.pval = paste("p = ",format.pval(man.pval,3),sep="")

png(str_c(out.path, "/total_length_brain_vs_nonbrain.png"), height = 3.75, width = 2.25, res = 600, units = "in")
ggplot(data=res2, aes(x=tissue_type, y=total_length/1000)) +
  geom_boxplot(fill="grey", lwd = 0.55, outlier.shape=NA, alpha=0.8) +
  geom_jitter(position=position_jitter(0.20), shape=21, size=1, alpha=0.7, fill="red") +
  theme_bw()+
  guides(fill=FALSE, color=FALSE)+
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.40,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        legend.position= c(0.83,0.91),
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
  theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_x_discrete(name = "Tissue type")+
  scale_y_continuous(name = "Total unannotated intergenic ER length (kb)", limits = c(700, 3500)) +
  geom_signif(comparisons = list(c("Brain", "Non-brain")),map_signif_level=TRUE, annotation = label.pval, size=0.4,textsize=3,tip_length=0.02, margin_top = 0.1, color = "black", fontface = 2)
dev.off()


##### numbers - brain vs non-brain

man.pval = wilcox.test(brain$n,nbrain$n)$p.value
label.pval = paste("p = ",format.pval(man.pval,3),sep="")

png(str_c(out.path, "/numbers_brain_vs_nonbrain.png"), height = 3.75, width = 2.25, res = 600, units = "in")
ggplot(data=res2, aes(x=tissue_type, y=n)) +
  geom_boxplot(fill="grey", lwd = 0.55, outlier.shape=NA, alpha=0.8) +
  geom_jitter(position=position_jitter(0.20), shape=21, size=1, alpha=0.7, fill="red") +
  theme_bw()+
  guides(fill=FALSE, color=FALSE)+
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.40,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        legend.position= c(0.83,0.91),
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
  theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_x_discrete(name = "Tissue type")+
  scale_y_continuous(name = "Number of unannotated intergenic ERs", limits = c(2500,13500)) +
  geom_signif(comparisons = list(c("Brain", "Non-brain")),map_signif_level=TRUE, annotation = label.pval, size=0.4,textsize=3,tip_length=0.02, margin_top = 0.1, color = "black", fontface = 2)
dev.off()


####### scatter plot - Numbers vs length #############
png(str_c(out.path, "/scatter_numbers_vs_totalLength.png"), height = 3.75, width = 5.25, res = 600, units = "in")
ggplot(data=res2, aes(x=n, y=total_length/1000, fill = tissue_type)) +
  geom_point(size=2, shape=21) +
  theme_bw()+
  geom_text_repel(aes(label=gtex_tissues_name_to_plot), size=1.8, segment.size=0.2, segment.alpha = 0.5) +
  scale_fill_manual(values = c("#EEEE00", "#CCFFFF")) +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.35,"cm"),
        legend.title = element_text(size=8.5),
        legend.text = element_text(size = 8),
        legend.position= c(0.90,0.10),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size = 9),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        panel.border = element_blank(),
        axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=11,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_x_continuous(name = "Number of intergenic ERs")+
  scale_y_continuous(name = "Total genomic space covered (kb)")
dev.off()
