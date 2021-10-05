args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
})

options(stringsAsFactors=F)

res_dir <- args[1]
gtf_file <- args[2]
outpath <- args[3]



gtf <- rtracklayer::import(gtf_file)
gtf = keepStandardChromosomes(gtf, species = "Homo_sapiens", pruning.mode="coarse")
gtf.df=as.data.frame(gtf,stringsAsFactor=F)
gtf.df.pc = gtf.df %>% dplyr::filter(gene_biotype %in% "protein_coding", transcript_biotype %in% "protein_coding")

# Extract 3 prime UTRs
utr.pc = gtf.df.pc %>% dplyr::filter(type %in% "three_prime_utr")
utr.pc = utr.pc %>% mutate(length = end - start)

utr_per_transcript = utr.pc %>%
  group_by(gene_id, transcript_id) %>%
  summarize(maximal_3utr_length = sum(length), .groups = "keep") %>%
  as.data.frame()


utr_per_gene = utr_per_transcript %>%
  group_by(gene_id) %>%
  summarize(max_3utr_length = max(maximal_3utr_length), .groups = "keep") %>%
  as.data.frame()


# loading predictions
bs.pos <- readRDS(str_c(res_dir, "/brainSpec.pos.rds")) %>%
  dplyr::left_join(utr_per_gene, by = c("associated_gene" = "gene_id")) %>%
  mutate(er_vs_utr_length = round(width/max_3utr_length, 2), group = "Highly brain-specific") %>%
  dplyr::select(associated_gene, width, max_3utr_length, er_vs_utr_length, group)


ts.pos <- readRDS(str_c(res_dir, "/tisSpec.pos.rds")) %>%
  dplyr::left_join(utr_per_gene, by = c("associated_gene" = "gene_id")) %>%
  mutate(er_vs_utr_length = round(width/max_3utr_length, 2), group = "Absolute tissue-specific") %>%
  dplyr::select(associated_gene, width, max_3utr_length, er_vs_utr_length, group)


s.pos <- readRDS(str_c(res_dir, "/shared.pos.rds")) %>%
  dplyr::left_join(utr_per_gene, by = c("associated_gene" = "gene_id")) %>%
  mutate(er_vs_utr_length = round(width/max_3utr_length, 2), group = "Shared") %>%
  dplyr::select(associated_gene, width, max_3utr_length, er_vs_utr_length, group)



cat(
  "Average length of brain-specific 3'UTR predictions: ", mean(bs.pos$width), "\n",
  "Average length fold change of brain-specific predictions vs 3'UTR maximal length of its associated genes: ", mean(bs.pos$er_vs_utr_length, na.rm=TRUE), "\n",
  "Average length of tissue-specific 3'UTR predictions: ", mean(ts.pos$width), "\n",
  "Average length fold change of tissue-specific predictions vs 3'UTR maximal length of its associated genes: ", mean(ts.pos$er_vs_utr_length, na.rm=TRUE), "\n",
  "Average length of shared 3'UTR predictions: ", mean(s.pos$width), "\n",
  "Average length fold change of shared predictions vs 3'UTR maximal length of its associated genes: ", mean(s.pos$er_vs_utr_length, na.rm=TRUE), "\n",
  file = str_c(outpath, "/predicted_threePrime_vs_utr_length_FC.txt")
)


df <- rbind(bs.pos, ts.pos, s.pos)


file1 = str_c(outpath, "/predicted_threePrime_vs_utr_length_ecdf.png")
png(file1,bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
ggplot(df, aes(er_vs_utr_length, colour = group)) +
  stat_ecdf(lwd = 0.6) +
  scale_color_manual(values= c("#e41a1c", "#377eb8", "#4daf4a")) +
  theme_bw()  +
  labs(color = "Predicted 3'UTRs")+
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.35,"cm"),
        legend.title = element_text(size=10),
        legend.text = element_text(size = 10),
        legend.position= c(0.77,0.17),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size = 9),
        panel.border = element_rect(colour="BLACK",size=0.4),
        #panel.border = element_blank(),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA),
        strip.background = element_rect(fill="white"),
        strip.text.x = element_text(face = "bold", size=9)
  )+
  #theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_x_continuous(name = "Fold change (3'UTR prediction length vs\nmaximal 3'UTR length of gene)", trans="sqrt")+
  scale_y_continuous(name = "CDF")
dev.off()
