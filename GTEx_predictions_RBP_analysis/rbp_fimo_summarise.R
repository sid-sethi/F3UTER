args <- commandArgs(TRUE) # 1 = fimo output file (.tsv); 2 = output file prefix; 3 = pval threshold
library(tidyverse)
library(matrixStats)
options(stringsAsFactors = FALSE)

i = args[1]
out = args[2]
pval <- as.numeric(args[3])

# read FIMO data
x = read.table(i, header=TRUE) %>%
  mutate(motif_len = (stop-start)+1) %>%
  filter(p.value < pval)

# group results by per sequence
x.group <- x %>%
  group_by(sequence_name, motif_id) %>%
  summarise(count = n(),
            motif_len = unique(motif_len))

# calculate enrichment score per sequence
x.df <- x.group %>%
  separate(sequence_name, c(NA, NA, "start", "end", NA, NA), sep=":", remove = FALSE, convert = TRUE) %>%
  mutate(sequence_len = (end-start)+1,
         enrich_score = count/(sequence_len/100)) %>%
  dplyr::select(-c("start", "end"))


x.wide <- pivot_wider(x.df, id_cols = sequence_name, names_from = motif_id, values_from = enrich_score)


x.wide$enrich_sum <- rowSums(x.wide[2:ncol(x.wide)], na.rm=TRUE)
x.rowSum <- x.wide %>%
  dplyr::select(sequence_name, enrich_sum)

write.table(x.rowSum, file = str_c(out,"_seq_sum.txt"), quote=FALSE, row.names=FALSE, sep="\t")
