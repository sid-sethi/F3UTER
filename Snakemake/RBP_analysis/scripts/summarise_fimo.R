args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(Biostrings)
})
options(stringsAsFactors = FALSE)

fimo_file <- args[1]
pval <- as.numeric(args[2])
outpath <- args[3]
prefix <- args[4]
fasta <- args[5]


seq = readDNAStringSet(fasta)
all_ids = data.frame(sequence_name = names(seq))

x = read.table(fimo_file, header=TRUE) %>%
  mutate(motif_len = (stop-start)+1) %>%
  dplyr::filter(p.value < pval)

x.group <- x %>%
  group_by(sequence_name) %>%
  summarise(count = n(), .groups = "keep")


x.group.final = dplyr::left_join(all_ids, x.group, by = "sequence_name") %>%
  replace_na(list(count = 0))


x.df <- x.group.final %>%
  separate(sequence_name, c(NA, NA, "start", "end", NA), sep=":", remove = FALSE, convert = TRUE) %>%
  mutate(sequence_len = (end-start)+1,
         #enrich_score = count/((sequence_len - motif_len)+1)) %>%
         enrich_score = count/(sequence_len/100)) %>%
  dplyr::select(-c("start", "end"))


write.table(x.df, file = str_c(outpath, "/", prefix, "_seq_sum.txt"), quote=FALSE, row.names=FALSE, sep="\t")
