library(tidyverse)
options(stringsAsFactors = FALSE)

# load ATTRACT database file
x = read.table("/Attract/ATtRACT_db.txt", header=TRUE, sep="\t")

# filter motifs
h = x %>% dplyr::filter(Organism %in% "Homo_sapiens", Len > 6, Score %in% "1.000000**") %>%
  group_by(Gene_id) %>%
  slice(which.max(Len))

write.table(h, "/Attract/ATtRACT_db_human.txt", quote = FALSE, row.names=FALSE, sep="\t")
