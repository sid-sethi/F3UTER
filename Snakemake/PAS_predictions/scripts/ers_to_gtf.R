args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
})

er_file <- args[1]
resPath <- args[2]
prefix <- args[3]

ers= read.table(er_file, header=TRUE, sep ="\t", stringsAsFactor = FALSE)
seqlevelsStyle(ers$seqnames) <- "NCBI"


# Ensembl format
# er_transcript <- data.frame(
#   seqnames = ers$seqnames,
#   source = "ER",
#   type = "transcript",
#   start = ers$start - 200,
#   end = ers$end,
#   score = ".",
#   strand = "+",
#   tmp1 = ".",
#   gene_id = str_c('gene_id "', ers$ER_id, '"'),
#   transcript_id = str_c('transcript_id "', ers$ER_id, '"'),
#   gene_name = str_c('gene_name "', ers$ER_id, '"'),
#   transcript_name = str_c('transcript_name "', ers$ER_id, '"')
# ) %>%
#   tidyr::unite("unite", gene_id:transcript_name, sep="; ",remove = TRUE) %>%
#   mutate(unite = str_c(unite, ";"))
#
#
# tmp_exon <- data.frame(
#   seqnames = ers$seqnames,
#   source = "ER",
#   type = "exon",
#   start = ers$start - 200,
#   end = ers$start - 100,
#   score = ".",
#   strand = "+",
#   tmp1 = ".",
#   gene_id = str_c('gene_id "', ers$ER_id, '"'),
#   transcript_id = str_c('transcript_id "', ers$ER_id, '"'),
#   exon_number = 'exon_number "1"',
#   gene_name = str_c('gene_name "', ers$ER_id, '"'),
#   transcript_name = str_c('transcript_name "', ers$ER_id, '"')
# ) %>%
#   tidyr::unite("unite", gene_id:transcript_name, sep="; ",remove = TRUE) %>%
#   mutate(unite = str_c(unite, ";"))
#
#
# er_exon <- data.frame(
#   seqnames = ers$seqnames,
#   source = "ER",
#   type = "exon",
#   start = ers$start,
#   end = ers$end,
#   score = ".",
#   strand = "+",
#   tmp1 = ".",
#   gene_id = str_c('gene_id "', ers$ER_id, '"'),
#   transcript_id = str_c('transcript_id "', ers$ER_id, '"'),
#   exon_number = 'exon_number "2"',
#   gene_name = str_c('gene_name "', ers$ER_id, '"'),
#   transcript_name = str_c('transcript_name "', ers$ER_id, '"')
# ) %>%
#   tidyr::unite("unite", gene_id:transcript_name, sep="; ",remove = TRUE) %>%
#   mutate(unite = str_c(unite, ";"))
#



# StringTie format
er_transcript <- data.frame(
  seqnames = ers$seqnames,
  source = "StringTie",
  type = "transcript",
  start = ers$start - 200,
  end = ers$end,
  score = ".",
  strand = "+",
  tmp1 = ".",
  gene_id = str_c('gene_id "STRG.', 1:nrow(ers), '"'),
  transcript_id = str_c('transcript_id "STRG.', 1:nrow(ers), '.', 1, '"'),
  ref_id = str_c('reference_id "', ers$ER_id, '"'),
  ref_gene_id = str_c('ref_gene_id "', ers$ER_id, '"'),
  ref_gene_name = str_c('ref_gene_name "', ers$ER_id, '"')
) %>%
  tidyr::unite("unite", gene_id:ref_gene_name, sep="; ",remove = TRUE) %>%
  mutate(unite = str_c(unite, ";"))


tmp_exon <- data.frame(
  seqnames = ers$seqnames,
  source = "StringTie",
  type = "exon",
  start = ers$start - 200,
  end = ers$start - 100,
  score = ".",
  strand = "+",
  tmp1 = ".",
  gene_id = str_c('gene_id "STRG.', 1:nrow(ers), '"'),
  transcript_id = str_c('transcript_id "STRG.', 1:nrow(ers), '.', 1, '"'),
  exon_number = 'exon_number "1"',
  ref_id = str_c('reference_id "', ers$ER_id, '"'),
  ref_gene_id = str_c('ref_gene_id "', ers$ER_id, '"'),
  ref_gene_name = str_c('ref_gene_name "', ers$ER_id, '"')
) %>%
  tidyr::unite("unite", gene_id:ref_gene_name, sep="; ",remove = TRUE) %>%
  mutate(unite = str_c(unite, ";"))


er_exon <- data.frame(
  seqnames = ers$seqnames,
  source = "StringTie",
  type = "exon",
  start = ers$start,
  end = ers$end,
  score = ".",
  strand = "+",
  tmp1 = ".",
  gene_id = str_c('gene_id "STRG.', 1:nrow(ers), '"'),
  transcript_id = str_c('transcript_id "STRG.', 1:nrow(ers), '.', 1, '"'),
  exon_number = 'exon_number "2"',
  ref_id = str_c('reference_id "', ers$ER_id, '"'),
  ref_gene_id = str_c('ref_gene_id "', ers$ER_id, '"'),
  ref_gene_name = str_c('ref_gene_name "', ers$ER_id, '"')
) %>%
  tidyr::unite("unite", gene_id:ref_gene_name, sep="; ",remove = TRUE) %>%
  mutate(unite = str_c(unite, ";"))


er_gtf <- rbind(er_transcript, tmp_exon, er_exon) %>%
    dplyr::arrange(seqnames, start)


write.table(er_gtf, file = str_c(resPath, "/", prefix, "_ers.gtf"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
