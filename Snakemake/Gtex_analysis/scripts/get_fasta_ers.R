args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(BSgenome)
  library(BSgenome.Hsapiens.UCSC.hg38.masked)
  library(universalmotif)
})


all_tissues <- args[1]
brain_spec <- args[2]
tis_spec <- args[3]
shared <- args[4]
known_utrs <- args[5]
outpath <- args[6]


get_fasta_seq <- function(data){

  # converting to Granges object
  gr = makeGRangesFromDataFrame(data, keep.extra.columns = TRUE)

  # adding "chr" in front of seqnames
  newStyle <- mapSeqlevels(seqlevels(gr), "UCSC")
  gr <- renameSeqlevels(gr, newStyle)

  # Extract sequences
  seq = BSgenome::getSeq(Hsapiens, gr)
  seq@ranges@NAMES <- gr$id

  return(seq)

}


# loading all tissue predictions
all <- readRDS(all_tissues) %>%
  dplyr::select(-c(strand))

all.pos <- all %>%  dplyr::filter(predicted.prob > 0.60)
all.neg <- all %>%  dplyr::filter(predicted.prob <= 0.60)

all.pos.seq <- get_fasta_seq(all.pos)
all.neg.seq <- get_fasta_seq(all.neg)

all.pos.shuf <- shuffle_sequences(all.pos.seq, k = 2, method = "euler", nthreads=2)

writeXStringSet(all.pos.seq, str_c(outpath, "/all.pos.fasta"), format="fasta")
writeXStringSet(all.neg.seq, str_c(outpath, "/all.neg.fasta"), format="fasta")
writeXStringSet(all.pos.shuf, str_c(outpath, "/all.pos.shuf.fasta"), format="fasta")


# brain specific
bs <- readRDS(brain_spec)

bs.pos <- bs %>% dplyr::filter(predicted.prob > 0.60)
bs.neg <- bs %>% dplyr::filter(predicted.prob <= 0.60)

bs.pos.seq <- get_fasta_seq(bs.pos)
bs.neg.seq <- get_fasta_seq(bs.neg)

bs.pos.shuf <- shuffle_sequences(bs.pos.seq, k = 2, method = "euler", nthreads=2)

writeXStringSet(bs.pos.seq, str_c(outpath, "/brainSpec.pos.fasta"), format="fasta")
writeXStringSet(bs.neg.seq, str_c(outpath, "/brainSpec.neg.fasta"), format="fasta")
writeXStringSet(bs.pos.shuf, str_c(outpath, "/brainSpec.pos.shuf.fasta"), format="fasta")

# tissue specific
ts <- readRDS(tis_spec)

ts.pos <- ts %>% dplyr::filter(predicted.prob > 0.60)
ts.neg <- ts %>% dplyr::filter(predicted.prob <= 0.60)

ts.pos.seq <- get_fasta_seq(ts.pos)
ts.neg.seq <- get_fasta_seq(ts.neg)

ts.pos.shuf <- shuffle_sequences(ts.pos.seq, k = 2, method = "euler", nthreads=2)

writeXStringSet(ts.pos.seq, str_c(outpath, "/tisSpec.pos.fasta"), format="fasta")
writeXStringSet(ts.neg.seq, str_c(outpath, "/tisSpec.neg.fasta"), format="fasta")
writeXStringSet(ts.pos.shuf, str_c(outpath, "/tisSpec.pos.shuf.fasta"), format="fasta")

# shared
s <- readRDS(shared)

s.pos <- s %>% dplyr::filter(predicted.prob > 0.60)
s.neg <- s %>% dplyr::filter(predicted.prob <= 0.60)

s.pos.seq <- get_fasta_seq(s.pos)
s.neg.seq <- get_fasta_seq(s.neg)

s.pos.shuf <- shuffle_sequences(s.pos.seq, k = 2, method = "euler", nthreads=2)

writeXStringSet(s.pos.seq, str_c(outpath, "/shared.pos.fasta"), format="fasta")
writeXStringSet(s.neg.seq, str_c(outpath, "/shared.neg.fasta"), format="fasta")
writeXStringSet(s.pos.shuf, str_c(outpath, "/shared.pos.shuf.fasta"), format="fasta")

# known utrs
load(known_utrs)
k <- three_prime %>% plyr::rename(c("three_prime_utr_id" = "id"))

k.seq <- get_fasta_seq(k)

writeXStringSet(k.seq, str_c(outpath, "/known_utrs.fasta"), format="fasta")
