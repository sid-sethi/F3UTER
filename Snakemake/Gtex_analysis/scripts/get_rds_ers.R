args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
})


all_tissues <- args[1]
brain_spec <- args[2]
tis_spec <- args[3]
shared <- args[4]
known_utrs <- args[5]
utr_tissue_specificity <- args[6]
outpath <- args[7]


# loading all tissue predictions
all <- readRDS(all_tissues)

all.pos <- all %>%  dplyr::filter(predicted.prob > 0.60)
all.neg <- all %>%  dplyr::filter(predicted.prob <= 0.60)

saveRDS(all.pos, file = str_c(outpath, "/all.pos.rds"))
saveRDS(all.neg, file = str_c(outpath, "/all.neg.rds"))


# brain specific
bs <- readRDS(brain_spec)

bs.pos <- bs %>% dplyr::filter(predicted.prob > 0.60)
bs.neg <- bs %>% dplyr::filter(predicted.prob <= 0.60)

saveRDS(bs.pos, file = str_c(outpath, "/brainSpec.pos.rds"))
saveRDS(bs.neg, file = str_c(outpath, "/brainSpec.neg.rds"))

# tissue specific
ts <- readRDS(tis_spec)

ts.pos <- ts %>% dplyr::filter(predicted.prob > 0.60)
ts.neg <- ts %>% dplyr::filter(predicted.prob <= 0.60)

saveRDS(ts.pos, file = str_c(outpath, "/tisSpec.pos.rds"))
saveRDS(ts.neg, file = str_c(outpath, "/tisSpec.neg.rds"))

# shared
s <- readRDS(shared)

s.pos <- s %>% dplyr::filter(predicted.prob > 0.60)
s.neg <- s %>% dplyr::filter(predicted.prob <= 0.60)

saveRDS(s.pos, file = str_c(outpath, "/shared.pos.rds"))
saveRDS(s.neg, file = str_c(outpath, "/shared.neg.rds"))

# known utrs
load(known_utrs)
utrs <- three_prime %>% plyr::rename(c("three_prime_utr_id" = "id"))

utr_ts <- read.table(utr_tissue_specificity, header=TRUE, sep="\t")
utrs <- dplyr::left_join(utrs, utr_ts, by = "id")

utr_all <- utrs %>% dplyr::filter(! outcome %in% "Not expressed")

saveRDS(utr_all, file = str_c(outpath, "/known_utrs.rds"))
