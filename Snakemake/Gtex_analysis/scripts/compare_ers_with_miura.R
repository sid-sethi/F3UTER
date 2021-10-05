args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
})

options(stringsAsFactor = FALSE, warn=-1)

gtf_file <- args[1]
er_file <- args[2]
outpath <- args[3]
prefix <- args[4]
dataset_name <- args[5]
tissue_mapping <- args[6]

mapping = read_delim(tissue_mapping, delim = ",") %>%
  drop_na()


er <- readRDS(er_file)
er <- left_join(er, mapping, by = c("tissue" =  "OMIM_gtex_name"))

er.gr <- er %>% makeGRangesFromDataFrame(., keep.extra.columns = TRUE)


######################################
############ Refseq processing ##########
######################################

ref <- read.table(gtf_file, header=FALSE, sep="\t")

# ref_mod <- ref %>% tidyr::separate(new_annotation, c("seqnames", "cord"), sep=":", remove=TRUE) %>%
#   tidyr::separate(cord, c("start", "end"), sep="-", remove=TRUE, convert =TRUE)

colnames(ref) <- c("seqnames", "start", "end", "tissue")

ref <- ref %>% dplyr::filter(!tissue %in% c("breast", "lymph_node", "ovary", "prostate", "testes"))

ref$seqnames <- as.character(ref$seqnames)
seqlevelsStyle(ref$seqnames) <- "UCSC"

utr.gr <- ref %>% makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

###################################


hits = IRanges::findOverlaps(er.gr, utr.gr, ignore.strand = TRUE, minoverlap = 0L)
intersect = pintersect(er.gr[queryHits(hits)], utr.gr[subjectHits(hits)])
percentOverlap <- width(intersect) / width(er.gr[queryHits(hits)])
hits <- hits[percentOverlap >= 0.5]

er.df <- er.gr[queryHits(hits)] %>% as.data.frame()
utr.df <- utr.gr[subjectHits(hits)] %>% as.data.frame() %>% rename_with(~str_c("miura_", .x))
merge.df <- cbind(er.df, utr.df) %>%
  mutate(
    same_tissue = ifelse(miura_tissue_group == miura_tissue, "yes", "no")
  ) %>%
  dplyr::filter(same_tissue %in% "yes")


er.hits <- merge.df %>% distinct(id, .keep_all = TRUE)
#er.hits <- er.gr[unique(queryHits(hits))] %>% as.data.frame()

  df <- data.frame(
    count = nrow(er.hits),
    total = nrow(er)
  ) %>% mutate(
    perc = (count/total)*100,
    dataset = dataset_name
  )


  write.table(er.hits, str_c(outpath, "/", prefix, "_ers_miura.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
  write.table(df, str_c(outpath, "/", prefix, "_overlap_miura.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
