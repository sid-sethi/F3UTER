args<-commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
})

options(stringsAsFactors=F)


erFile = args[1]
outPath = args[2]
tissue = args[3]
repeats_data = args[4] # repeatmasker file

tissues = c(tissue)


# repeats data
d = read.table(repeats_data, header=FALSE)
colnames(d) = c("seqnames", "start", "end", "repeat_type")

repeats = d %>% dplyr::filter(!repeat_type %in% c("RNA", "rRNA", "scRNA", "snRNA", "srpRNA", "tRNA", "Unknown", "Simple_repeat", "Satellite/telo", "Low_complexity", "Satellite", "Satellite/acro", "Satellite/centr"))

rep.gr = makeGRangesFromDataFrame(repeats, keep.extra.columns = TRUE)
rep.gr = keepStandardChromosomes(rep.gr, species = "Homo_sapiens", pruning.mode="coarse")



for(tissue in tissues){

  print(tissue)

  #ER data
  er = read.table(erFile, header=T, sep="\t")

  # converting to Granges object
  er.gr = makeGRangesFromDataFrame(er, keep.extra.columns = TRUE)

  # overlap
  hits = IRanges::findOverlaps(er.gr, rep.gr, ignore.strand = TRUE)
  ovr = pintersect(er.gr[queryHits(hits)], rep.gr[subjectHits(hits)])
  ovr$fracOverlap <- round( width(ovr) / width(er.gr[queryHits(hits)]) , 3)

  overlaps <- ovr %>% as.data.frame() %>%
    group_by(ER_id) %>%
    summarise(fracOverlap = sum(fracOverlap))

  er.export <- dplyr::left_join(er, overlaps, by = "ER_id") %>%
    mutate(fracOverlap = replace_na(fracOverlap, 0)) %>%
    mutate(fracOverlap = ifelse(fracOverlap > 1, 1, fracOverlap)) %>%
    mutate(tissue = tissue) %>%
    dplyr::select(c(ER_id, fracOverlap, tissue)) %>%
    plyr::rename(c("ER_id" = "id"))

  file_name = str_c(outPath, "/", tissue, "_repeats.txt")
  write.table(er.export, file = file_name , row.names = FALSE, quote = FALSE, sep="\t")

  rm(overlaps, er.export, hits, er.gr, ovr)

}
