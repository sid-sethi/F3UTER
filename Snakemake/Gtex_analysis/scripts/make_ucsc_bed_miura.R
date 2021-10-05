args<-commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
})

miura_dir <- args[1]
outpath <- args[2]

options(stringsAsFactors=F)


files <- list.files(miura_dir, full.names=TRUE, pattern = "list_hg38.bed")

bed <- data.frame()

for(f in files){

  print(str_c(Sys.time(), " - processing ", basename(f)))

  ref <- read.table(f, header=FALSE, sep="\t")
  colnames(ref) <- c("seqnames", "start", "end", "tissue")
  ref <- ref %>% dplyr::filter(tissue %in% "brain")
  ref$seqnames <- as.character(ref$seqnames)
  seqlevelsStyle(ref$seqnames) <- "UCSC"


  ref_mod <- ref %>%
    mutate(
      tissue = str_c(seqnames, start, end, tissue, sep=":"),
      score = 1000,
      strand = ".",
      thickstart = start,
      thickend = end,
      rgb = "153,0,0"
    )

    bed <- rbind(bed, ref_mod)

}


write.table(bed,
          file = paste(outpath, "/miura_brain.bed", sep=""),
          row.names = FALSE, col.names = FALSE, quote = FALSE, sep="\t")
