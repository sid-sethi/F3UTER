args<-commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
})

er_data <- args[1]
outpath <- args[2]
tissue <- args[3]

options(stringsAsFactors=F)

  
print(str_c(Sys.time(), " - processing: ", tissue))

er = readRDS(er_data)

bed = er %>% dplyr::select(seqnames, start, end, id, predicted.prob) %>%
mutate(
  strand = ".",
  thickstart = start,
  thickend = end,
  rgb = ifelse(predicted.prob > 0.60, "0,175,187", "231,184,0")
)

write.table(bed,
          file = paste(outpath, "/", tissue, ".bed", sep=""),
          row.names = FALSE, col.names = FALSE, quote = FALSE, sep="\t")




## trackline row to add afterwards
## track name="ItemRGBDemo" visibility=1 itemRgb="On"
