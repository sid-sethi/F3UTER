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


er <- readRDS(er_file)
er.gr <- er %>% makeGRangesFromDataFrame(., keep.extra.columns = TRUE)


######################################
############ GTF processing ##########
######################################

gtf <- rtracklayer::import(gtf_file)
gtf = keepStandardChromosomes(gtf, species = "Homo_sapiens", pruning.mode="coarse")
gtf.df=as.data.frame(gtf,stringsAsFactor=F)
gtf.df$seqnames <- as.character(gtf.df$seqnames)
seqlevelsStyle(gtf.df$seqnames) <- "UCSC"

##############################################################
################# Extract 3 prime UTRs #######################
##############################################################

utr_all = gtf.df %>% dplyr::filter(type %in% "UTR", ! exon_number %in% c("1", "2"))
utr.gr <- utr_all %>% makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

###################################


hits = IRanges::findOverlaps(er.gr, utr.gr, ignore.strand = TRUE, minoverlap = 0L)
intersect = pintersect(er.gr[queryHits(hits)], utr.gr[subjectHits(hits)])
percentOverlap <- width(intersect) / width(er.gr[queryHits(hits)])
hits <- hits[percentOverlap >= 0.5]


er.hits <- er.gr[unique(queryHits(hits))] %>% as.data.frame()

  df <- data.frame(
    count = nrow(er.hits),
    total = nrow(er)
  ) %>% mutate(
    perc = (count/total)*100,
    dataset = dataset_name
  )


  er.gr.hits <- er.gr[unique(queryHits(hits))]
  utr.hits = mergeByOverlaps(er.gr.hits, utr.gr, ignore.strand = TRUE, minoverlap = 0L) %>% as.data.frame()

  if(nrow(utr.hits)>0){

  #utr.hits$transcript_support_level = str_replace(utr.hits$transcript_support_level, "\\s*\\([^\\)]+\\)","")

  biotype <- table(utr.hits$transcript_type) %>% prop.table() %>% as.data.frame() %>%
    mutate(
      Freq = round(Freq*100, 2),
      dataset = dataset_name
    )

  tsl <- table(utr.hits$transcript_support_level) %>% prop.table() %>% as.data.frame() %>%
    mutate(
      Freq = round(Freq*100, 2),
      dataset = dataset_name
    )

  }else{
    biotype <- data.frame(Var1 = c("nonsense_mediated_decay", "protein_coding"), Freq = 0, dataset= dataset_name)
    tsl <- data.frame(Var1 = c(1:5, "NA"), Freq = 0, dataset= dataset_name)
  }

  write.table(er.hits, str_c(outpath, "/", prefix, "_ers_gencode.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
  write.table(df, str_c(outpath, "/", prefix, "_overlap_gencode.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
  write.table(biotype, str_c(outpath, "/", prefix, "_biotype_gencode.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
  write.table(tsl, str_c(outpath, "/", prefix, "_tsl_gencode.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")


  # matching gene names
  er.df <- er.gr[queryHits(hits)] %>% as.data.frame()
  utr.df <- utr.gr[subjectHits(hits)] %>% as.data.frame() %>% rename_with(~str_c("gtf_", .x)) %>% mutate(gtf_gene_id = str_split_fixed(gtf_gene_id, "\\.", 2)[,1])
  merge.df <- cbind(er.df, utr.df) %>%
    mutate(
      same_gene = ifelse(associated_gene == gtf_gene_id, "yes", "no")
    ) %>%
    distinct(id, .keep_all = "TRUE")
  

  same_gene.df <- table(merge.df$same_gene) %>%
    as.data.frame() %>%
    mutate(total = nrow(merge.df),
    perc = Freq/total*100
    )

  write.table(same_gene.df, str_c(outpath, "/", prefix, "_geneConcordance_gencode.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
