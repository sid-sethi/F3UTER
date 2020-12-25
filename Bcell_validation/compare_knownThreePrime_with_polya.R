library(tidyverse)

# loading the known three prime utr regions
load("/GS.RData") # utr

# converting to Granges object
utr.gr = makeGRangesFromDataFrame(utr, keep.extra.columns = TRUE)
newStyle <- mapSeqlevels(seqlevels(utr.gr), "UCSC")
utr.gr <- renameSeqlevels(utr.gr, newStyle)


print(str_c(Sys.time(), " - processing GCB1 poly(A) data"))

# polya site data
pa = read.table("/GSE111310/GCB1/Polya_clusters/GCB1.clusters.hg38.2-0.bed",header=F,sep="\t")
colnames(pa) = c("chr","start","end","id","exp_thisSample","strand", "perc_sample_support", "number_protocol_support", "avg_exp_allSamples", "cluster_annotation","pas")
pac = makeGRangesFromDataFrame(pa, keep.extra.columns = TRUE)
pac = keepStandardChromosomes(pac, species = "Homo_sapiens", pruning.mode="coarse")
newStyle <- mapSeqlevels(seqlevels(pac), "UCSC")
pac <- renameSeqlevels(pac, newStyle)

utr.ov = IRanges::subsetByOverlaps(utr.gr, pac, ignore.strand = FALSE, maxgap= -1, type = "any") %>% length()
utr.ov.perc <- round((utr.ov/length(utr.gr))*100,2)
# 87.9%


print(str_c(Sys.time(), " - processing GCB2 poly(A) data"))

# polya site data
pa2 = read.table("/GSE111310/GCB2/Polya_clusters/GCB2.clusters.hg38.2-0.bed",header=F,sep="\t")
colnames(pa2) = c("chr","start","end","id","exp_thisSample","strand", "perc_sample_support", "number_protocol_support", "avg_exp_allSamples", "cluster_annotation","pas")
pac2 = makeGRangesFromDataFrame(pa2, keep.extra.columns = TRUE)
pac2 = keepStandardChromosomes(pac2, species = "Homo_sapiens", pruning.mode="coarse")
newStyle <- mapSeqlevels(seqlevels(pac2), "UCSC")
pac2 <- renameSeqlevels(pac2, newStyle)

utr.ov = IRanges::subsetByOverlaps(utr.gr, pac2, ignore.strand = FALSE, maxgap= -1, type = "any") %>% length()
utr.ov.perc <- round((utr.ov/length(utr.gr))*100,2)
# 87.9%
