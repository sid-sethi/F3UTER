args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(DOSE)
  library(org.Hs.eg.db)
})

options(stringsAsFactors=F)

pos_list <- args[1]
neg_list <- args[2]
universe <- args[3]
outpath <- args[4]

pos <- read.table(pos_list, header=TRUE, sep="\t")
neg <- read.table(neg_list, header=TRUE, sep="\t")
uni <- read.table(universe, header=TRUE, sep="\t")

pos.entrez <- bitr(pos$associated_gene, fromType="ENSEMBL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
neg.entrez <- bitr(neg$associated_gene, fromType="ENSEMBL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
uni.entrez <- bitr(uni$Ensembl.gene.ID, fromType="ENSEMBL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")

pos_prefix <- basename(pos_list) %>% str_replace(., ".txt", "")
neg_prefix <- basename(neg_list) %>% str_replace(., ".txt", "")

pos.go <- enrichGO(gene          = pos.entrez$ENTREZID,
                universe      = uni.entrez$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 1,
                readable      = TRUE)

pos.go.res <- pos.go@result
write.table(pos.go.res, str_c(outpath, "/", pos_prefix, "_GO.txt"), row.names=FALSE, quote=FALSE, sep="\t")

neg.go <- enrichGO(gene          = neg.entrez$ENTREZID,
                universe      = uni.entrez$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 1,
                readable      = TRUE)

neg.go.res <- neg.go@result
write.table(neg.go.res, str_c(outpath, "/", neg_prefix, "_GO.txt"), row.names=FALSE, quote=FALSE, sep="\t")


list = list(Predicted_3_UTRs = pos.entrez$ENTREZID, Predicted_non_3_UTRs = neg.entrez$ENTREZID)

ck <- compareCluster(geneCluster = list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  universe = uni.entrez$ENTREZID,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 1,
  readable = TRUE)

ccRes <- ck@compareClusterResult
n = nrow(ccRes)
write.table(ccRes, str_c(outpath, "/compareClusterRes.txt"), row.names=FALSE, quote=FALSE, sep="\t")

png(str_c(outpath, "/compareClusterRes.png"), bg="transparent", units="in",width = 6.25, height= 6.00 , res=600)
dotplot(ck,
  showCategory = n,
  font.size = 10
) +
theme(axis.text.x = element_text(angle=45, hjust=1))
dev.off()
