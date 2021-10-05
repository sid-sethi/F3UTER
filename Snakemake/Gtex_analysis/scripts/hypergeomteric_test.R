args<-commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
})

file_a <- args[1]
file_b <- args[2]
universe <- args[3]


list_a <- read.table(file_a, header=TRUE, sep="\t")

bg_list <- read.table(universe, header=TRUE, sep="\t")

# GE panel list
list_b <- read.table("/tmp_mnt/filer1/bioinformatics/Sid/Data/GenomicsEngland/NeurodegenerativeDisorders_adultOnset_greenGenes.txt", header=TRUE, sep="\t") %>%
  #dplyr::select(Ensembl_id) %>%
  dplyr::inner_join(bg_list, by = c("Ensembl_id" = "gene_id"))





## OMIM data
# list_b <- read_delim("/tmp_mnt/filer1/bioinformatics/Sid/OMIM_data/OMIM_data_from_David.csv", delim = ",") %>%
#   dplyr::filter(!is.na(neurologic)) %>%
#   dplyr::select(ensembl_gene_id) %>%
#   dplyr::inner_join(background_list, by = c("ensembl_gene_id" = "Gene.stable.ID"))





## TARDBP iclip targets
# list_b <- read_csv("/tmp_mnt/filer1/bioinformatics/Sid/Features_training_all/Ml_evaluate_23_Jan_20/MakePredictions/Gtex/Intergenic_3/PredictionApp/RBP_scan_3/TARDBP_target_overlap/POSTAR_tardbp_iclip.csv") %>%
#   dplyr::filter(`Target gene type` %in% "protein_coding") %>%
#   dplyr::select(`Target gene ID`) %>%
#   dplyr::inner_join(background_list, by = c(`Target gene ID` = "Gene.stable.ID"))



overlap <- dplyr::semi_join(list_a, list_b, by = c("associated_gene" = "Target gene ID"))

write.table(overlap, file = str_c(out.path, "/", out, ".overlap.txt"), quote=FALSE, row.names=FALSE, sep="\t")


#########################################################
# list A
A = list_a$associated_gene %>% unique() %>% length()
# list B
B = list_b$`Target gene ID` %>% unique() %>% length()
# overlap
o = overlap$associated_gene %>% unique() %>% length()
# total genes in background
total = background_list$Gene.stable.ID %>% unique() %>% length()
########################################################


#######################################
############ FISHER TEST ##############
#######################################

inA_inB <- o
inA_notB <- A - o
inB_notA <- B - o
notA_notB <- total - (inA_notB + B)

cont.table <- matrix(c(notA_notB,inB_notA,inA_notB,inA_inB), nrow = 2)
rownames(cont.table) <- c("notB", "inB")
colnames(cont.table) <- c("notA", "inA")


x <- fisher.test(cont.table) %>% unlist() %>% as.data.frame()
y <- fisher.test(cont.table, alternative = "g") %>% unlist() %>% as.data.frame()
z <- fisher.test(cont.table, alternative = "l") %>% unlist() %>% as.data.frame()

fisherRes <- cbind(x,y,z)
write.table(fisherRes, file=str_c(out.path, "/", out, "_fisherRes.txt"), quote=FALSE, col.names=FALSE, sep="\t")

#######################################
######### HYPERGEOMETRIC TEST #########
#######################################


p.value <- phyper(q=o -1, m=B, n= total - B, k=A, lower.tail=FALSE)

## Random expectation
marked.proportion <- B / total
exp.x <- A * marked.proportion
## Fold enrichment, as computed by David
fold.enrichment <-  (o / A ) / (B / total)

res1 <- paste("hypergeometric p-value = ", p.value, sep="")
res2 <- paste("fold enrichment = ", fold.enrichment, sep="")
hyperRes <- paste(res1, res2, sep="\n")

write.table(hyperRes, file=str_c(out.path, "/", out, "_hyperRes.txt"), quote=FALSE, col.names=FALSE, row.names=FALSE)

### drawing the hypergeometric distribution

x.range <- 0:min(A,B)
# Compute the distribution of density P(X=x)
dens <- dhyper(x=x.range, m=B, n= total - B, k=A)


## Compute the distribution of P-values (for all possible values of x).
p.values <- phyper(q=x.range -1, m=B, n=total - B, k=A, lower.tail=FALSE)

## Plot the P-value distribution on linear scales
png(str_c(out.path, "/", out,"_hyperRes.png"),bg="transparent",units="in",width = 3.75, height= 4.25 ,res=600)

plot(x.range, p.values, type="l", lwd=2, col="violet", main="Hypergeometric P-value", xlab="x = number of overlapping elements in the selection", ylab="P-value = P(X>=x)", ylim=c(0, 1), panel.first=grid(), cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.6)
## We can plot the density below the P-value
lines(x.range, dens, type="h", col="blue", lwd=2)
legend("topright", legend=c("P-value", "density"), lwd=2, col=c("violet", "blue"), bg="white", bty="o", cex=0.6)

## Arrow indicating observed and expected value
arrows(exp.x, max(dens)*1.15, exp.x, max(dens)*1.05, lwd=2, col='darkgreen', angle=30, length=0.05, code=2)
text(A*.8, 0.7, labels=paste("exp=", round(exp.x, digits=2)), col="darkgreen", font=2, cex=0.6)
arrows(o, max(dens)*1.35, o, max(dens)*1.1, lwd=2, col='red', angle=30, length=0.05, code=2)
text(A*.8, 0.6, labels=paste("x=", o, "; p-val=", signif(digits=1, p.value), sep=""), col="red", font=2, cex=0.6)

dev.off()
