args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
})

filea <- args[1]
fileb <- args[2]
universe <- args[3]
outpath <- args[4]
prefix <- args[5]


background_list <- read.table(universe, header=TRUE, sep="\t") %>%
  plyr::rename(c("Ensembl.gene.ID" = "gene_id"))
  #plyr::rename(c("associated_gene" = "gene_id"))

list_a <- read.table(filea, header=TRUE, sep="\t") %>%
  plyr::rename(c("associated_gene" = "gene_id"))


# GE panel list
list_b <- read.table(fileb, header=TRUE, sep="\t") %>%
  dplyr::select(Ensembl_id) %>%
  plyr::rename(c("Ensembl_id" = "gene_id")) %>%
  dplyr::inner_join(background_list, by = "gene_id")


overlap <- dplyr::semi_join(list_a, list_b, by = "gene_id")

write.table(overlap, file = str_c(outpath, "/", prefix, ".overlap.txt"), quote=FALSE, row.names=FALSE, sep="\t")


#########################################################
# list A
A = list_a$gene_id %>% unique() %>% length()
# list B
B = list_b$gene_id %>% unique() %>% length()
# overlap
o = overlap$gene_id %>% unique() %>% length()
# total genes in background
total = background_list$gene_id %>% unique() %>% length()
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
write.table(fisherRes, file=str_c(outpath, "/", prefix, "_fisherRes.txt"), quote=FALSE, col.names=FALSE, sep="\t")

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

write.table(hyperRes, file=str_c(outpath, "/", prefix, "_hyperRes.txt"), quote=FALSE, col.names=FALSE, row.names=FALSE)

### drawing the hypergeometric distribution

x.range <- 0:min(A,B)
# Compute the distribution of density P(X=x)
dens <- dhyper(x=x.range, m=B, n= total - B, k=A)


## Compute the distribution of P-values (for all possible values of x).
p.values <- phyper(q=x.range -1, m=B, n=total - B, k=A, lower.tail=FALSE)

## Plot the P-value distribution on linear scales
png(str_c(outpath, "/", prefix,"_hyperRes.png"), bg="transparent",units="in",width = 3.75, height= 4.25 ,res=600)

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
