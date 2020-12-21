args<-commandArgs(TRUE) # 1st argument: ML output directory

library(rlist)
library(matrixStats)
library(e1071)
library(tidyverse)

outPath = args[1]

system(paste("mkdir -m a=rwx ",outPath, "/Summary", sep=""))


overall <- list.files(path = outPath, pattern = "\\.overall.txt", recursive = TRUE, full.names = TRUE) %>%
  map(read.table, header=F, row.names = 1) %>%
  list.cbind() %>%
  rownames_to_column("stats") %>%
  mutate(mean = rowMeans(.[,-1]), sd = rowSds(as.matrix(.[,-1])))
file1 = paste(outPath, "/Summary/overall_stats.txt",sep="")
write.table(overall, file=file1, sep="\t",quote=F,row.names=F, col.names=T)

fp <- list.files(path = outPath, pattern = "\\.fp.txt", recursive = TRUE, full.names = TRUE) %>%
  map(read.table, header=F, row.names = 1) %>%
  list.cbind() %>%
  rownames_to_column("reference") %>%
  mutate(mean = rowMeans(.[,-1]), sd = rowSds(as.matrix(.[,-1])))
file2 = paste(outPath, "/Summary/fp.txt",sep="")
write.table(fp, file=file2, sep="\t",quote=F,row.names=F, col.names=T)

fn <- list.files(path = outPath, pattern = "\\.fn.txt", recursive = TRUE, full.names = TRUE) %>%
  map(read.table, header=F, row.names = 1) %>%
  list.cbind() %>%
  rownames_to_column("reference") %>%
  mutate(mean = rowMeans(.[,-1]), sd = rowSds(as.matrix(.[,-1])))
file3 = paste(outPath, "/Summary/fn.txt",sep="")
write.table(fn, file=file3, sep="\t",quote=F,row.names=F, col.names=T)


utr <- list.files(path = outPath, pattern = "\\.utr.txt", recursive = TRUE, full.names = TRUE) %>%
  map(read.table, header=F, row.names = 1, sep="\t") %>%
  list.cbind() %>%
  rownames_to_column("stats") %>%
  mutate(mean = rowMeans(.[,-1]), sd = rowSds(as.matrix(.[,-1])))
file4 = paste(outPath, "/Summary/utr.txt",sep="")
write.table(utr, file=file4, sep="\t",quote=F,row.names=F, col.names=T)


five <- list.files(path = outPath, pattern = "\\.five.txt", recursive = TRUE, full.names = TRUE) %>%
  map(read.table, header=F, row.names = 1, sep="\t") %>%
  list.cbind() %>%
  rownames_to_column("stats") %>%
  mutate(mean = rowMeans(.[,-1]), sd = rowSds(as.matrix(.[,-1])))
file5 = paste(outPath, "/Summary/five.txt",sep="")
write.table(five, file=file5, sep="\t",quote=F,row.names=F, col.names=T)


lncrna <- list.files(path = outPath, pattern = "\\.lncrna.txt", recursive = TRUE, full.names = TRUE) %>%
  map(read.table, header=F, row.names = 1, sep="\t") %>%
  list.cbind() %>%
  rownames_to_column("stats") %>%
  mutate(mean = rowMeans(.[,-1]), sd = rowSds(as.matrix(.[,-1])))
file6 = paste(outPath, "/Summary/lncrna.txt",sep="")
write.table(lncrna, file=file6, sep="\t",quote=F,row.names=F, col.names=T)


ncrna <- list.files(path = outPath, pattern = "\\.ncrna.txt", recursive = TRUE, full.names = TRUE) %>%
  map(read.table, header=F, row.names = 1, sep="\t") %>%
  list.cbind() %>%
  rownames_to_column("stats") %>%
  mutate(mean = rowMeans(.[,-1]), sd = rowSds(as.matrix(.[,-1])))
file7 = paste(outPath, "/Summary/ncrna.txt",sep="")
write.table(ncrna, file=file7, sep="\t",quote=F,row.names=F, col.names=T)


pseudo <- list.files(path = outPath, pattern = "\\.pseudo.txt", recursive = TRUE, full.names = TRUE) %>%
  map(read.table, header=F, row.names = 1, sep="\t") %>%
  list.cbind() %>%
  rownames_to_column("stats") %>%
  mutate(mean = rowMeans(.[,-1]), sd = rowSds(as.matrix(.[,-1])))
file8 = paste(outPath, "/Summary/pseudo.txt",sep="")
write.table(pseudo, file=file8, sep="\t",quote=F,row.names=F, col.names=T)


ice <- list.files(path = outPath, pattern = "\\.ice.txt", recursive = TRUE, full.names = TRUE) %>%
  map(read.table, header=F, row.names = 1, sep="\t") %>%
  list.cbind() %>%
  rownames_to_column("stats") %>%
  mutate(mean = rowMeans(.[,-1]), sd = rowSds(as.matrix(.[,-1])))
file9 = paste(outPath, "/Summary/ice.txt",sep="")
write.table(ice, file=file9, sep="\t",quote=F,row.names=F, col.names=T)
