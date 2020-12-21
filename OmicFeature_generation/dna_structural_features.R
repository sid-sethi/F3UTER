library(tidyverse)
library(scales)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38.masked)

options(stringsAsFactors=F)


# conversion table path
pathTable = "/DNA_structural_properties/Conversion_tables"

# loading the UTR training regions
# this code can be modified to be used on other training regions and on ERs
load("/GS.RData")

# converting to Granges object
utr.gr = makeGRangesFromDataFrame(utr, keep.extra.columns = TRUE)

# adding "chr" in front of seqnames
newStyle <- mapSeqlevels(seqlevels(utr.gr), "UCSC")
utr.gr <- renameSeqlevels(utr.gr, newStyle)


# Extract sequences (using the masked version)
utr.seq = BSgenome::getSeq(Hsapiens, utr.gr)
utr.seq@ranges@NAMES <- utr.gr$three_prime_utr_id


################################################
############ Conversion tables #################
################################################

aphilicity <- read.table(paste(pathTable, "/aphilicity.tsv", sep=""), header=FALSE)
colnames(aphilicity) = c("nuc", "value")
aphilicity$nuc = stringr::str_to_upper(aphilicity$nuc)

basestacking <- read.table(paste(pathTable, "/basestacking.tsv", sep=""), header=FALSE)
colnames(basestacking) = c("nuc", "value")
basestacking$nuc = stringr::str_to_upper(basestacking$nuc)

bdnatwist <- read.table(paste(pathTable, "/bdnatwist.tsv", sep=""), header=FALSE)
colnames(bdnatwist) = c("nuc", "value")
bdnatwist$nuc = stringr::str_to_upper(bdnatwist$nuc)

bendability <- read.table(paste(pathTable, "/bendability.tsv", sep=""), header=FALSE)
colnames(bendability) = c("nuc", "value")
bendability$nuc = stringr::str_to_upper(bendability$nuc)

bendingstiffness <- read.table(paste(pathTable, "/bendingstiffness.tsv", sep=""), header=FALSE)
colnames(bendingstiffness) = c("nuc", "value")
bendingstiffness$nuc = stringr::str_to_upper(bendingstiffness$nuc)

cpgislands <- read.table(paste(pathTable, "/cpgislands.tsv", sep=""), header=FALSE)
colnames(cpgislands) = c("nuc", "value")
cpgislands$nuc = stringr::str_to_upper(cpgislands$nuc)

cpnpgcpgislands <- read.table(paste(pathTable, "/cpnpgcpgislands.tsv", sep=""), header=FALSE)
colnames(cpnpgcpgislands) = c("nuc", "value")
cpnpgcpgislands$nuc = stringr::str_to_upper(cpnpgcpgislands$nuc)

cpnpgislands <- read.table(paste(pathTable, "/cpnpgislands.tsv", sep=""), header=FALSE)
colnames(cpnpgislands) = c("nuc", "value")
cpnpgislands$nuc = stringr::str_to_upper(cpnpgislands$nuc)

dnadenaturation <- read.table(paste(pathTable, "/dnadenaturation.tsv", sep=""), header=FALSE)
colnames(dnadenaturation) = c("nuc", "value")
dnadenaturation$nuc = stringr::str_to_upper(dnadenaturation$nuc)

duplexstabilitydisruptenergy <- read.table(paste(pathTable, "/duplexstabilitydisruptenergy.tsv", sep=""), header=FALSE)
colnames(duplexstabilitydisruptenergy) = c("nuc", "value")
duplexstabilitydisruptenergy$nuc = stringr::str_to_upper(duplexstabilitydisruptenergy$nuc)

duplexstabilityfreeenergy <- read.table(paste(pathTable, "/duplexstabilityfreeenergy.tsv", sep=""), header=FALSE)
colnames(duplexstabilityfreeenergy) = c("nuc", "value")
duplexstabilityfreeenergy$nuc = stringr::str_to_upper(duplexstabilityfreeenergy$nuc)

nucleosomeposition <- read.table(paste(pathTable, "/nucleosomeposition.tsv", sep=""), header=FALSE)
colnames(nucleosomeposition) = c("nuc", "value")
nucleosomeposition$nuc = stringr::str_to_upper(nucleosomeposition$nuc)

propellortwist <- read.table(paste(pathTable, "/propellortwist.tsv", sep=""), header=FALSE)
colnames(propellortwist) = c("nuc", "value")
propellortwist$nuc = stringr::str_to_upper(propellortwist$nuc)

proteindeformation <- read.table(paste(pathTable, "/proteindeformation.tsv", sep=""), header=FALSE)
colnames(proteindeformation) = c("nuc", "value")
proteindeformation$nuc = stringr::str_to_upper(proteindeformation$nuc)

proteindnatwist <- read.table(paste(pathTable, "/proteindnatwist.tsv", sep=""), header=FALSE)
colnames(proteindnatwist) = c("nuc", "value")
proteindnatwist$nuc = stringr::str_to_upper(proteindnatwist$nuc)

zdna <- read.table(paste(pathTable, "/zdna.tsv", sep=""), header=FALSE)
colnames(zdna) = c("nuc", "value")
zdna$nuc = stringr::str_to_upper(zdna$nuc)

################################################
################################################


###### main function ############
seq_to_score <- function(sequence, conversion_table){

  vector = c()
  for(i in 1:nrow(conversion_table)){
    nuc = conversion_table[i,]$nuc
    value = conversion_table[i,]$value
    regx = paste0("(?<=(", nuc, "))")
    count.list <- gregexpr(regx,sequence, perl=TRUE)
    if(count.list[[1]][1] > 0){
      times = sapply(count.list, length)
      vector = c(vector, rep(value, times))
    }
  }

  #return(median(vector))
  return(mean(vector))

}
##############################

# old code,  slow (but accurate!)
# seq_to_score2 <- function(sequence, conversion_table, step){ # step = di-nucleotide(2) or tri-nucleotide(3)
#   width = length(sequence)
#   end = width - (step-1)
#
#   res = data.frame()
#   for(s in 1:end){
#     window = stringr::str_sub(sequence, s, s + (step-1))
#     score = conversion_table %>% dplyr::filter(nuc %in% window)
#     res = rbind(res, score)
#   }
#
#   return(median(res$value))
#
# }


tab1 = sapply(utr.seq, function(x) seq_to_score(x, aphilicity)) %>%
  as.data.frame() %>%
  plyr::rename(c("." = "aphilicity")) %>%
  rownames_to_column("id")

tab2 = sapply(utr.seq, function(x) seq_to_score(x, basestacking)) %>%
  as.data.frame() %>%
  plyr::rename(c("." = "basestacking"))

tab3 = sapply(utr.seq, function(x) seq_to_score(x, bdnatwist)) %>%
  as.data.frame() %>%
  plyr::rename(c("." = "bdnatwist"))

tab4 = sapply(utr.seq, function(x) seq_to_score(x, bendability)) %>%
  as.data.frame() %>%
  plyr::rename(c("." = "bendability"))

tab5 = sapply(utr.seq, function(x) seq_to_score(x, bendingstiffness)) %>%
  as.data.frame() %>%
  plyr::rename(c("." = "bendingstiffness"))

tab6 = sapply(utr.seq, function(x) seq_to_score(x, cpgislands)) %>%
  as.data.frame() %>%
  plyr::rename(c("." = "cpgislands"))

tab7 = sapply(utr.seq, function(x) seq_to_score(x, cpnpgcpgislands)) %>%
  as.data.frame() %>%
  plyr::rename(c("." = "cpnpgcpgislands"))

tab8 = sapply(utr.seq, function(x) seq_to_score(x, cpnpgislands)) %>%
  as.data.frame() %>%
  plyr::rename(c("." = "cpnpgislands"))

tab9 = sapply(utr.seq, function(x) seq_to_score(x, dnadenaturation)) %>%
  as.data.frame() %>%
  plyr::rename(c("." = "dnadenaturation"))

tab10 = sapply(utr.seq, function(x) seq_to_score(x, duplexstabilitydisruptenergy)) %>%
  as.data.frame() %>%
  plyr::rename(c("." = "duplexstabilitydisruptenergy"))

tab11 = sapply(utr.seq, function(x) seq_to_score(x, duplexstabilityfreeenergy)) %>%
  as.data.frame() %>%
  plyr::rename(c("." = "duplexstabilityfreeenergy"))

tab12 = sapply(utr.seq, function(x) seq_to_score(x, nucleosomeposition)) %>%
  as.data.frame() %>%
  plyr::rename(c("." = "nucleosomeposition"))

tab13 = sapply(utr.seq, function(x) seq_to_score(x, propellortwist)) %>%
  as.data.frame() %>%
  plyr::rename(c("." = "propellortwist"))

tab14 = sapply(utr.seq, function(x) seq_to_score(x, proteindeformation)) %>%
  as.data.frame() %>%
  plyr::rename(c("." = "proteindeformation"))

tab15 = sapply(utr.seq, function(x) seq_to_score(x, proteindnatwist)) %>%
  as.data.frame() %>%
  plyr::rename(c("." = "proteindnatwist"))

tab16 = sapply(utr.seq, function(x) seq_to_score(x, zdna)) %>%
  as.data.frame() %>%
  plyr::rename(c("." = "zdna"))

utr.export = bind_cols(tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8, tab9, tab10, tab11, tab12, tab13, tab14, tab15, tab16) %>%
  mutate(class = "3-UTR")



write.table(utr.export, file = "/Features/structural_feat_mean.txt", row.names = FALSE, quote = FALSE, sep="\t")
