library(tidyverse)
library(rtracklayer)
library(scales)
library(BSgenome)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38.masked)

options(stringsAsFactors=F)


# building the position frequency matrix for PAS
pfm <- PFMatrix(ID="Motif_all", name="PolyA_signals",
                bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                profileMatrix=matrix(c(9L,  9L, 0L,  11L,  9L,  12L,
                                       1L, 1L,  0L, 0L,  1L,  0L,
                                       1L,  1L,  1L,  1L, 1L,  0L,
                                       1L,  1L,  11L,  0L,  1L, 0L),
                                     byrow=TRUE, nrow=4,
                                     dimnames=list(c("A", "C", "G", "T"))
                )
)

# conversion to position weight matrix
pwm <- TFBSTools::toPWM(pfm, type="log2probratio", pseudocounts=0.8,
                        bg=c(A=0.25, C=0.25, G=0.25, T=0.25))


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


# searching sequences with 95% match to PWM (time intensive commands)
utr.sites <- searchSeq(pwm, utr.seq, min.score="95%", strand="*")
save(utr.sites, file = "/PASignals/TfbsSTools.RData")

#loading pre-computed results of the above commands
#load("/PASignals/TfbsSTools.RData")


utr.relscore = relScore(utr.sites)
utr.hits.len = sapply(utr.relscore, length) %>% as.data.frame()
utr.hits.len$id = rownames(utr.hits.len)
utr.export = utr.hits.len %>% mutate(polyA_signal = as.factor(ifelse(. > 0, "1", .)), class = "3-UTR") %>%
  select(-c(.))


write.table(utr.export, file = "/Features/polyA_signal.txt", row.names = FALSE, quote = FALSE, sep="\t")
