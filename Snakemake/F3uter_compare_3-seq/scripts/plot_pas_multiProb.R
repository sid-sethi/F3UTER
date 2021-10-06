args <- commandArgs(TRUE)
suppressPackageStartupMessages({
  library(tidyverse)
  library(plotrix)
})


multiProb_res <- args[1] # "_f3uter_compare_pas_multiProb.txt"
respath <- args[2]
tissue <- args[3]


x <- read.table(multiProb_res, header=TRUE, sep="\t")

ppv <- x %>% dplyr::filter(metric %in% "PPV")

ymax = max(ppv$count) + 50

png(str_c(respath, "/", tissue, "_ppv_multiProb.png"), bg = "transparent", width = 5.25, height = 4.75, res = 600, units = "in")

par(las=0)
twoord.plot(ppv$cutoff, ppv$count,
            ppv$cutoff, ppv$value,
            xlab="3'UTR prediction probability",
            ylab="Number of 3'UTR predictions",
            rylab="3'UTR PPV (% value)",
            lylim=range(0, ymax),
            rylim=range(0,100),
            type=c("bar","b"),halfwidth=0.03,
            xtickpos = seq(0.1, 0.9, 0.1),
            xticklab = seq(0.1, 0.9, 0.1),
            lcol= 2,rcol= 4, lwd = 3,
            rpch = 17,
            main = str_c("3'UTR Positive Predictive Value - ", tissue)
            )

dev.off()



FOR <- x %>% dplyr::filter(metric %in% "FOR")

ymax <- max(FOR$count) + 50

png(str_c(respath, "/", tissue, "_FOR_multiProb.png"), bg = "transparent", width = 5.25, height = 4.75, res = 600, units = "in")

par(las=0)
twoord.plot(FOR$cutoff, FOR$count,
            FOR$cutoff, FOR$value,
            xlab="3'UTR prediction probability",
            ylab="Number of 3'UTR predictions",
            rylab="3'UTR FOR (% value)",
            lylim=range(0, ymax),
            rylim=range(0,100),
            type=c("bar","b"),halfwidth=0.03,
            xtickpos = seq(0.1, 0.9, 0.1),
            xticklab = seq(0.1, 0.9, 0.1),
            lcol= 2,rcol= "#762a83", lwd = 3,
            rpch = 19,
            main = str_c("3'UTR False Omission Rate - ", tissue)
            )

dev.off()
