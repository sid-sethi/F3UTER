library(tidyverse)

# false-positives

fp <- read.table("fp.txt", header=TRUE, sep="\t")

fp.df2 <- fp %>% dplyr::select(reference, mean) %>%
  mutate(group = "Predicted 3'UTRs")

fp.df2$reference <- recode(fp.df2$reference, five_prime_UTR = "5-UTR", three_prime_UTR = "3-UTR")

png("fp_barplot.png",bg="transparent",units="in",width = 3.25, height= 3.75 ,res=600)
ggplot(data=fp.df2, aes(x=group, y=mean, fill = reference)) +
  geom_bar(stat="identity", color = "black", alpha=0.8) +
  scale_fill_manual(values = c("#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999", "#00AFBB")) +
  theme_bw() +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.40,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        #legend.position= c(0.85,0.90),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=11, face = "bold"),
        axis.text.y = element_text(size=10),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=11,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_y_continuous(name = "Actual\n(% composition)")+
  scale_x_discrete(labels = c("Predicted\n3'UTRs"))
dev.off()



# false-negatives

fn <- read.table("fn.txt", header=TRUE, sep="\t")

fn.df2 <- fn %>% dplyr::select(reference, mean) %>%
  mutate(group = "Actual 3'UTRs")

fn.df2$reference <- recode(fn.df2$reference, five_prime_UTR = "5-UTR", three_prime_UTR = "3-UTR")

png("fn_barplot.png",bg="transparent",units="in",width = 3.25, height= 3.75 ,res=600)
ggplot(data=fn.df2, aes(x=group, y=mean, fill = reference)) +
  geom_bar(stat="identity", color = "black", alpha=0.8) +
  scale_fill_manual(values = c("#404040", "#D1A700", "#CC00CC", "#009900", "#FF9999", "#00AFBB")) +
  theme_bw() +
  theme(plot.title = element_text(size=13, hjust=0.5),
        legend.key.size = unit(0.40,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        #legend.position= c(0.85,0.90),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=11, face = "bold"),
        axis.text.y = element_text(size=10),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=11,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_y_continuous(name = "Predicted\n(% composition)")+
  scale_x_discrete(labels = c("Actual\n3'UTRs"))
dev.off()
