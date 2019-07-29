library(tidyverse)
library(data.table)

#setwd("E:/Priv/Dropbox/Dissertation/Programme")

# Read files
DE <- fread("Daten/NN_Bigram_Frequencies_DE.csv")
PL <- fread("Daten/NN_Bigram_Frequencies_PO.csv")

# Get top 50 DE bigrams
names <- head(DE[order(DE$freq, decreasing=TRUE), ],50)

# Make factors for plotting, reorder by bigram count
names$bigram <- factor(names$bigram)
names$bigram <- reorder(names$bigram, names$freq, mean)

# Scientific notation disabled
options(scipen=10000)

# counts for DE
p1 <- ggplot(data=names, aes(x=bigram, y=freq)) + geom_bar(stat="identity", fill="#004c93") + coord_flip()  +  theme_bw() + 
  scale_y_continuous(expand = c(0,0), limits = c(0, (max(names$freq) + round(max(names$freq) * 0.05 ) ))) +  xlab("Bigram") + ylab("\nCount")
p1 <- p1 + ggtitle("German Nationals\n")
p1   + theme_bw() +  theme(axis.text.x = element_text(angle = 0, vjust = 0.1), axis.text.y = element_text(size=20, hjust=0.5),
                             panel.spacing = unit(0.85, "lines"), text = element_text(size=35),
                             legend.justification=c(0.9,0.05), legend.key = element_blank(),  plot.title = element_text(size=22),
                             legend.key.size = unit(5,"char"),  strip.background = element_rect(colour='white',fill = 'white')) # look
# Save plot
ggsave("results/Bigram_Counts_DE.pdf", 
       height=16, width=10, units="in", device=cairo_pdf()) # Fig 4.3a



# Get top 50 for PL
names <- head(PL[order(PL$freq, decreasing=TRUE), ],50)


# Make factors for plotting, reorder by bigram count
names$bigram <- factor(names$bigram)
names$bigram <- reorder(names$bigram, names$freq, mean)

# Scientific notation disabled
options(scipen=10000)

# counts for pL
p1 <- ggplot(data=names, aes(x=bigram, y=freq)) + geom_bar(stat="identity", fill="#004c93") + coord_flip()  +  theme_bw() + 
  scale_y_continuous(expand = c(0,0), limits = c(0, (max(names$freq) + round(max(names$freq) * 0.05 ) ))) +  xlab("Bigram") + ylab("\nCount")
p1 <- p1 + ggtitle("Polish Nationals\n")
p1   + theme_bw() +  theme(axis.text.x = element_text(angle = 0, vjust = 0.1), axis.text.y = element_text(size=20, hjust=0.5),
                           panel.spacing = unit(0.85, "lines"), text = element_text(size=35),
                           legend.justification=c(0.9,0.05), legend.key = element_blank(),  plot.title = element_text(size=22),
                           legend.key.size = unit(5,"char"),  strip.background = element_rect(colour='white',fill = 'white')) # look
#Save plot
ggsave("results/Bigram_Counts_PL.pdf", 
       height=16, width=10, units="in", device=cairo_pdf()) # Fig 4.3b

