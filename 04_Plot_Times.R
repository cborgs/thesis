library(PPRL)
library(fastdigest)
library(pryr)
library(ggplot2)
library(data.table)
library(tidyverse)

# Time results in days for 5-20 mio records (Mail Adrian Brown 28/01/16)
dat <- data.frame(size = c(5, 6.834999, 10, 15, 20), time = c(0.9326388889, 1.6979166667, 4.0888888889, 19.0833333333, 31.7104166667))

# Plot time in hours
p1 <- ggplot(dat, aes(x=size, y=time/60, color="white")) + geom_line(size=1, color="black") + geom_point(size=8, color="black")  
p1 <- p1  + xlab("Filesize (Mio)") + ylab("Time in Hours\n") + guides(colour=FALSE) 
p1 + theme_bw()  +  theme(panel.margin = unit(0.7, "lines"), text = element_text(size=36),  legend.justification=c(0.9,0.05), legend.key = element_blank(), legend.key.size = unit(3,"char"))
#p1 + theme_bw() + scale_colour_solarized("blue") 
#  theme_solarized(light = TRUE)  + scale_colour_solarized("blue") 
ggsave("AUS_Times_Q2016.pdf", height=12, width=14, units="in", device=cairo_pdf)


### Times 1 mio vs 1 mio

# Data: Curtin university evaluation of MBT with different thresholds
mbtTimes <- tribble(
  ~threshold, ~time,
  1,	0.161666666666664,
  0.95,	0.604500000000007,
  0.9,	1.90033333333333,
  0.85,	4.35499999999999,
  0.8,	20.8125
)

# Plot
p1 <- ggplot(data = mbtTimes, aes(x = threshold, y = (time)))
p1 <- p1 + geom_point(size=8) #  points
p1 <- p1 + geom_line(size=1) #  lines
p1 <- p1 + xlab("Tanimoto threshold") + ylab("Time in minutes") # labels
p1 <- p1 + scale_y_continuous(breaks = seq(0,22,2)) # set year breaks
p1 <- p1 + scale_x_continuous(breaks = seq(0.8,1.0,0.05)) # set year breaks

p1   + theme_bw() +  theme(axis.text.x = element_text(angle = 0, vjust = 0.1), 
                           panel.spacing = unit(0.85, "lines"), text = element_text(size=35),
                           legend.justification=c(0.9,0.05), legend.key = element_blank(), 
                           legend.key.size = unit(5,"char"), 
                           strip.background = element_rect(colour='white',fill = 'white')) # look
ggsave("./Time_Threshold2.pdf", 
       height=12, width=12, units="in", device=cairo_pdf()) # Fig 4.4


### Test: MBT variant and l  and threshold vs Time and RAM



# Variants to test
Variants <- c("mtan",
              "utan",
              "mtanp",
              "utanp")
# Set number of test cases and BF lengths and thresholds
N <- c(50000, 100000, 200000, 300000, 400000, 500000)
L <- c(250,500,1000)
tanimotoThresholds <- c(0.8,0.85,0.9,0.95)
# Base parameters
q <- 2
baseK <- 20
basel <- 1000
cores <- 6
leaflimit <- 3
set.seed(69)

# Counter for progress
passes <- length(N) * length(L) * length(Variants) * length(tanimotoThresholds)
i <- 0


# Init result
res <- matrix(NA, nrow=passes, ncol=6)


# Loop over sizes
for (n in N){
  
  
  dat <- sapply(rep(NA, n), function(x) fastdigest(sample(LETTERS, 512, replace = TRUE)))
  
  
  # 20% errors
  datB <- dat
  datB[(n - round(0.2 * n)):n] <- gsub("2","y", datB[(n - round(0.2 * n)):n])
  #datB[(n - round(0.2 * n)):n] <- gsub("1","x", datB[(n - round(0.2 * n)):n])
  
  for (l in L){
    
    # Make equivalent k dependent on length
    k <- round((baseK/basel) * l)
    
    
    # Gen BF
    testDataBF <- CreateBF(ID = as.character(1:n), dat,
                           k = k, padding = 0, q = 2, l = l, password = "(H]$6Uh04q")
    
    
    
    testDataBFB <- CreateBF(ID = as.character(1:n), datB,
                            k = k, padding = 0, q = 2, l = l, password = "(H]$6Uh04q")
    
    
    
    # Write Tree data
    write.table(testDataBF, "mbtSearch_2.0.1/TreeA.csv", sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
    write.table(testDataBFB, "mbtSearch_2.0.1/TreeB.csv", sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
    
    for (variant in Variants){
      
      ## Loop over thresholds
      for (Tani in tanimotoThresholds){
        
        # Start timer
        starttime <- proc.time()
        
        #memuse <- memory.size(TRUE)
        
        if (variant == "mtanp" | variant == "utanp"){
          # Kill the "p"
          variant <- substr(variant,1,4)
          # Run MBT
          system(paste("cd mbtSearch_2.0.1/ ;./mbtSearch -i TreeA.csv -m ", Tani, " -o MBT_Outfile3.csv -t ",cores," -a ",variant," -p -l ",leaflimit," TreeB.csv", sep=""))
          # Put "p" back in
          variant <- paste0(variant,"p")
        } else {
          # Just run it
          system(paste("cd mbtSearch_2.0.1/ ;./mbtSearch -i TreeA.csv -m ", Tani, " -o MBT_Outfile3.csv -t ",cores," -a ",variant," -l ",leaflimit," TreeB.csv", sep=""))
        }
        
        # Used memory
        memuse <- as.numeric(mem_used())/1000/1000
        
        # Time taken in minutes
        time <- proc.time() - starttime
        time <- as.numeric(time[3])/60
        
        # Read result file
        #result <- fread("mbtSearch_2.0.1/MBT_Outfile3.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
        
        i <- i + 1
        
        cat("\nPass ", i, "of ",passes, "finished (",round(i/passes, digits = 2),"%, time/RAM:",round(time),"/", round(memuse)," min/MB)", format(Sys.time()," (%a, %d-%m-%Y, %H:%M:%S)"))
        
        res[i,] <- c(variant, n, l, Tani, time, memuse)
      }    
    }
  }
}


# Make df
res <- as.data.frame(res, stringsAsFactors = FALSE)
# Names
names(res) <- c("Variant", "n", "l", "Threshold", "time", "memuse")
# Make numeric
res[,2:6] <- sapply(res[,2:6],as.numeric)

# Write outfile
write.table(res, "Timing_l_Tani.csv", sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)

#### Plotting

res <- fread("Timing_l_Tani.csv")

p5 <- ggplot(res, aes(x=as.numeric(n), y=time, color=Variant, group=Variant, shape=Variant)) + 
  geom_point(size=5) + geom_line(size=1) + 
  scale_x_continuous(breaks = seq(0,500000,50000)) +
  scale_shape_manual(name = "Tree type", labels= c("Multibit tree w/o Symdex","Union bit tree w/o Symdex","Multibit tree + Symdex","Union bit tree + Symdex"), values=c(15:18)) +
  scale_colour_manual(name = "Tree type", labels= c("Multibit tree w/o Symdex","Union bit tree w/o Symdex","Multibit tree + Symdex","Union bit tree + Symdex"), values=c("#004c93","#56B1F7", "#000000","#c43f4b")) +
  xlab("File Size") + ylab("Time in Minutes") +
  theme_bw() + facet_grid(l ~ Threshold, labeller = label_both) +
  theme(panel.spacing = unit(1, "lines"), text = element_text(size=27), axis.text.x = element_text(size = 21, angle = 90, vjust = 0.5),  axis.title=element_text(size=30),  legend.justification=c(0.9,0.05), legend.key = element_blank(), legend.key.size = unit(5,"char"), strip.background = element_rect(colour='white',fill = 'white'))
show(p5)
ggsave("Timing_l_Tani.pdf",height=13, width=21, units="in", device=cairo_pdf)




# Time: l vs RAM and l vs Time
subs <- subset(res, Threshold == 0.8 & Variant =="mtan")
subs$l <- factor(subs$l, levels = c("1000","500","250"))

p5 <- ggplot(subs, aes(x = as.numeric(n), y=memuse, color=l, group=l)) + 
  geom_point(size=5) + geom_line(size=1) + 
  scale_x_continuous(breaks = seq(0,500000,50000)) +
  scale_colour_manual(name = "Length l", values=c("#004c93","#56B1F7", "#000000","#c43f4b")) +
  ylab("Memory used in MB") + xlab("File Size") +
  theme_bw() +
  theme(panel.spacing = unit(1, "lines"), text = element_text(size=27), axis.text.x = element_text(size = 21, angle = 90, vjust = 0.5),  axis.title=element_text(size=30),  legend.justification=c(0.9,0.05), legend.key = element_blank(), legend.key.size = unit(5,"char"), strip.background = element_rect(colour='white',fill = 'white'))
show(p5)
ggsave("RAM_l500.pdf",height=13, width=21, units="in", device=cairo_pdf) # Fig 4.1


p5 <- ggplot(subs, aes(x = as.numeric(n), y=time, color=l, group=l)) + 
  geom_point(size=5) + geom_line(size=1) + 
  scale_x_continuous(breaks = seq(0,500000,50000)) +
  scale_colour_manual(name = "Length l", values=c("#004c93","#56B1F7", "#000000","#c43f4b")) +
  xlab("File Size") + ylab("Time in Minutes") +
  theme_bw() +
  theme(panel.spacing = unit(1, "lines"), text = element_text(size=27), axis.text.x = element_text(size = 21, angle = 90, vjust = 0.5),  axis.title=element_text(size=30),  legend.justification=c(0.9,0.05), legend.key = element_blank(), legend.key.size = unit(5,"char"), strip.background = element_rect(colour='white',fill = 'white'))
show(p5)
ggsave("Time_l500.pdf",height=13, width=21, units="in", device=cairo_pdf) # Fig 4.2




# Threshold vs Time vs RAM

subs <- subset(res, l == 500 & Variant =="mtan")
subs$Threshold <- factor(subs$Threshold, levels = c("0.8","0.85","0.9","0.95"), labels=c("0.80","0.85","0.90","0.95"))
  
p5 <- ggplot(subs, aes(x=as.numeric(n), y=time, color=Threshold, group=Threshold)) + 
  geom_point(size=5) + geom_line(size=1) + 
  scale_x_continuous(breaks = seq(0,500000,50000)) +
  scale_colour_manual(name = "Tanimoto threshold", values=c("#004c93","#56B1F7", "#000000","#c43f4b")) +
  xlab("File Size") + ylab("Time in Minutes") +
  theme_bw() +
  theme(panel.spacing = unit(1, "lines"), text = element_text(size=27), axis.text.x = element_text(size = 21, angle = 90, vjust = 0.5),  axis.title=element_text(size=30),  legend.justification=c(0.9,0.05), legend.key = element_blank(), legend.key.size = unit(5,"char"), strip.background = element_rect(colour='white',fill = 'white'))
show(p5)
ggsave("Time_Threshold.pdf",height=13, width=21, units="in", device=cairo_pdf) # Fig 4.5



# Tree type vs runtime
subs <- subset(res, l == 1000 & Threshold == 0.9)

p5 <- ggplot(subs, aes(x=as.numeric(n), y=time, color=Variant, group=Variant, shape=Variant)) + 
  geom_point(size=5) + geom_line(size=1) + 
  scale_x_continuous(breaks = seq(0,500000,50000)) +
  scale_shape_manual(name = "Tree type", labels= c("Multibit tree w/o Symdex","Union bit tree w/o Symdex","Multibit tree + Symdex","Union bit tree + Symdex"), values=c(15:18)) +
  scale_colour_manual(name = "Tree type", labels= c("Multibit tree w/o Symdex","Union bit tree w/o Symdex","Multibit tree + Symdex","Union bit tree + Symdex"), values=c("#004c93","#56B1F7", "#000000","#c43f4b")) +
  xlab("File Size") + ylab("Time in Minutes") +
  theme_bw() + facet_grid(l ~ Threshold, labeller = label_both) +
  theme(panel.spacing = unit(1, "lines"), text = element_text(size=27), axis.text.x = element_text(size = 21, angle = 90, vjust = 0.5),  axis.title=element_text(size=30),  legend.justification=c(0.9,0.05), legend.key = element_blank(), legend.key.size = unit(5,"char"), strip.background = element_rect(colour='white',fill = 'white'))
show(p5)
ggsave("Times_MBT_type.pdf",height=13, width=21, units="in", device=cairo_pdf)

# Times vs Tree type for l = 500/l=1000
subs <- subset(res, l == 500 | l == 1000)

p5 <- ggplot(subs, aes(x=as.numeric(n), y=time, color=Variant, group=Variant, shape=Variant)) + 
  geom_point(size=5) + geom_line(size=1) +
  scale_x_continuous(breaks = seq(0,500000,50000)) +
  scale_shape_manual(name = "Tree type", labels= c("Multibit tree w/o Symdex","Union bit tree w/o Symdex","Multibit tree + Symdex","Union bit tree + Symdex"), values=c(15:18)) +
  scale_colour_manual(name = "Tree type", labels= c("Multibit tree w/o Symdex","Union bit tree w/o Symdex","Multibit tree + Symdex","Union bit tree + Symdex"), values=c("#004c93","#56B1F7", "#000000","#c43f4b")) +
  xlab("File Size") + ylab("Time in Minutes") +
  theme_bw() + facet_grid(l ~ Threshold, labeller = label_both) +
  theme(panel.spacing = unit(1, "lines"), text = element_text(size=27), axis.text.x = element_text(size = 21, angle = 90, vjust = 0.5),  axis.title=element_text(size=30),  legend.justification=c(0.9,0.05), legend.key = element_blank(), legend.key.size = unit(5,"char"), strip.background = element_rect(colour='white',fill = 'white'))
show(p5)
ggsave("Times_MBT_type.pdf",height=13, width=21, units="in", device=cairo_pdf) # Fig 3.10


