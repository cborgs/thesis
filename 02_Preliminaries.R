
#setwd("E:/Priv/Dropbox/Dissertation/Programme")

library(PPRL)
library(fastdigest)
library(xtermStyle)
library(tidyverse)
library(reshape2)
library(multibitTree)
library(gmodels)


# Pre-tests for properties of Bloom Filters



# Set number of test cases and BF length
n <- 1000
l <- 256
q <- 2
set.seed(69)

# Sim data by hashing sampled letters
dat <- replicate(n ,fastdigest(sample(LETTERS, 10000, replace = TRUE)))





# 20% errors
datB <- dat
datB[(n - round(0.2 * n)):n] <- gsub("2","y", datB[(n - round(0.2 * n)):n])
#datB[(n - round(0.2 * n)):n] <- gsub("1","x", datB[(n - round(0.2 * n)):n])



# Init result matrix
results <- NULL

# Look over k
for (k in seq(from = 1, to = 30)){
  
  # Gen BF
  testDataBF <- CreateBF(ID = as.character(1:n), dat,
                         k = k, padding = 0, q = 2, l = l, password = "(H]$6Uh04q")
  
  testDataBFB <- CreateBF(ID = as.character(1:n), datB,
                         k = k, padding = 0, q = 2, l = l, password = "(H]$6Uh04q")
  
  
  write.table(testDataBF, "BFA.csv", sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)
  write.table(testDataBFB, "BFB.csv", sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)
  
  # Empirical fill and k/l
  fill <- mean(nchar(gsub("0","",testDataBF$CLKs))) / l
  quotient <- k / l
  elements <- k * (32/q)

  
  # Link
  multibitTree.load("BFB.csv",threads=7, leafLimit = 3)
  res <- multibitTree.searchFile("BFA.csv", 0.85)
  

  
  # Calc quality measures
  tp <- sum(as.character(res$query) == as.character(res$fingerprint))
  fp <- sum(as.character(res$query) != as.character(res$fingerprint))
  fn <- n - tp
  recall <-  (tp/(tp+fn))
  precision <-  (tp/(tp+fp))
  fscore <- (recall + precision) / 2
  
  # Print result
  cat("\n",style("\t    PPRL Test\t\t", bg = "black", fg = "white", font="bold"),
      "\n Positives:\t\t", n,"\n k:\t\t\t", k,"\n Fill:\t\t\t", fill,
      "\n True positives:\t", tp, 
      "\n False positives:\t", fp, "\n False Negatives:\t", fn,
      "\n Recall:\t\t", (tp/(tp+fn)),"\n Precision:\t\t", (tp/(tp+fp)),
      "\n",style("\t\t   \t\t", bg = "black", fg = "white", font="bold"))
  
  # Append results
  results <- rbind(results, data.frame(data="Synthetic",k, q, l, fill,elements, quotient, tp, fp, fn, precision, recall, fscore))
}

# Plotting
#plotdata <- results
write.table(results, "./results/fill_k_quality.csv", sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)


plotdata <- read.csv("./results/fill_k_quality.csv", head=TRUE, sep = "\t", stringsAsFactors = FALSE)

# No. of elements hashed. Fastdigest returns 32 bit hashes, divided into q-grams, hashes k times.


library(ggplot2)

p1 <- ggplot(data = plotdata, aes(y = round(fill * 100), x = k))
p1 <- p1 + geom_smooth(se=F, method="loess", color="#004c93") #  method="loess"
p1 <- p1 + geom_point(size=6) #  
p1 <- p1 + xlab("\nHash functions k") + ylab("% of bit positions set to one\n")
p1 <- p1 + scale_y_continuous(breaks = seq(0,100,10), limits=  c(1,100), expand = c(0.02, 0))
p1 <- p1 + scale_x_continuous(breaks = seq(0,30,5))
p1   + theme_bw() +  theme(panel.spacing = unit(0.85, "lines"), text = element_text(size=35),legend.justification=c(0.9,0.05), legend.key = element_blank(), legend.key.size = unit(5,"char"), strip.background = element_rect(colour='white',fill = 'white'))
ggsave("./results/k_vs_fill.pdf", height=11, width=15, units="in", device=cairo_pdf) # Fig 3.15

p2 <- ggplot(data = plotdata, aes(y = fscore, x = round(fill * 100)))
p2 <- p2 + geom_smooth(se=F, method="loess", color="#004c93", span=.25) #  method="loess"
p2 <- p2 + geom_point(size=6) #  methos="loess"
p2 <- p2 + geom_vline(xintercept=50, colour="#004c93") +
  geom_text(aes(x=50, label="\nTheoretical optimum", y=0.8), colour="#004c93", angle=90, size=6)
p2 <- p2 + scale_x_continuous(breaks = seq(0,100,10), limits=  c(1,100), expand = c(0.02, 0))
p2 <- p2 + xlab("\n% of bit positions set to one") + ylab("F-Score\n")
p2   + theme_bw() +  theme(panel.spacing = unit(0.85, "lines"), text = element_text(size=35),legend.justification=c(0.9,0.05), legend.key = element_blank(), legend.key.size = unit(5,"char"), strip.background = element_rect(colour='white',fill = 'white'))
ggsave("./results/fill_vs_f.pdf", height=11, width=15, units="in", device=cairo_pdf) # Fig  3.18



p2 <- ggplot(data = plotdata, aes(y = precision, x = round(fill * 100)))
p2 <- p2 + geom_smooth(se=F, method="loess", color="#004c93", span=.25) #  method="loess"
p2 <- p2 + geom_point(size=6) #  methos="loess"
p2 <- p2 + geom_vline(xintercept=50, colour="#004c93") +
  geom_text(aes(x=50, label="\nTheoretical optimum", y=0.8), colour="#004c93", angle=90, size=6)
p2 <- p2 + scale_x_continuous(breaks = seq(0,100,10), limits=  c(1,100), expand = c(0.02, 0))
p2 <- p2 + xlab("\n% of bit positions set to one") + ylab("Precision\n")
p2   + theme_bw() +  theme(panel.spacing = unit(0.85, "lines"), text = element_text(size=35),legend.justification=c(0.9,0.05), legend.key = element_blank(), legend.key.size = unit(5,"char"), strip.background = element_rect(colour='white',fill = 'white'))
ggsave("./results/fill_vs_prec.pdf", height=11, width=15, units="in", device=cairo_pdf) # Fig 3.16




p2 <- ggplot(data = plotdata, aes(y = recall, x = round(fill * 100)))
p2 <- p2 + geom_smooth(se=F, method="loess", color="#004c93", span=.25) #  method="loess"
p2 <- p2 + geom_point(size=6) #  methos="loess"
#p2 <- p2 + geom_vline(xintercept=50, colour="#004c93") +
#  geom_text(aes(x=50, label="\nTheoretical optimum", y=0.8), colour="#004c93", angle=90, size=6)
p2 <- p2 + scale_x_continuous(breaks = seq(0,100,10), limits=  c(1,100), expand = c(0.02, 0))
p2 <- p2 + xlab("\n% of bit positions set to one") + ylab("Recall\n")
p2   + theme_bw() +  theme(panel.spacing = unit(0.85, "lines"), text = element_text(size=35),legend.justification=c(0.9,0.05), legend.key = element_blank(), legend.key.size = unit(5,"char"), strip.background = element_rect(colour='white',fill = 'white'))
ggsave("./results/fill_vs_rec.pdf", height=11, width=15, units="in", device=cairo_pdf) # Fig 3.17

p3 <- ggplot(data = plotdata, aes(y = fscore, x = k))
p3 <- p3 + geom_smooth(se=F, method="loess", color="#004c93", span=.25) #  method="loess"
p3 <- p3 + geom_point(size=6)  
p3 <- p3 + geom_point(aes(y = 0.976, x = 11), color="#004c93", size=6)  
p3 <- p3 + geom_text(aes(y = 1.02, x = 11, label="Optimal k by Smith (2012)"), color="#004c93", size=6)  
p3 <- p3 + xlab("\nHash functions k") + ylab("F-Score\n")
p3   + theme_bw() +  theme(panel.spacing = unit(0.85, "lines"), text = element_text(size=35),legend.justification=c(0.9,0.05), legend.key = element_blank(), legend.key.size = unit(5,"char"), strip.background = element_rect(colour='white',fill = 'white'))
ggsave("./results/k_vs_f.pdf", height=11, width=15, units="in", device=cairo_pdf) # Fig 3.19


p4 <- ggplot(data = plotdata, aes(y = fscore, x = elements))
p4 <- p4 + geom_smooth(se=F, method="loess", color="#004c93", span=.25) #  method="loess"
p4 <- p4 + geom_point(size=6)  
p4 <- p4 + xlab("\nNo. of Elements hashed") + ylab("F-Score\n")
p4   + theme_bw() +  theme(panel.spacing = unit(0.85, "lines"), text = element_text(size=35),legend.justification=c(0.9,0.05), legend.key = element_blank(), legend.key.size = unit(5,"char"), strip.background = element_rect(colour='white',fill = 'white'))
ggsave("./results/elements_vs_f.pdf", height=11, width=15, units="in", device=cairo_pdf)



# AMount of hash functions for 50% fill acc to smith 2012:
(log(1-0.5) / log(1-(1/l))) / 16
# [1] 11.06868



# Init result matrix
results <- NULL
l = 500

set.seed(1)

# Bootstrap with pws
for (k in c(5,10, 15)){
  for (i in seq(from = 1, to = 30)){
  

  PWs <- replicate(length(1), fastdigest(sample(LETTERS, 10000, replace = TRUE)))
  # Gen BF
  testDataBF <- CreateBF(ID = as.character(1:n), dat,
                         k = k, padding = 0, q = 2, l = l, password = PWs)
  
  testDataBFB <- CreateBF(ID = as.character(1:n), datB,
                          k = k, padding = 0, q = 2, l = l, password = PWs)
  
  write.table(testDataBF, "BFA.csv", sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)
  write.table(testDataBFB, "BFB.csv", sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)
  
  # Empirical fill and k/l
  fill <- mean(nchar(gsub("0","",testDataBF$CLKs))) / l
  quotient <- k / l
  elements <- k * (32/q)
  
  
  # Link
  multibitTree.load("BFB.csv",threads=7, leafLimit = 3)
  res <- multibitTree.searchFile("BFA.csv", 0.9)

  # Calc quality measures
  tp <- sum(as.character(res$query) == as.character(res$fingerprint))
  fp <- sum(as.character(res$query) != as.character(res$fingerprint))
  fn <- n - tp
  recall <-  (tp/(tp+fn))
  precision <-  (tp/(tp+fp))
  fscore <- (recall + precision) / 2
  
  # Print result
  cat("\n",style("\t    PPRL Test\t\t", bg = "black", fg = "white", font="bold"),
      "\n Positives:\t\t", n,"\n k/i:\t\t\t", k, "/", i,"\n Fill:\t\t\t", fill,
      "\n True positives:\t", tp, 
      "\n False positives:\t", fp, "\n False Negatives:\t", fn,
      "\n Recall:\t\t", (tp/(tp+fn)),"\n Precision:\t\t", (tp/(tp+fp)),
      "\n",style("\t\t   \t\t", bg = "black", fg = "white", font="bold"))
  
  # Append results
  results <- rbind(results, data.frame(data="Synthetic",k,i, fill,elements, quotient, tp, fp, fn, precision, recall, fscore))
 }
}

write.table(results,"./results/Result_Variation_F.csv", row.names=FALSE, sep=";")

results <- read.csv("./results/Result_Variation_F.csv", head=TRUE, sep = ";", stringsAsFactors = FALSE)

aggregate(results$fscore, by=list(results$k), ci)



ggplot(results, aes(fscore)) + geom_dotplot(method="dotdensity", binwidth=0.0005, dotsize=1.1, stackratio = 1.5) +   facet_grid(k ~ ., labeller = label_both) +
  geom_vline(data=filter(results, k==5), aes(xintercept=mean(results$fscore[results$k==5])), colour="#004c93", size=2) + 
  geom_vline(data=filter(results, k==10), aes(xintercept=mean(results$fscore[results$k==10])), colour="#004c93", size=2) + 
  geom_vline(data=filter(results, k==15), aes(xintercept=mean(results$fscore[results$k==15])), colour="#004c93", size=2) + 
  ggtitle("MPR by k") + xlab("Mean prec./rec.") + ylab("Fraction") +
  theme_bw() +  theme(panel.spacing = unit(0.85, "lines"), text = element_text(size=35),legend.justification=c(0.9,0.05), legend.key = element_blank(), legend.key.size = unit(5,"char"), strip.background = element_rect(colour='white',fill = 'white'))
ggsave("./results/Variation_F.pdf", height=18, width=13, units="in", device=cairo_pdf) # Fig 3.11



cairo_pdf(file = "Variation_Normality.pdf", width=9, height=15) # Fig 3.12

par(mfrow=c(3,1), mar=c(5,5,3,1))

qqnorm(results$fscore[results$k==5], main = "QQ-Plot, k=5",   cex=2,  cex.axis=2,cex.lab=2,cex.main=2.2)
qqline(results$fscore[results$k==5])
qqnorm(results$fscore[results$k==10], main = "QQ-Plot, k=10", cex=2,  cex.axis=2,cex.lab=2,cex.main=2.2)
qqline(results$fscore[results$k==10])
qqnorm(results$fscore[results$k==15], main = "QQ-Plot, k=15", cex=2,  cex.axis=2,cex.lab=2,cex.main=2.2)
qqline(results$fscore[results$k==15])

dev.off()



#### Real data ####

#### FEBRL

# Sim data by hashing sampled letters
dat <- read.csv("./Daten/FEBRL_10000_20_A.csv", head=FALSE, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")
datB <- read.csv("./Daten/FEBRL_10000_20_B.csv", head=FALSE, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")

# Generate full string to create BF
dat$linkkey <- paste0(dat$V2, dat$V3, dat$V10)
dat$linkkey <- StandardizeString(dat$linkkey)
datB$linkkey <- paste0(datB$V2, datB$V3, datB$V10)
datB$linkkey <- StandardizeString(datB$linkkey)

l <- 256
q <- 2
n <- nrow(dat)


# Init result matrix
#results <- NULL

# Look over k
for (k in seq(from = 1, to = 30)){
  
  print(k)

  # Gen BF
  testDataBF <- CreateBF(ID = dat$V1, dat$linkkey,
                         k = k, padding = 0, q = 2, l = l, password = "(H]$6Uh04q")
  print("First BF")
  
  testDataBFB <- CreateBF(ID = datB$V1, datB$linkkey,
                          k = k, padding = 0, q = 2, l = l, password = "(H]$6Uh04q")
  print("sec BF")
  
  write.table(testDataBF, "BFA.csv", sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)
  write.table(testDataBFB, "BFB.csv", sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)
  
  # Empirical fill and k/l
  fill <- mean(nchar(gsub("0","",testDataBF$CLKs))) / l
  quotient <- k / l
  elements <- k * (mean(nchar(dat$linkkey))/q)
  
  
  # Link
  multibitTree.load("BFB.csv",threads=7, leafLimit = 3)
  res <- multibitTree.searchFile("BFA.csv", 0.85)
  
  
  
  # Calc quality measures
  tp <- sum(as.character(res$query) == as.character(res$fingerprint))
  fp <- sum(as.character(res$query) != as.character(res$fingerprint))
  fn <- n - tp
  recall <-  (tp/(tp+fn))
  precision <-  (tp/(tp+fp))
  fscore <- (recall + precision) / 2
  print("measures")
  
  rm(res)
  gc()
  
  # Print result
  cat("\n",style("\t    PPRL Test\t\t", bg = "black", fg = "white", font="bold"),
      "\n Positives:\t\t", n,"\n k:\t\t\t", k,"\n Fill:\t\t\t", fill,
      "\n True positives:\t", tp, 
      "\n False positives:\t", fp, "\n False Negatives:\t", fn,
      "\n Recall:\t\t", (tp/(tp+fn)),"\n Precision:\t\t", (tp/(tp+fp)),
      "\n",style("\t\t   \t\t", bg = "black", fg = "white", font="bold"))
  
  # Append results
  results <- rbind(results, data.frame(data="FEBRL",k, q, l, fill, elements, quotient, tp, fp, fn, precision, recall, fscore))
  print("res + next loop")
  
}



#### TelCD


# Sim data by hashing sampled letters
dat <- read.csv("./Daten/A_10000_python_R0_C0_O100_percent.csv", head=FALSE, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")
datB <- read.csv("./Daten/B_10000_python_R20_C50_O100_percent.csv", head=FALSE, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")

# Generate full string to create BF
dat$linkkey <- paste0(dat$V3, dat$V4, dat$V5, dat$V7, dat$V8)
dat$linkkey <- StandardizeString(dat$linkkey)
datB$linkkey <- paste0(datB$V3, datB$V4, datB$V5, datB$V7, datB$V8)
datB$linkkey <- StandardizeString(datB$linkkey)

l <- 256
q <- 2
n <- nrow(dat)


# Init result matrix
#results <- NULL

# Look over k
for (k in seq(from = 1, to = 30)){
  
  
  # Gen BF
  testDataBF <- CreateBF(ID = dat$V1, dat$linkkey,
                         k = k, padding = 0, q = 2, l = l, password = "(H]$6Uh04q")
  
  testDataBFB <- CreateBF(ID = datB$V1, datB$linkkey,
                          k = k, padding = 0, q = 2, l = l, password = "(H]$6Uh04q")
  
  
  write.table(testDataBF, "BFA.csv", sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)
  write.table(testDataBFB, "BFB.csv", sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)
  
  # Empirical fill and k/l
  fill <- mean(nchar(gsub("0","",testDataBF$CLKs))) / l
  quotient <- k / l
  elements <- k * (mean(nchar(dat$linkkey))/q)
  
  
  # Link
  multibitTree.load("BFB.csv",threads=7, leafLimit = 3)
  res <- multibitTree.searchFile("BFA.csv", 0.85)
  
  
  
  # Calc quality measures
  tp <- sum(as.character(res$query) == as.character(res$fingerprint))
  fp <- sum(as.character(res$query) != as.character(res$fingerprint))
  fn <- n - tp
  recall <-  (tp/(tp+fn))
  precision <-  (tp/(tp+fp))
  fscore <- (recall + precision) / 2
  
  # Print result
  cat("\n",style("\t    PPRL Test\t\t", bg = "black", fg = "white", font="bold"),
      "\n Positives:\t\t", n,"\n k:\t\t\t", k,"\n Fill:\t\t\t", fill,
      "\n True positives:\t", tp, 
      "\n False positives:\t", fp, "\n False Negatives:\t", fn,
      "\n Recall:\t\t", (tp/(tp+fn)),"\n Precision:\t\t", (tp/(tp+fp)),
      "\n",style("\t\t   \t\t", bg = "black", fg = "white", font="bold"))
  
  # Append results
  results <- rbind(results, data.frame(data="TelCD",k, q, l, fill, elements, quotient, tp, fp, fn, precision, recall, fscore))
}


# Plotting
#plotdata <- results
write.table(results, "./results/fill_k_quality_full.csv", sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)

plotdata <- read.csv("./results/fill_k_quality_full.csv", head=TRUE, sep = "\t", stringsAsFactors = FALSE)

# No. of elements hashed. Fastdigest returns 32 bit hashes, divided into q-grams, hashes k times.
plotdata$elements <- plotdata$k * (32/plotdata$q)


library(ggplot2)

p1 <- ggplot(data = plotdata, aes(y = round(fill * 100), x = k))
p1 <- p1 + geom_smooth(se=F, method="loess", color="#004c93") #  method="loess"
p1 <- p1 + geom_point(size=6) #  
p1 <- p1 + xlab("\nHash functions k") + ylab("% of bit positions set to one\n")
p1 <- p1 + scale_y_continuous(breaks = seq(0,100,10), limits=  c(1,100), expand = c(0.02, 0))
p1 <- p1 + scale_x_continuous(breaks = seq(0,30,5))
p1 +  facet_grid(data ~ ., labeller = label_both)   + theme_bw() +  theme(panel.spacing = unit(0.85, "lines"), text = element_text(size=35),legend.justification=c(0.9,0.05), legend.key = element_blank(), legend.key.size = unit(5,"char"), strip.background = element_rect(colour='white',fill = 'white'))
ggsave("./results/k_vs_fill_data.pdf", height=11, width=15, units="in", device=cairo_pdf) # Fig 3.20

p2 <- ggplot(data = plotdata, aes(y = fscore, x = round(fill * 100)))
p2 <- p2 + geom_smooth(se=F, method="loess", color="#004c93", span=.25) #  method="loess"
p2 <- p2 + geom_point(size=6) #  methos="loess"
p2 <- p2 + geom_vline(xintercept=50, colour="#004c93") +
  geom_text(aes(x=50, label="\nTheoretical optimum", y=0.8), colour="#004c93", angle=90, size=6)
p2 <- p2 + scale_x_continuous(breaks = seq(0,100,10), limits=  c(1,100), expand = c(0.02, 0))
p2 <- p2 + xlab("\n% of bit positions set to one") + ylab("F-Score\n")
p2  +  facet_grid(data ~ ., labeller = label_both)  + theme_bw() +  theme(panel.spacing = unit(0.85, "lines"), text = element_text(size=35),legend.justification=c(0.9,0.05), legend.key = element_blank(), legend.key.size = unit(5,"char"), strip.background = element_rect(colour='white',fill = 'white'))
ggsave("./results/fill_vs_f_data.pdf", height=11, width=15, units="in", device=cairo_pdf) # Fig 3.21

p3 <- ggplot(data = plotdata, aes(y = fscore, x = k))
p3 <- p3 + geom_smooth(se=F, method="loess", color="#004c93", span=.25) #  method="loess"
p3 <- p3 + geom_point(size=6)  
#p3 <- p3 + geom_point(aes(y = 0.976, x = 8), color="#004c93", size=6)  
#p3 <- p3 + geom_text(aes(y = 1.02, x = 8, label="Optimal k by Smith (2012)"), color="#004c93", size=6)  
p3 <- p3 + xlab("\nHash functions k") + ylab("F-Score\n")
p3 +  facet_grid(data ~ ., labeller = label_both)   + theme_bw() +  theme(panel.spacing = unit(0.85, "lines"), text = element_text(size=35),legend.justification=c(0.9,0.05), legend.key = element_blank(), legend.key.size = unit(5,"char"), strip.background = element_rect(colour='white',fill = 'white'))
ggsave("./results/k_vs_f_data.pdf", height=11, width=15, units="in", device=cairo_pdf)


p4 <- ggplot(data = plotdata, aes(y = fscore, x = elements))
p4 <- p4 + geom_smooth(se=F, method="loess", color="#004c93", span=.25) #  method="loess"
p4 <- p4 + geom_point(size=6)  
p4 <- p4 + xlab("\nNo. of Elements hashed") + ylab("F-Score\n")
p4 +  facet_grid(data ~ ., labeller = label_both)   + theme_bw() +  theme(panel.spacing = unit(0.85, "lines"), text = element_text(size=35),legend.justification=c(0.9,0.05), legend.key = element_blank(), legend.key.size = unit(5,"char"), strip.background = element_rect(colour='white',fill = 'white'))
ggsave("./results/elements_vs_f_data.pdf", height=11, width=15, units="in", device=cairo_pdf)

# AMount of hash functions for 50% fill acc to smith 2012:
(log(1-0.5) / log(1-(1/256))) / 16
# [1] 11.06868

# Opt hash functions for 50% fill acc to smith 2012:
#(l/n) * log(2)
(256/16) * log(2)
#[1] 11.09035

## TelCD
#mean(nchar(dat$linkkey))
# 20.0304

# AMount of hash functions for 50% fill acc to smith 2012:
(log(1-0.5) / log(1-(1/256))) / 20
# [1] 8.854944

# Opt hash functions for 50% fill acc to smith 2012:
#(l/n) * log(2)
(256/20) * log(2)
#[1] 8.872284


## FEBRL
#mean(nchar(dat$linkkey))
# 20.4278

# AMount of hash functions for 50% fill acc to smith 2012:
(log(1-0.5) / log(1-(1/256))) / 20
# [1] 8.854944

# Opt hash functions for 50% fill acc to smith 2012:
#(l/n) * log(2)
(256/20) * log(2)
#[1] 8.872284


#### Random forests vs. linear regression plot

library(ggplot2)
library(randomForest)

# Random numbers
x = runif(1000, min = 0, max = 100)
# Log of x
y = log(x) + runif(1000, min=-0.2, max=0.2)

# Function for plotting
logfun <- function(x) log(x)

# Linear prediction
lin <- lm(y ~ x)
linPred <- predict(lin)

# RF prediction
r2 <- randomForest(y ~ x, ntree=3000, nodesize=5, importance=TRUE, nPerm=5, keep.forest=TRUE)
rfPred <- predict(r2)

# Put all into DF
df <- data.frame(x,y, linPred, rfPred)

# Plot and save
ggplot(df, aes(x=x, y=y)) +                 
  geom_point(aes(color="orange"), size = 3, alpha=1) +  
  stat_function(fun = logfun, color = "black", alpha = 1) +
  geom_point(aes(x=x, y=linPred, color="blue"), size = 4, alpha=.2) +
  geom_point(aes(x=x, y=rfPred, color="red"), size = 4, alpha=.4) +
  xlab("Uniform random numbers") + ylab("log(x)") +
  scale_color_manual(name="Prediction", values=c("#991607", "#d38e0e","#004c93"), labels = c("Linear prediction", "True values","Random Forest prediction")) +
  theme_bw() +   theme(panel.spacing = unit(1, "lines"), text = element_text(size=27), 
                       axis.text.x = element_text(size = 21, angle = 0, vjust = 0.5),  axis.title=element_text(size=30),  
                       legend.justification=c(0.9,0.05), legend.key = element_blank(), legend.key.size = unit(5,"char"), 
                       strip.background = element_rect(colour='white',fill = 'white'))
ggsave("linear_vs_RF.pdf",height=13, width=18, units="in", device=cairo_pdf) # Fig 4.11


