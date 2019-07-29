###################### PLOTTING ##########################################


# Load libs
library(ggplot2)
library(grid)
library(gridExtra)
library(colorspace)
library(RColorBrewer)
library(rsm)
library(boot)
library(randomForest)
library(stargazer)
library(xtable)
library(dplyr)
library(effects)
library(cvTools)
library(car)
library(gmodels)
library(MASS)
library(devtools)
install_github("markwestcott34/stargazer-booktabs")

#setwd("E:/Priv/Dropbox/Dissertation/Programme")


## Bootstrap mean function
boot.mean = function(x,B,confidence = 0.95,binwidth = NULL){
  require(boot)
  n = length(x)
  boot.samples = matrix( sample(x,size = n*B,replace = TRUE), B, n)
  boot.statistics = apply(boot.samples,1,mean)
  se = sd(boot.statistics)
  require(ggplot2)
  if ( is.null(binwidth) )
    binwidth = diff(range(boot.statistics))/30
  p = ggplot(data.frame(x = boot.statistics),aes(x = x)) +geom_histogram(aes(y=..density..),binwidth = binwidth) + geom_density(color = "red")
  #plot(p)
  ciLow = mean(x) + -1*qnorm(1-(confidence/2))*se
  ciHigh = mean(x) + 1*qnorm(1-(confidence/2))*se
 # print( ciLow, ciHigh )
  return(data.frame(ciLow = ciLow, ciHigh = ciHigh, se = se) )
}

# Datasets list
datasets <- c("German Mortality", "NC Voter", "German Telephone CD 20% errors", "German Telephone CD 0% errors", "FEBRL", "FEBRL WA", "Mortality Test Data")



# # Read data
plotdata <- as.data.frame(read.csv("./Results/R_Metainfo_Testbett_CompareOPC_vs_Train.csv", sep = "\t"), stringsAsFactors = FALSE)
 

# Make sure central measures are numeric
plotdata$recall <- as.numeric(plotdata$recall)
plotdata$precision <- as.numeric(plotdata$precision)
plotdata$Fmeasure <- as.numeric(plotdata$Fmeasure)


# Read training file for model building
traindata <- read.csv("Complete_Results.csv", sep = "\t")

# Subset to fixed tanimoto threshold
fixed <- subset(traindata, Tani == 0.9)


# Build new complete data set with evaluation data results for rapid model comparisons
trainresultsFull <- rbind(read.csv("Trainresults_FEBRLWA.csv", sep = "\t"),read.csv("Trainresults_Mortality.csv", sep = "\t"))
trainresultsFull$l <- 500
trainresultsFull$missingamount.1 <- NULL

trainresultsFull <- rbind(trainresultsFull, read.csv("Trainresults_NCVoter.csv", sep = "\t"))
trainresultsFull$type <- NULL

kill <- c("tp", "fp",	"fn",	"tn",	"recall",	"precision")
fulldata <- fixed[,(!names(fixed) %in% kill)]
fulldata$variant <- NULL
fulldata$positives <- NULL
fulldata$i <- 1

fulldata <- rbind(fulldata,trainresultsFull)
fixed <- subset(fulldata, Tani == 0.9)


# RSM Plot
# RSM: Preproc data
RSMdata <- subset(fulldata, Tani == 0.9)[,c("k", "q", "u", "meanengrams", "hammingweight", "Fmeasure")]
RSMdata_unscaled <- subset(fulldata, Tani == 0.9)[,c("k", "q", "u", "meanengrams", "hammingweight", "Fmeasure")]

# Scale min/max
maxi = apply(RSMdata , 2, function(x) max(as.numeric(x)))
mini = apply(RSMdata,  2, function(x) min(as.numeric(x)))
# Scaling
RSMdata = as.data.frame(scale(RSMdata, center = mini, scale = maxi - mini))
# Train model
rsm.train <- rsm(Fmeasure ~   SO(k, q) + TWI(q, k, hammingweight) + u, data = RSMdata)



# Fixed model calls
fullmod <- lm(Fmeasure ~  k + (q) + errorestimation +  skew + uniqueness + meanentropy + log(u) + ginicoefficient +  uniquepatternsA + errorestimation +  k * q *  hammingweight +  meanengrams*hammingweight + 
              meanmissing +  log(m) + u*m, 
              data = fixed)

minmod <- lm(Fmeasure ~ k + (q), data = fixed)



stepTest <- stepAIC(fullmod, minmod, direction = "both",  k = log(nrow(fixed))) # bic
#stepTest2 <- stepAIC(fullmod, minmod, direction = "both",  k = 2) # aic
BIC2 <- stepTest


# LM
LM <- lm(Fmeasure ~ k + (q) + errorestimation * skew + uniquepatternsA * 
           uniqueness + ginicoefficient + k *  hammingweight + 
           hammingweight * meanengrams + uniquepatternsA + log(u),
         data = fixed)


# Random forest
r1 <- randomForest(Fmeasure ~ k  +  q + errorestimation +  skew + uniqueness + ginicoefficient + hammingweight + 
                     meanengrams + u + uniquepatternsA, data = fixed, 
                   mtry= 3, ntree = 500, nodesize = 20, 
                   importance = TRUE, nPerm = 5, keep.forest = TRUE)


fileConn<-file("Models.tex")

# Model regression tables
writeLines(stargazer(LM, BIC2, title = "Regression models",align = TRUE, omit.stat = c("f"), no.space = TRUE, type = "latex",
          column.labels = c("Simple linear model", "BIC-selected linear model"), digits = 3, digits.extra = 0,
          ci = TRUE,
          ci.separator = ";",
          star.char = c("+","*"),
          star.cutoffs = c(0.05,0.001),
          notes = c("$^{+}$p$<$0.05", "$^{*}$p$<$0.001"), 
          notes.append = FALSE,
covariate.labels = c("k", "q", "errorestimation",
                     "skewness", "uniquepattA", "uniqueness", "meanentropy", "gini",
                     "hammingweight", "meanngr", "meanmissing", "log(m)", "k*q", "log(u)",
                     "errorest*skewness", "uniquepattA*uniqueness", "k*hammingweight",
                     "hammingweight*meanngr", "q*hammingweight", "k*q*hammingweight"),
dep.var.caption  = "Mean of precision and recall",
dep.var.labels   = "MPR")[17:63], fileConn)
close(fileConn)


# RSM Surfaces
cairo_pdf(file = "SO_RSM.pdf", width = 20, height = 10) # Fig 4.8
par(mfrow = c(1,2))
contour(rsm.train, ~ hammingweight + k, image = TRUE, main = "Second-order RS contour plot",  xlabs = c("k Hash functions", "q-grams"), cex.lab = 1.3, cex.axis = 1.2, labcex = 1.4)
persp(rsm.train, hammingweight ~ k, xlabs = c("\nStandard. k Hash functions", "\nStandard. q-grams"), 
      zlab = "\nPredicted Standard. MPR", main = "Second-order response surface", theta = 150, phi = 30, cex.lab = 1.6, cex.axis = 1.2, contours = TRUE)
dev.off()


cairo_pdf(file = "testResult_BIC_Multidata.pdf", width=20, height=10) # Fig 5.2
par( mfrow=c(1, 2) )
fullmod <- lm(Fmeasure ~  k + (q) + Tani + errorestimation +  skew + uniqueness + meanentropy + log(u) + ginicoefficient +  uniquepatternsA + errorestimation +  k * q *  hammingweight +  meanengrams*hammingweight + 
                meanmissing +  log(m) + u*m, 
              data = subset(traindata, data=="FEBRL"))

persp(fullmod, k ~ Tani, zlab = "\n\nMPR", theta = -120, phi = 25,
      col="gray40", border="gray10", main = "FEBRL",
      xlabs =rev(c("\nTanimoto threshold","\nHash functions")))

fullmod <- lm(Fmeasure ~  k + (q) + Tani + errorestimation +  skew + uniqueness + meanentropy + log(u) + ginicoefficient +  uniquepatternsA + errorestimation +  k * q *  hammingweight +  meanengrams*hammingweight + 
                meanmissing +  log(m) + u*m, 
              data = subset(traindata, data=="German Mortality"))
persp(fullmod, k ~ Tani , zlab = "\n\nMPR", theta = -120, phi = 25,
      col="gray40", border="gray10", main = "Mortality & Hospital data",
      xlabs =rev(c("\nTanimoto threshold","\nHash functions")))
dev.off()



# Fit plot BIC
newdata <- subset(fulldata, Tani == 0.9 & data == "NC Voter" & q == 1)
newfittedBIC <- data.frame(newdata, predicted = predict(BIC2, newdata, type = "response"))

cairo_pdf(file = "Predicted_vs_F.pdf", width = 18, height = 15)  # Fig 5.7
par(mgp = c(3,1,5), mar = c(4.5,4.5,1,1))
plot(newfittedBIC$Fmeasure ~ newfittedBIC$predicted,  ylab = "\nF-Measure", xlab = "Predicted F-measure", cex = 2, xlim = c(0.84,0.92), ylim = c(0.84,0.92), 
     cex.axis = 2,cex.lab = 2,cex.main = 2.2) # xlim = c(0.84,0.95), ylim = c(0.84,0.95)
segments(0,0,1,1, lty = "dotted")
dev.off()


# Fit plot RF
newdata <- subset(fulldata, Tani == 0.9 & data == "NC Voter")
newdata$fitForest <- predict(r1, newdata)
ggplot(newdata, aes(y = Fmeasure, x = fitForest, color = factor(q), group = q)) + 
  xlab("RF-predicted mean prec./rec.") + ylab("Real mean prec./rec") +
  scale_colour_manual(name = "q", labels= as.character(1:4), values = c("#004c93", "#56B1F7", "#c43f4b", "#000000")) +
  geom_point(size = 6) + geom_line(size = 1) + 
  # coord_fixed(xlim = c(0.82, 0.97), ylim = c(0.82,0.97)) + 
  theme_bw() + theme(panel.spacing = unit(0.25, 'lines'), 
        text = element_text(size = 22), 
        axis.title = element_text(size = 28),  
        legend.justification = c(0.9,0.05), 
        legend.key = element_blank(), 
        legend.key.size = unit(5,"char"), 
        strip.background = element_rect(colour='white',fill = 'white')) + 
  ggtitle(paste("Data:", unique(newdata$data)))
ggsave("RF_predict_vs_fit.pdf", width = 15, height = 12, unit = "in", device=cairo_pdf)  # Fig 5.8






# Crossvalidate
crossvalids <- NULL
for (d in datasets){
  fixedCV <- subset(fixed, data == d)
  fixedCV$kx <- fixedCV$k
  fixedCV$hw <- fixedCV$hammingweight
  cvBICs <- cvFit (BIC2, data = fixedCV, y = fixedCV$Fmeasure,
                   cost = mape, K = 5, R = 30, includeSE = TRUE)
  cvLMs  <- cvFit (LM, data = fixedCV, y = fixedCV$Fmeasure,
                   cost = mape, K = 5, R = 30, includeSE = TRUE)
  cvRSMs <- cvFit (rsm.train, data = fixedCV, y = fixedCV$Fmeasure,
                   cost = mape, K = 5, R = 30, includeSE = TRUE)  
  cvRFs <- cvFit (r1, data = fixedCV, y = fixedCV$Fmeasure,
                  cost = mape, K = 5, R = 30, includeSE = TRUE)  
  crossvalids <- rbind(crossvalids, data.frame(data = d, 
                        cvBIC = cvBICs[["cv"]][["CV"]], seBIC = cvBICs[["se"]][["CV"]], 
                        cvLM = cvLMs[["cv"]][["CV"]], seLM = cvLMs[["se"]][["CV"]], 
                        cvRSM = cvRSMs[["cv"]][["CV"]], seRSM = cvRSMs[["se"]][["CV"]],
                        cvRF = cvRFs[["cv"]][["CV"]], seRF = cvRFs[["se"]][["CV"]]))
}

# Reorder
crossvalids <- crossvalids[c(4,3,5,1,2,7,6),]

# Write all diagnostics to tex file
fileConn<-file("Diagnostics.tex")

# LM
sink(fileConn)

print(xtable(
  data.frame(data = aggregate(sub$Fmeasure[sub$variant == "Simple linear model"], by = list(sub$data[sub$variant == "Simple linear model"]), 
             mean, na.action = na.exclude)$Group.1, meanF = aggregate(sub$Fmeasure[sub$variant == "Simple linear model"], 
             by = list(sub$data[sub$variant == "Simple linear model"]), 
             mean, na.action= na.exclude)$x,  cv = crossvalids$cvLM, se = crossvalids$seLM, 
             RMSE = aggregate(sub$RMSE[sub$variant == "Simple linear model"], 
             by = list(sub$data[sub$variant == "Simple linear model"]), 
             mean, na.action= na.exclude)$x),
  digits = 3, display = c("s", "s", "f", "e", "e", "f")),math.style.exponents = TRUE, booktabs = TRUE, include.rownames = FALSE)

sink(fileConn, append = TRUE)

# BIC
print(xtable(
  data.frame(data = aggregate(sub$Fmeasure[sub$variant == "BIC-selected linear model"], 
             by = list(sub$data[sub$variant == "BIC-selected linear model"]), mean, na.action= na.exclude)$Group.1,
             meanF = aggregate(sub$Fmeasure[sub$variant == "BIC-selected linear model"], 
             by = list(sub$data[sub$variant == "BIC-selected linear model"]), mean, na.action= na.exclude)$x,  
             cv = crossvalids$cvBIC, se = crossvalids$seBIC, RMSE = aggregate(sub$RMSE[sub$variant == "BIC-selected linear model"], 
             by = list(sub$data[sub$variant == "BIC-selected linear model"]), 
             mean, na.action= na.exclude)$x),
  digits = 3, display = c("s", "s", "f", "e", "e", "f")),math.style.exponents = TRUE, booktabs = TRUE, include.rownames = FALSE)

sink(fileConn, append = TRUE)

# RSM
print(xtable(
  data.frame(data = aggregate(sub$Fmeasure[sub$variant == "Optimal Parameter choice (RSM)"], 
            by = list(sub$data[sub$variant == "Optimal Parameter choice (RSM)"]), mean, na.action= na.exclude)$Group.1,
            meanF = aggregate(sub$Fmeasure[sub$variant == "Optimal Parameter choice (RSM)"], 
            by = list(sub$data[sub$variant == "Optimal Parameter choice (RSM)"]), mean, na.action= na.exclude)$x, 
            cv = crossvalids$cvRSM, se = crossvalids$seRSM, 
            RMSE = aggregate(sub$RMSE[sub$variant == "Optimal Parameter choice (RSM)"], 
            by = list(sub$data[sub$variant == "Optimal Parameter choice (RSM)"]), mean, na.action= na.exclude)$x),
  digits = 3, display = c("s", "s", "f", "e", "e", "f")),math.style.exponents = TRUE, booktabs = TRUE, include.rownames = FALSE)


sink(fileConn, append = TRUE)





# RF
print(xtable(
  data.frame(data = aggregate(sub$Fmeasure[sub$variant == "Optimal Parameter choice (Forests)"], 
              by = list(sub$data[sub$variant == "Optimal Parameter choice (Forests)"]), 
              mean, na.action= na.exclude)$Group.1,meanF = aggregate(sub$Fmeasure[sub$variant == "Optimal Parameter choice (Forests)"], 
              by = list(sub$data[sub$variant == "Optimal Parameter choice (Forests)"]), 
              mean, na.action= na.exclude)$x,  cv = crossvalids$cvRF, se = crossvalids$seRF, 
              RMSE = aggregate(sub$RMSE[sub$variant == "Optimal Parameter choice (Forests)"], 
              by = list(sub$data[sub$variant == "Optimal Parameter choice (Forests)"]), mean, na.action= na.exclude)$x),
  digits = 3, display = c("s", "s", "f", "e", "e", "f")),math.style.exponents = TRUE, booktabs = TRUE, include.rownames = FALSE)

sink()
close(fileConn)



# Plots by measure

# Set colors, preprocess data
set.seed(19)
cbbPalette <- (c("#000000", "#2643dc",
                 "#5ca03f",
                 "#344494",
                 "#c43f4b",
                 "#ce5d2c", "#555555"))


set.seed(19)
shapes <- sample(21:25, length(unique(plotdata$variant)), replace = TRUE)

plotdata$data <- factor(plotdata$data, 
                        levels = c("FEBRL WA", "Mortality Test Data", "FEBRL", "German Mortality", "German Telephone CD 0% errors", "German Telephone CD 20% errors", "NC Voter"), 
                        labels = c("FEBRL WA", "Mortality & Commercial", "FEBRL", "Mortality & Hospital", "Telephone CD 0% errors", "Telephone CD 20% errors", "NC Voter"))

plotdata$variant <- factor(plotdata$variant, 
                           levels = c(  "Optimal Parameter choice (RSM)",  "Optimal Parameter choice (RF)",  "Optimal Parameter choice (NN)", "Optimal Parameter choice (BIC)",  
                                        "Optimal Parameter choice (LM)", "50%-Rule Parameter choice", "Standard CLK k = 15" , "Standard CLK k = 10"  ,"Standard CLK k = 5"),
labels = c(  "Optimal Parameter choice (RSM)",  "Optimal Parameter choice (Forests)",  "Optimal Parameter choice (NN)", "BIC-selected linear model",  
             "Simple linear model", "50%-Rule Parameter choice", "Standard CLK k = 15" , "Standard CLK k = 10"   ,"Standard CLK k = 5"))

# No NN, no k = 15
plotdata <- subset(plotdata, plotdata$variant != "Optimal Parameter choice (NN)" & plotdata$variant != "Standard CLK k = 15")
plotdata$variant <- droplevels(plotdata$variant)


# Main Plot: F-Measure of all variants by Tanimoto-Threshold and errors
p1 <- ggplot(subset(plotdata, data == "Mortality & Commercial" | data =="NC Voter" | data == "FEBRL WA" ))
p1 <- p1 + geom_line(aes(y = Fmeasure, x = Tani,group = variant, colour = variant), size = 1) 
p1 <- p1 + geom_point(aes(y = Fmeasure,x = Tani,group = variant, shape = variant, colour = variant, fill = variant), size = 7)
p1 <- p1 + xlab("\nTanimoto threshold") + ylab("Mean Prec./Rec.\n") 
p1 <- p1 + scale_colour_manual(name = "Variant", labels= levels(plotdata$variant), values = cbbPalette) 
p1 <- p1 + scale_fill_manual(name = "Variant", labels= levels(plotdata$variant), values = cbbPalette) 
p1 <- p1 + scale_shape_manual(name = "Variant", labels= levels(plotdata$variant), values = shapes)
p1 <- p1 + ggtitle("Test data")

p1  + facet_grid(. ~ data) + theme_bw() +  theme(panel.spacing = unit(0.85, "lines"), 
       text = element_text(size = 27), axis.text.x = element_text(size = 23, angle = 90, vjust = 0.5), 
       axis.title = element_text(size = 30),  legend.justification = c(0.9,0.05), 
       legend.key = element_blank(), legend.key.size = unit(5,"char"),
       strip.background = element_rect(colour='white',fill = 'white'))
# Save Plot
ggsave("./results/Testbett_EncryptionResults_F-Score_Testdata_experimental.pdf", height = 11, width = 21, units = "in", device=cairo_pdf) #  Fig 5.14


# Main Plot: F-Measure of all variants by Tanimoto-Threshold and errors
p1 <- ggplot(subset(plotdata, data != "Mortality & Commercial" & data !="NC Voter" & data != "FEBRL WA" ))
p1 <- p1 + geom_line(aes(y = Fmeasure, x = Tani,group = variant, colour = variant), size = 1) 
p1 <- p1 + geom_point(aes(y = Fmeasure,x = Tani,group = variant, shape = variant, colour = variant, fill = variant), size = 7)
p1 <- p1 + xlab("\nTanimoto threshold") + ylab("Mean Prec./Rec.\n") 
p1 <- p1 + scale_colour_manual(name = "Variant", labels= levels(plotdata$variant), values = cbbPalette) 
p1 <- p1 + scale_fill_manual(name = "Variant", labels= levels(plotdata$variant), values = cbbPalette) 
p1 <- p1 + scale_shape_manual(name = "Variant", labels= levels(plotdata$variant), values = shapes)
p1 <- p1 + ggtitle("Training data")
p1  + facet_grid(. ~ data) + theme_bw() +  theme(panel.spacing = unit(0.85, "lines"), text = element_text(size = 27), axis.text.x = element_text(size = 23, angle = 90, vjust = 0.5),  
       axis.title = element_text(size = 30),  legend.justification = c(0.9,0.05), 
       legend.key = element_blank(), legend.key.size = unit(5,"char"), strip.background = element_rect(colour='white',fill = 'white'))
# Save Plot
ggsave("./results/Testbett_EncryptionResults_F-Score_experimental.pdf", height = 11, width = 21, units = "in", device=cairo_pdf)  #  Fig 5.13


# Sub Plot: Recall of all variants by Tanimoto-Threshold and errors
p2 <-  ggplot(subset(plotdata, data != "Mortality & Commercial" & data !="NC Voter" & data != "FEBRL WA" ))
p2 <- p2 + geom_line(aes(y = recall, x = Tani,group = variant, colour = variant), size = 1) 
p2 <- p2 + geom_point(aes(y = recall,x = Tani,group = variant, shape = variant, colour = variant, fill = variant), size = 7)
p2 <- p2 + xlab("\nTanimoto threshold") + ylab("Recall\n") 
p2 <- p2 + scale_colour_manual(name = "Variant", labels= levels(plotdata$variant), values = cbbPalette) 
p2 <- p2 +  scale_fill_manual(name = "Variant", labels= levels(plotdata$variant), values = cbbPalette) 
p2 <- p2 + scale_shape_manual(name = "Variant", labels= levels(plotdata$variant), values = shapes)
p2 <- p2 + ggtitle("Training data")
p2  + facet_grid(. ~ data) + theme_bw() +  theme(panel.spacing = unit(0.85, "lines"), text = element_text(size = 27), 
         axis.text.x = element_text(size = 23, angle = 90, vjust = 0.5),  axis.title = element_text(size = 30),  
         legend.justification = c(0.9,0.05), legend.key = element_blank(), legend.key.size = unit(5,"char"), 
         strip.background = element_rect(colour='white',fill = 'white'))
# Save Plot
ggsave("./results/Testbett_EncryptionResults_Recall_experimental.pdf", height = 11, width = 21, units = "in", device=cairo_pdf)  #  Fig 5.9


# Sub Plot: Recall of all variants by Tanimoto-Threshold and errors
p2 <-  ggplot(subset(plotdata, data == "Mortality & Commercial" | data == "NC Voter" | data == "FEBRL WA" ))
p2 <- p2 + geom_line(aes(y = recall, x = Tani,group = variant, colour = variant), size = 1) 
p2 <- p2 + geom_point(aes(y = recall,x = Tani,group = variant, shape = variant, colour = variant, fill = variant), size = 7)
p2 <- p2 + xlab("\nTanimoto threshold") + ylab("Recall\n") 
p2 <- p2 + scale_colour_manual(name = "Variant", labels= levels(plotdata$variant), values = cbbPalette) 
p2 <- p2 +  scale_fill_manual(name = "Variant", labels= levels(plotdata$variant), values = cbbPalette) 
p2 <- p2 + scale_shape_manual(name = "Variant", labels= levels(plotdata$variant), values = shapes)
p2 <- p2 + ggtitle("Test data")
p2  + facet_grid(. ~ data) + theme_bw() +  theme(panel.spacing = unit(0.85, "lines"), text = element_text(size = 27), 
        axis.text.x = element_text(size = 23, angle = 90, vjust = 0.5),  axis.title = element_text(size = 30),  
        legend.justification = c(0.9,0.05), legend.key = element_blank(), legend.key.size = unit(5,"char"), 
        strip.background = element_rect(colour='white',fill = 'white'))
# Save Plot
ggsave("./results/Testbett_EncryptionResults_Recall_experimental_Testdata.pdf", height = 11, width = 21, units = "in", device=cairo_pdf) #  Fig 5.10


# Sub Plot: Precision of all variants by Tanimoto-Threshold and errors
p3 <- ggplot(subset(plotdata, data != "Mortality & Commercial" & data !="NC Voter" & data != "FEBRL WA" ))
p3 <- p3 + geom_line(aes(y = precision, x = Tani,group = variant, colour = variant), size = 1) 
p3 <- p3 + geom_point(aes(y = precision,x = Tani,group = variant, shape = variant, colour = variant, fill = variant), size = 7)
p3 <- p3 + xlab("\nTanimoto threshold") + ylab("Precision\n") 
p3 <- p3 + scale_colour_manual(name = "Variant", labels= levels(plotdata$variant), values = cbbPalette) 
p3 <- p3 +  scale_fill_manual(name = "Variant", labels= levels(plotdata$variant), values = cbbPalette) 
p3 <- p3 + scale_shape_manual(name = "Variant", labels= levels(plotdata$variant), values = shapes)
p3 <- p3 + ggtitle("Training data")
p3  + facet_grid(. ~ data) + theme_bw() +  theme(panel.spacing = unit(0.85, "lines"), 
    text = element_text(size = 27), axis.text.x = element_text(size = 23, angle = 90, vjust = 0.5),  
    axis.title = element_text(size = 30), legend.justification = c(0.9,0.05), 
    legend.key = element_blank(), legend.key.size = unit(5,"char"), 
    strip.background = element_rect(colour='white',fill = 'white'))
# Save Plot
ggsave("./results/Testbett_EncryptionResults_Precision_experimental.pdf", height = 11, width = 21, units = "in", device=cairo_pdf) #  Fig 5.11


# Sub Plot: Precision of all variants by Tanimoto-Threshold and errors
p3 <- ggplot(subset(plotdata, data == "Mortality & Commercial" | data == "NC Voter" | data == "FEBRL WA" ))
p3 <- p3 + geom_line(aes(y = precision, x = Tani,group = variant, colour = variant), size = 1) 
p3 <- p3 + geom_point(aes(y = precision,x = Tani,group = variant, shape = variant, colour = variant, fill = variant), size = 7)
p3 <- p3 + xlab("\nTanimoto threshold") + ylab("Precision\n") 
p3 <- p3 + scale_colour_manual(name = "Variant", labels= levels(plotdata$variant), values = cbbPalette) 
p3 <- p3 +  scale_fill_manual(name = "Variant", labels= levels(plotdata$variant), values = cbbPalette) 
p3 <- p3 + scale_shape_manual(name = "Variant", labels= levels(plotdata$variant), values = shapes)
p3 <- p3 + ggtitle("Test data")
p3  + facet_grid(. ~ data) + theme_bw() +  theme(panel.spacing = unit(0.85, "lines"), text = element_text(size = 27), 
        axis.text.x = element_text(size = 23, angle = 90, vjust = 0.5),  axis.title = element_text(size = 30),  
        legend.justification = c(0.9,0.05), legend.key = element_blank(), legend.key.size = unit(5,"char"), 
        strip.background = element_rect(colour='white',fill = 'white'))
# Save Plot
ggsave("./results/Testbett_EncryptionResults_Precision_experimental_Testdata.pdf", height = 11, width = 21, units = "in", device=cairo_pdf) #  Fig 5.12




## Central dot chart: MPR by method at T = 0.9

sub <- subset(plotdata, Tani == 0.9) # & data != "German Telephone CD 0% errors")
sub$data <- droplevels(sub$data)
sub$variant <- reorder(sub$variant, sub$Fmeasure, mean)
sub$data <- reorder(sub$data, -sub$Fmeasure, min)

sub <- sub[order(sub$Fmeasure),]

# Central Dotchart with ggplot
p4 <- ggplot(sub, aes(x = variant, y = Fmeasure, colour = data))
p4 <- p4 +  geom_point(shape = 16, size = 7) + facet_wrap( ~ data, ncol = 1) + coord_flip() 
p4 <- p4 + xlab("") + ylab("\nMean Prec./Rec.") +    guides(colour = FALSE) + 
  labs(caption = "All values at a Tanimoto threshold of T = 0.9")
p4 <- p4 +  theme_bw() + theme(panel.spacing = unit(0.25, 'lines'), 
                     text = element_text(size = 22), 
                     axis.title = element_text(size = 28),  
                     legend.justification = c(0.9,0.05), 
                     legend.key = element_blank(), 
                     legend.key.size = unit(5,"char"), 
                     strip.background = element_rect(colour='white',fill = 'white'),
                     plot.caption = element_text(size = 16, hjust = -.7, vjust = 0.1),
                     panel.grid.major.x = element_blank() ,
                     panel.grid.minor.x = element_blank() ,
                     axis.text.x = element_text(colour = "#111111"),
                     axis.text.y = element_text(size = 22, colour = "#111111"),
                     ## explicitly set the horizontal lines (or they will disappear too)
                     panel.grid.major.y = element_line(color = "#444444", linetype = 3))
show(p4)
#Save the plot
ggsave("./results/Dotchart_F-Score_OptTrain_experimental.pdf", height = 19, width = 13, units = "in", device=cairo_pdf) # Fig 5.15

 


## Diagnostic plots


# Random forest diagnostics plot

modelname <- r1
fitted <- r1$predicted
truescore <- fixed$Fmeasure
resid <- truescore - fitted
stdresid <- resid / sd(resid)


cooksd <- rep(0, length(fitted))
for (i in 1:length(fitted)){
  
  cooksd[i] <- sum((fitted[-i] - predict(r1, newdata = fixed[-i,]))^2) / var(r1$mse) * 6
  
  
}

model <- data.frame(fitted, resid, stdresid, truescore, cooksd, data = fixed$data, q = factor(fixed$q))

p1 <- ggplot(model, aes(fitted, resid))+geom_point(size = 2, aes(color = q))
p1 <- p1 + stat_smooth(method = "loess")+geom_hline(yintercept = 0, col = "red", linetype = "dashed")
p1 <- p1 + xlab("Fitted values") + ylab("Residuals")
p1 <- p1 + ggtitle("Residual vs Fitted Plot") + theme_bw() + facet_grid(data ~ .) + 
          theme(panel.spacing = unit(0.25, 'lines'), 
          text = element_text(size = 22), 
          axis.title = element_text(size = 28),  
          legend.justification = c(0.9,0.05), 
          legend.key = element_blank(), 
          legend.key.size = unit(5,"char"), 
          strip.text = element_text(size = 8),
          strip.background = element_rect(colour='white',fill = 'white'))


p2 <- ggplot(model, aes(sample = stdresid)) +
  stat_qq(size = 2) +
  stat_qq_line(size = 1) + xlab("Theoretical Quantiles")+ylab("Standardized Residuals") 
p2 <- p2 + ggtitle("Normal Q-Q")+theme_bw() + theme(panel.spacing = unit(0.25, 'lines'), 
                                                text = element_text(size = 22), 
                                                axis.title = element_text(size = 28),  
                                                legend.justification = c(0.9,0.05), 
                                                legend.key = element_blank(), 
                                                legend.key.size = unit(5,"char"), 
                                                strip.background = element_rect(colour='white',fill = 'white'))

p3 <- ggplot(model, aes(fitted, sqrt(abs(stdresid))))+geom_point(aes(color = q),size = 2,  na.rm = TRUE)
p3 <- p3 + stat_smooth(method = "loess", na.rm = TRUE)+xlab("Fitted Value")
p3 <- p3 + ylab(expression(sqrt("|Standardized residuals|")))
p3 <- p3 + ggtitle("Scale-Location")+theme_bw() + facet_grid(data ~ .) + 
      theme(panel.spacing = unit(0.25, 'lines'), 
      text = element_text(size = 22), 
      axis.title = element_text(size = 28),  
      legend.justification = c(0.9,0.05), 
      legend.key = element_blank(), 
      legend.key.size = unit(5,"char"), 
      strip.text = element_text(size = 8),
      strip.background = element_rect(colour='white',fill = 'white'))

subs <-   data.frame(variable = names(importance(modelname)[,2]), importance = importance(modelname)[,2])
subs$variable <-  reorder(subs$variable, subs$importance, mean)

# Dotchart with ggplot
p4 <- ggplot(subs, aes(x = variable, y = importance)) +ggtitle("Variable Importance")+
  geom_point(shape = 16, size = 4) + coord_flip() + xlab("") + ylab("\nImportance")  +
  theme_bw() + theme(panel.spacing = unit(0.25, 'lines'), 
                     text = element_text(size = 22), 
                     axis.title = element_text(size = 28),  
                     legend.justification = c(0.9,0.05), 
                     legend.key = element_blank(), 
                     legend.key.size = unit(5,"char"), 
                     strip.background = element_rect(colour='white',fill = 'white'),
                     plot.caption = element_text(size = 16, hjust = -.7, vjust = 0.1),
                     panel.grid.major.x = element_blank() ,
                     panel.grid.minor.x = element_blank() ,
                     axis.text.x = element_text(colour = "#111111"),
                     axis.text.y = element_text(size = 22, colour = "#111111"),
                     ## explicitly set the horizontal lines (or they will disappear too)
                     panel.grid.major.y = element_line(color = "#444444", linetype = 3))

# Save
g1 <- arrangeGrob(p1,p3, layout_matrix = rbind(c(1),c(2)))
ggsave("./results/Diagnostic_RF.pdf", height = 19, width = 10, units = "in", g1, device=cairo_pdf) # Fig 5.6

# Variable importance
ggsave("./results/RF_VarImp.pdf", height = 11, width = 11, units = "in", p4, device=cairo_pdf) # Fig 4.12


# Function for all other diagnostic plots
diagPlot<-function(model){
  p1 <- ggplot(model, aes(.fitted, .resid))+geom_point()
  p1 <- p1 + stat_smooth(method = "loess")+geom_hline(yintercept = 0, col = "red", linetype = "dashed")
  p1 <- p1 + xlab("Fitted values")+ylab("Residuals")
  p1 <- p1 + ggtitle("Systematic Residuals vs Fitted?") + theme_bw() +  theme(panel.spacing = unit(0.25, 'lines'), 
     text = element_text(size = 18), 
     axis.title = element_text(size = 28),  
     legend.justification = c(0.9,0.05), 
     legend.key = element_blank() )

  
  p2 <- ggplot(model, aes(sample=.stdresid)) + stat_qq() + stat_qq_line() +
    xlab("Theoretical Quantiles")+ylab("Standardized Residuals")
  p2 <- p2 + ggtitle("Q-Q Plot: Normally distributed?") + theme_bw() +theme(panel.spacing = unit(0.25, 'lines'), 
    text = element_text(size = 18), 
    axis.title = element_text(size = 28),  
    legend.justification = c(0.9,0.05), 
    legend.key = element_blank() )
  
  p3 <- ggplot(model, aes(.fitted, sqrt(abs(.stdresid)))) + geom_point(na.rm = TRUE)
  p3 <- p3 + stat_smooth(method = "loess", na.rm = TRUE) + xlab("Fitted Value")
  p3 <- p3 + ylab(expression(sqrt("|Standardized residuals|")))
  p3 <- p3 + ggtitle("Visual test for Homoscedacity") + theme_bw() +  theme(panel.spacing = unit(0.25, 'lines'), 
      text = element_text(size = 18), 
      axis.title = element_text(size = 28),  
      legend.justification = c(0.9,0.05), 
      legend.key = element_blank() )
  
  p4 <- ggplot(model, aes(seq_along(.cooksd), .cooksd)) + geom_bar(stat = "identity", position = "identity")
  p4 <- p4 + xlab("Obs. Number") + ylab("Cook's distance")
  p4 <- p4 + ggtitle("Cook's distance") + theme_bw() +  theme(panel.spacing = unit(0.25, 'lines'), 
      text = element_text(size = 18), 
      axis.title = element_text(size = 28),  
      legend.justification = c(0.9,0.05), 
      legend.key = element_blank() )
  
  p5 <- ggplot(model, aes(.hat, .stdresid))+geom_point(aes(size=.cooksd), na.rm = TRUE)
  p5 <- p5 + stat_smooth(method = "loess", na.rm = TRUE)
  p5 <- p5 + xlab("Leverage")+ylab("Standardized Residuals")
  p5 <- p5 + ggtitle("Leverage plot + Cooks distance")
  p5 <- p5 + scale_size_continuous("Cook's Distance", range = c(1,5))
  p5 <- p5 + theme_bw() + theme(legend.position = "bottom")+  theme(panel.spacing = unit(0.25, 'lines'), 
      text = element_text(size = 18), 
      axis.title = element_text(size = 28),  
      legend.justification = c(0.9,0.05), 
      legend.key = element_blank() )
  
  p6 <- ggplot(model, aes(.hat, .cooksd))+geom_point(na.rm = TRUE)+stat_smooth(method = "loess", na.rm = TRUE)
  p6 <- p6 + xlab("Leverage hii")+ylab("Cook's Distance")
  p6 <- p6 + ggtitle("Cook's dist vs Leverage hii/(1-hii)")
  p6 <- p6 + geom_abline(slope = seq(0,3,0.5), color = "gray", linetype = "dashed")
  p6 <- p6 + theme_bw() +  theme(panel.spacing = unit(0.25, 'lines'), 
      text = element_text(size = 18), 
      axis.title = element_text(size = 28),  
      legend.justification = c(0.9,0.05), 
      legend.key = element_blank() )
  
  return(list(rvfPlot = p1, qqPlot = p2, sclLocPlot = p3, cdPlot = p4, rvlevPlot = p5, cvlPlot = p6))
}

diagPlts <- diagPlot(BIC2)
g1 <- arrangeGrob(diagPlts$rvfPlot, diagPlts$qqPlot, diagPlts$sclLocPlot, diagPlts$rvlevPlot, layout_matrix = rbind(c(1,2),c(3,4)))

ggsave("./results/Diagnostic_BIC.pdf", height = 18, width = 15, units = "in", g1, device=cairo_pdf) # Fig 5.4



diagPlts<-diagPlot(LM)
g1 <- arrangeGrob(diagPlts$rvfPlot, diagPlts$qqPlot, diagPlts$sclLocPlot, diagPlts$rvlevPlot, layout_matrix = rbind(c(1,2),c(3,4)))

ggsave("./results/Diagnostic_LM.pdf", height = 18, width = 15, units = "in", g1, device=cairo_pdf) # Fig 4.6/5.3



diagPlts<-diagPlot(rsm.train)
g1 <- arrangeGrob(diagPlts$rvfPlot, diagPlts$qqPlot, diagPlts$sclLocPlot, diagPlts$rvlevPlot, layout_matrix = rbind(c(1,2),c(3,4)))

ggsave("./results/Diagnostic_RSM.pdf", height = 18, width = 15, units = "in", g1, device=cairo_pdf) # Fig 5.5



# VIFs
(vif(LM))
(vif(BIC2))
(vif(r1))
(vif(rsm.train))



# Preprocess for T = 0.9 and right ordering of factors
sub <- subset(plotdata, Tani == 0.9) # & data != "German Telephone CD 0% errors")
sub$data <- droplevels(sub$data)
sub$variant <- reorder(sub$variant, sub$Fmeasure, mean)
sub$data <- reorder(sub$data, -sub$Fmeasure, min)

sub <- sub[order(sub$Fmeasure),]

cairo_pdf("./results/Dotchart_F-Score_OptTrain.pdf", height = 18, width = 15)
dotchart(sub$Fmeasure, groups = sub$data, labels = sub$variant, cex = 2.2, 
         xlab = "Mean Prec./Rec.", main = "Mean Prec./Rec. at T = 0.9", gcolor = "#004c93", pt.cex = 3, pch = 19, lheight = 22)
dev.off()


# Print measures
print(aggregate(sub$Fmeasure, by = list(sub$data), mean))
print(aggregate(as.numeric(sub$tp), by = list(sub$data,sub$variant), mean))
print(aggregate(as.numeric(sub$fp), by = list(sub$data,sub$variant), mean))


# Table q and k by method
sub$data <- reorder(sub$data, -sub$Fmeasure, mean)
tabs <- subset(sub, !variant %in% c("Optimal Parameter choice (NN)", "Standard CLK k = 5", "Standard CLK k = 10" ))
tabs$variant <- droplevels(tabs$variant)

fileConn<-file("Parameter_Choices.tex")
writeLines(print(xtable(cbind(xtabs(formula = k~data + variant, data = tabs), xtabs(formula = q ~ data + variant, data = tabs))), booktabs = TRUE), fileConn)
close(fileConn)


# Table F-Measure + bootstrapped CIs by method
tabs <- subset(sub, !variant %in% c("Optimal Parameter choice (NN)"))
tab3 <- data.frame(data = aggregate(tabs$Fmeasure, by = list(tabs$variant), mean)$Group.1,
      meanF = aggregate(tabs$Fmeasure, by = list(tabs$variant), mean)$x,ciLow = aggregate(tabs$Fmeasure, 
      by = list(tabs$variant), function(x) (boot.mean(x, 1000)$ciLow))$x,ciHigh = aggregate(tabs$Fmeasure, 
      by = list(tabs$variant), function(x) (boot.mean(x, 1000)$ciHigh))$x)
tab3 <- tab3[order(tab3$meanF, decreasing = TRUE),]

fileConn<-file("MPR_by_Variant.tex")
writeLines(print(xtable(tab3, digits = 4), booktabs = TRUE,  include.rownames = FALSE), fileConn)
close(fileConn)


# Table meta data by data set
tabs <- subset(fulldata, q==2 & k ==1 & Tani == 0.9)
tabs <- tabs[c(2,7,4,3,1,5,6),]
posits <-  subset(plotdata, q==2 & Tani == 0.9 & variant == "Standard CLK k = 5")
posits <- posits[c(7,2,3,4,5,6,1),]
tabs <- cbind(tabs, posits$realpositives)
tabs <- tabs[,c("data", "filesizeA", "filesizeB", "posits$realpositives", "meanentropy", "meanengrams", "skew", "ginicoefficient", "errorestimation", "uniqueness", "m", "u")]
tabs <- tabs[c(7,3,4,5,2,6,1),]

fileConn<-file("Metadata_datasets.tex")
writeLines(print(xtable(tabs, digits = 3, booktabs = TRUE, display = c("s", "d",  "d", "d", "d", "f", "f","f", "f", "f", "f", "g", "g"), math.style.exponents = TRUE), include.rownames = FALSE), fileConn)
close(fileConn)

