#########################################################################
#
# Full test routine testing Optimal parameter choices against best-practice solutions
# 
#########################################################################



# Clear workspace, set working dir
#rm(list=ls())
#setwd("E:/Priv/Dropbox/Dissertation/Programme")


# Libs used
library(PPRL)
library(fastdigest)
library(stringdist)
library(entropy)
library(data.table)
library(grid)
library(gridExtra)
library(fastLink)
library(ggplot2)
library(sjstats)
library(gmodels)
library(rsm)
library(boot)
library(neuralnet)
library(randomForest)
library(MASS)
library(multibitTree)

#### File names ####
ewodataFN <- "Daten/Abgleich_UKE_SCHUFA_EMA.csv"
ncvoterAFN <- "Daten/ncvoter-20140619-temporal-balanced-ratio-1to1-a.csv"
ncvoterBFN <- "Daten/ncvoter-20140619-temporal-balanced-ratio-1to1-b.csv"
telCDAFN <- "Daten/A_10000_python_R0_C0_O100_percent.csv"
telCDBFN <- "Daten/B_10000_python_R20_C50_O100_percent.csv"
FebrlAFN <- "Daten/FEBRL_10000_20_A.csv"
FebrlBFN <- "Daten/FEBRL_10000_20_B.csv"
FebrlWAAFN <- "Daten/Adrian_Perfect_100k.csv"
FebrlWABFN <- "Daten/Adrian_10percent_100k.csv"
inputFile <- "Daten/EWMA.csv"
inputFileB <- "Daten/HOSPITAL.csv"


# Set data sets for loop
datasets <- c("German Mortality", "NC Voter", "German Telephone CD 20% errors", "German Telephone CD 0% errors", "FEBRL", "FEBRL WA", "Mortality Test Data")
variant <- c("Standard CLK k = 5", "Standard CLK k = 10", "Standard CLK k = 15", "Optimal Parameter choice (LM)", "Optimal Parameter choice (BIC)",  "50%-Rule Parameter choice", "Optimal Parameter choice (RSM)",  "Optimal Parameter choice (RF)", "Optimal Parameter choice (NN)") 
cores <- 7
leaflimit <- 3
tanimotoThresholds <- c(1.0, 0.95, 0.9, 0.85, 0.8)
fehler <-  c(0)
l <- 500
nNCV <- 50000

# For progress reports, calc number of runs needed
numCombinations <- length(tanimotoThresholds) * length(variant) * length(fehler) * 
  length(datasets) 

# Init Counter to count progress
counter <- 0


# Init result matrix
meta <- NULL


# Read training file for model building
traindata <- read.csv("Complete_Results.csv", sep="\t")

# Subset to fixed tanimoto threshold
fixed <- subset(traindata, Tani == 0.9)

# Build new complete data set with evaluation data results for rapid model comparisons
trainresultsFull <- rbind(read.csv("Trainresults_FEBRLWA.csv", sep="\t"),read.csv("Trainresults_Mortality.csv", sep="\t"), read.csv("Trainresults_NCVoter.csv", sep="\t"))

# Kill unnecessary cols
kill <- c("tp", "fp",	"fn",	"tn",	"recall",	"precision")
fulldata <- fixed[,(!names(fixed) %in% kill)]
fulldata$variant <- NULL
fulldata$positives <- NULL
fulldata$i <- 1

# Bind training to evaluation data
fulldata <- rbind(fulldata,trainresultsFull)

# Train models only on training data
fixed <- subset(traindata, Tani == 0.9)

cat("\n\n\nClassifications start!", format(Sys.time()," (%a, %d-%m-%Y, %H:%M:%S)"),"\n")



#### Function definitions ####

#### Standardization ####
unwanted_array = list('Š'='S', 'š'='s', 'Ž'='Z', 'ž'='z', 'À'='A', 'Á'='A', 'Â'='A', 'Ã'='A', 'Ä'='Ae', 'Å'='A', 'Æ'='A', 'Ç'='C', 'Ð'='D', 'È'='E', 'É'='E',
                      'Ê'='E', 'Ë'='E', 'Ì'='I', 'Í'='I', 'Î'='I', 'Ï'='I', 'Ñ'='N', 'Ö'='Oe', 'Ò'='O', 'Ó'='O', 'Ô'='O', 'Õ'='O', 'Ö'='O', 'Ø'='O', 'Ù'='U',
                      'Ú'='U', 'Û'='U', 'Ü'='Ue', 'Ý'='Y', 'Þ'='B', 'ß'='Ss', 'à'='a', 'á'='a', 'â'='a', 'ã'='a', 'ä'='ae', 'å'='a', 'æ'='a', 'ç'='c',
                      'è'='e', 'é'='e', 'ê'='e', 'ë'='e', 'ì'='i', 'í'='i', 'î'='i', 'ï'='i', 'ð'='o', 'ñ'='n', 'ò'='o', 'ó'='o', 'ô'='o', 'õ'='o',
                      'ö'='oe', 'ø'='o', 'ü'='ue', 'ù'='u', 'ú'='u', 'û'='u', 'ý'='y', 'ý'='y', 'þ'='b', 'ÿ'='y' )

normalize_string=function(string){
  string = enc2utf8(string)
  string=gsub("[]$*+.?[^{|(\\#%&~_/<=>'!,:;`\")}@]", "",string)
  string=gsub("-", "",string)
  string =chartr(paste(names(unwanted_array), collapse=''),
                 paste(unwanted_array, collapse=''),
                 string)
  string=gsub("[\t\n\r\f\v]", "",string)
  string=gsub("[0-9]", "",string)
  string=toupper(string)
  string=iconv(string, "utf8", "ASCII", sub="")
  return(string)
}

normalize_dobs=function(string){
  if(is.na(string)){
    return("")
  } else {
    string = as.character(string)
    if (nchar(string) == 1){
      string = paste0("0",string)
    } else if (nchar(string) == 0){
      string = ""
    }
    return(string)
  }
}


## Bootstrap mean 
boot.mean = function(x,B,binwidth=NULL){
  n = length(x)
  boot.samples = matrix( sample(x,size=n*B,replace=TRUE), B, n)
  boot.statistics = apply(boot.samples,1,mean)
  se = sd(boot.statistics)
  require(ggplot2)
  if ( is.null(binwidth) )
    binwidth = diff(range(boot.statistics))/30
  p = ggplot(data.frame(x=boot.statistics),aes(x=x)) +geom_histogram(aes(y=..density..),binwidth=binwidth) + geom_density(color="red")
  #plot(p)
  interval = mean(x) + c(-1,1)*2*se
  print( interval )
  return( list(boot.statistics = boot.statistics, interval=interval, se=se, plot=p) )
}




# Fixed model calls
# BIC

# Full model for BIC-selection
fullmod <- lm(Fmeasure ~  k + (q) + errorestimation +  skew + uniqueness + meanentropy + log(u) + ginicoefficient +  uniquepatternsA + errorestimation +  k * q *  hammingweight +  meanengrams*hammingweight + 
                meanmissing +  log(m) + u*m, 
              data=fixed)
# Minimal model
minmod <- lm(Fmeasure ~ k + (q), data = fixed)

# Use BIC-selection to get the model
stepTest <- stepAIC(fullmod, minmod, direction = "both",  k = log(nrow(fixed))) # bic
#stepTest2 <- stepAIC(fullmod, minmod, direction = "both",  k = 2) # aic
BIC2 <- stepTest


# LM
LM <- lm(Fmeasure ~ k + (q) + errorestimation * skew + uniquepatternsA * 
           uniqueness + ginicoefficient + k *  hammingweight + 
           hammingweight * meanengrams + uniquepatternsA + log(u),
         data=fixed)


# Tune mtry value
tune <- subset(traindata, Tani == 0.9)
tuneRF(tune[,c("k", "q", "errorestimation", "skew", "uniquepatternsA", "uniqueness", "ginicoefficient", "hammingweight", "meanengrams", "u", "m")], tune[,"Fmeasure"], ntreeTry=3000, stepFactor=1, plot=FALSE)

# Random forest
r1 <- randomForest(Fmeasure ~ k  +  q + errorestimation +  skew + uniqueness + ginicoefficient + hammingweight + 
                     meanengrams + u + uniquepatternsA, data=fixed, 
                   mtry= 3, ntree=500, nodesize=20, #nodesize=5,
                   importance=TRUE, nPerm=5, keep.forest=TRUE)


# RSM: Preproc data
RSMdata <- subset(fixed, Tani == 0.9)[,c("k", "q", "u", "meanengrams", "hammingweight", "Fmeasure")]
RSMdata_unscaled <- subset(fixed, Tani == 0.9)[,c("k", "q", "u", "meanengrams", "hammingweight", "Fmeasure")]
# Scale min/max
maxi = apply(RSMdata , 2, function(x) max(as.numeric(x)))
mini = apply(RSMdata,  2, function(x) min(as.numeric(x)))
# Scaling
RSMdata = as.data.frame(scale(RSMdata, center = mini, scale = maxi - mini))
# Train model
rsm.train <- rsm(Fmeasure ~   SO(k, q) + TWI(q, k, hammingweight) + (u), data = RSMdata)


# Neural network
NNata <- fixed[fixed$Tani == 0.9 ,c("k", "q", "skew", "uniqueness", "errorestimation", "hammingweight", "ginicoefficient", "meanengrams", "uniquepatternsA", "u", "Fmeasure")]
NNata_unscaled <- NNata
# Scale min/max
maxi = apply(NNata , 2, function(x) max(as.numeric(x)))
mini = apply(NNata,  2, function(x) min(as.numeric(x)))
# Scaling
NNata = as.data.frame(scale(NNata, center = mini, scale = maxi - mini))
# fit neural network
set.seed(42)
NN <- neuralnet(Fmeasure ~  k + q +  ginicoefficient +  hammingweight + 
                   uniquepatternsA + (u), NNata, hidden = 2, linear.output = T, rep = 1)


#### Loop over datasets ####
for (d in datasets){

  
  # RSM: Preproc data
  newdataRSM <- subset(fulldata, Tani == 0.9 & data == d)
  newdataRSM <- newdataRSM[,c("k", "q", "u", "meanengrams", "hammingweight", "Fmeasure")]
  
  maxi = apply(newdataRSM , 2, function(x) max(as.numeric(x)))
  mini = apply(newdataRSM,  2, function(x) min(as.numeric(x)))
  
  # Scale data
  newdataRSM = as.data.frame(scale(newdataRSM, center = mini, scale = maxi - mini))
  newdataRSM$u <- 0.5
  newdataRSM$meanengrams <- 0.5
  
  # Predict RSM
  RSM <- data.frame(newdataRSM, predicted = predict(rsm.train, newdataRSM, type = "response"))
  # Re-Scale to get actual values
  RSM$Fmeasure <- (RSM$Fmeasure  * (max(RSMdata_unscaled$Fmeasure) - min(RSMdata_unscaled$Fmeasure))) + min(RSMdata_unscaled$Fmeasure)
  RSM$k <- (RSM$k  * (max(RSMdata_unscaled$k) - min(RSMdata_unscaled$k))) + min(RSMdata_unscaled$k)
  RSM$q <- (RSM$q  * (max(RSMdata_unscaled$q) - min(RSMdata_unscaled$q))) + min(RSMdata_unscaled$q)
  
  
  # New data to predict from
  newdata <- subset(fulldata, Tani == 0.9 & data == d)# & q == 2)
  
  # Predict LM/BIC/RF
  newfitted <- data.frame(newdata, predicted = predict(LM, newdata, type = "response"))
  newfittedBIC <- data.frame(newdata, predicted = predict(BIC2, newdata, type = "response"))
  newfittedRF <- data.frame(newdata, predicted = predict(r1, newdata))
  
  # Get optimal prameter estimate
  parameters <- newfitted[which.max((newfitted$predicted)), c("Fmeasure", "Tani", "q", "k", "l", "predicted")]
  parameters2 <- newfittedBIC[which.max((newfittedBIC$predicted)), c("Fmeasure", "Tani", "q", "k", "l", "predicted")]
  parameters3 <- newfittedRF[which.max((newfittedRF$predicted)), c("Fmeasure", "Tani", "q", "k", "l", "predicted")]
  parametersRSM <- RSM[which.max((RSM$predicted)), c("Fmeasure", "q", "k", "predicted")]
  
  
  ## Prediction using neural network
  newdataNN <- fulldata[fulldata$Tani == 0.9 & fulldata$data == d ,colnames(NN$covariate)]
  newdataNN_unscaled <- newdataNN
  # Scale again
  maxi = apply(newdataNN , 2, function(x) max(as.numeric(x)))
  mini = apply(newdataNN,  2, function(x) min(as.numeric(x)))
  # Scaling
  newdataNN = as.data.frame(scale(newdataNN, center = mini, scale = maxi - mini))
  # Unscalable because constant values -> 0.5
  newdataNN[sapply(newdataNN, is.nan)] <- 0.5
  
  # Predict and rescale to original values
  predict_testNN = neuralnet::compute(NN, newdataNN)
  predict_testNN = (predict_testNN$net.result * (max(NNata_unscaled$Fmeasure) - min(NNata_unscaled$Fmeasure))) + min(NNata_unscaled$Fmeasure)
  # Prediction data frame
  newfitNN <- cbind(fulldata[fulldata$Tani == 0.9 & fulldata$data == d,], predict_testNN)
  
  # Get NN parameters
  kNN <- round(mean(newfitNN[which(newfitNN$predict_testNN %in% newfitNN$predict_testNN[newfitNN$predict_testNN >= (quantile(as.numeric((newfitNN$predict_testNN)), probs=seq(0,1,0.1))["90%"])]),"k"]))
  qNN <- round(mean(newfitNN[which(newfitNN$predict_testNN %in% newfitNN$predict_testNN[newfitNN$predict_testNN >= (quantile(as.numeric((newfitNN$predict_testNN)), probs=seq(0,1,0.1))["90%"])]),"q"]))
  
  # Get RSM parameters
  kRSM <-  parametersRSM$k
    #round(mean(RSM[which(RSM$predicted %in% RSM$predicted[RSM$predicted >= (quantile(as.numeric((RSM$predicted)), probs=seq(0,1,0.1))["90%"])]),"k"]))
  qRSM <- parametersRSM$q
    #round(mean(RSM[which(RSM$predicted %in% RSM$predicted[RSM$predicted >= (quantile(as.numeric((RSM$predicted)), probs=seq(0,1,0.1))["90%"])]),"q"]))
  
  
  # Write out parameters for use in loop
  write.table(parameters, "Optimal_Parameters.csv",sep="\t", row.names = FALSE)
  write.table(parameters2, "Optimal_Parameters_BIC.csv",sep="\t", row.names = FALSE)
  write.table(parameters3, "Optimal_Parameters_RF.csv",sep="\t", row.names = FALSE)
  
  
  #### Preprocess depending on data ####
  if (d =="NC Voter"){
    
    # Read
    clearTextB <- fread(ncvoterBFN, sep=",")
    clearTextA <- fread(ncvoterAFN, sep=",")
    
    # Subset identifiers
    clearTextA <- clearTextA[,c(1,2,4,5)]
    clearTextB <- clearTextB[,c(1,2,4,5)]
    
    # Subset overlapping records
    clearTextA <- clearTextA[(clearTextA$voter_id %in% clearTextB$voter_id),]
    clearTextB <- clearTextB[(clearTextB$voter_id %in% clearTextA$voter_id),]
    
    # Deduplicate
    clearTextA <- clearTextA[!duplicated(clearTextA$voter_id),]
    clearTextB <- clearTextB[!duplicated(clearTextB$voter_id),]
    
    # Draw sample
    set.seed(21)
    clearTextA <- clearTextA[sample(1:nrow(clearTextA), nNCV, replace=FALSE),]
    set.seed(21)
    clearTextB <- clearTextB[sample(1:nrow(clearTextB), nNCV, replace=FALSE),]
    
    clearTextA$Day <- ""
    clearTextA$Month <- ""
    clearTextA$Year <- clearTextA$age
    
    clearTextB$Day <- ""
    clearTextB$Month <- ""
    clearTextB$Year <- clearTextB$age
    clearTextA$age <- NULL
    clearTextB$age <- NULL
    
    ### names anpassen für alle datasets
    
  } else if (d == "German Telephone CD 20% errors"){
    
    clearTextB <- fread(telCDBFN, sep="\t")
    clearTextA <- fread(telCDAFN, sep="\t")
    
    # Reorder
    clearTextA <- clearTextA[,c(1,7,8,5,4,3)]
    clearTextB <- clearTextB[,c(1,7,8,5,4,3)]
    
    
  } else if (d == "German Telephone CD 0% errors"){
    
    clearTextB <- fread(telCDAFN, sep="\t")
    clearTextA <- fread(telCDAFN, sep="\t")
    
    # Reorder
    clearTextA <- clearTextA[,c(1,7,8,5,4,3)]
    clearTextB <- clearTextB[,c(1,7,8,5,4,3)]
    
    
  } else if (d == "German Mortality"){
    clearText <- fread(ewodataFN, sep=",")
    
    clearTextB <- subset(clearText, quelle=="uke")
    clearTextA <- subset(clearText, quelle=="einwo")
    clearTextA <- clearTextA[,c(14,2,3,5,6,7)]
    clearTextB <- clearTextB[,c(14,2,3,5,6,7)]
    
  } else if (d == "Mortality Test Data"){
    clearText <- fread(ewodataFN, sep=",")
    
    clearTextB <- subset(clearText, quelle=="schufa")
    clearTextA <- subset(clearText, quelle=="einwo")
    clearTextA <- clearTextA[,c(14,2,3,5,6,7)]
    clearTextB <- clearTextB[,c(14,2,3,5,6,7)]
    
  } else if (d == "FEBRL"){
    
    clearTextB <- fread(FebrlBFN, sep="\t")
    clearTextA <- fread(FebrlAFN, sep="\t")
    
    clearTextA <- clearTextA[,c(1,2,3,10)] 
    clearTextB <- clearTextB[,c(1,2,3,10)] 
    clearTextA$Day <- (substr(clearTextA$V10,7,8))
    clearTextA$Month <- (substr(clearTextA$V10,5,6))
    clearTextA$Year <- (substr(clearTextA$V10,1,4))
    clearTextA$V10 <- NULL
    
    clearTextB$Day <- (substr(clearTextB$V10,7,8))
    clearTextB$Month <- (substr(clearTextB$V10,5,6))
    clearTextB$Year <- (substr(clearTextB$V10,1,4))
    clearTextB$V10 <- NULL
    
    
  } else {
    clearTextA <- fread(FebrlWAAFN, colClasses = "character", stringsAsFactors = FALSE)
    clearTextB <- fread(FebrlWABFN, colClasses = "character", stringsAsFactors = FALSE)
    
    clearTextA$Day <- (substr(clearTextA$`Date Of Birth`,7,8))
    clearTextA$Month <- (substr(clearTextA$`Date Of Birth`,5,6))
    clearTextA$Year <- (substr(clearTextA$`Date Of Birth`,1,4))
    clearTextA$`Date Of Birth` <- NULL
    clearTextB$Day <- (substr(clearTextB$`Date Of Birth`,7,8))
    clearTextB$Month <- (substr(clearTextB$`Date Of Birth`,5,6))
    clearTextB$Year <- (substr(clearTextB$`Date Of Birth`,1,4))
    clearTextB$`Date Of Birth` <- NULL
    
    
    clearTextA <- clearTextA[,c(2,4,6,12,13,14)] 
    clearTextB <- clearTextB[,c(2,4,6,12,13,14)] 
    
    
    
    
  }
  
  # Read error-free file
  names(clearTextA) <- c("ID", "Vorname", "Nachname", "Day", "Month", "Year")
  names(clearTextB) <- c("ID", "Vorname", "Nachname", "Day", "Month", "Year")
  
  
  
  #### General preproc ####
  clearTextA$ID <-  as.character(clearTextA$ID)
  clearTextB$ID <-  as.character(clearTextB$ID)
  
  
  # Preprocess A
  
  ## Standardize data
  clearTextA$Vorname <- normalize_string(clearTextA$Vorname)
  clearTextA$Nachname <- normalize_string(clearTextA$Nachname)
  
  clearTextA$Day <- sapply(clearTextA$Day, normalize_dobs)
  clearTextA$Month <- sapply(clearTextA$Month, normalize_dobs)
  clearTextA$Year <- as.character(ifelse(is.na(clearTextA$Year), "", clearTextA$Year))
  
  # Generate full DOB
  clearTextA$DOB <-  paste(clearTextA$Day,clearTextA$Month,clearTextA$Year, sep="")
  #clearTextA$DOB <- substr(clearTextA$DOB, 3,6)
  
  
  # Preprocess B
  
  ## Standardize data
  clearTextB$Vorname <- normalize_string(clearTextB$Vorname)
  clearTextB$Nachname <- normalize_string(clearTextB$Nachname)
  
  clearTextB$Day <- sapply(clearTextB$Day, normalize_dobs)
  clearTextB$Month <- sapply(clearTextB$Month, normalize_dobs)
  clearTextB$Year <- as.character(ifelse(is.na(clearTextB$Year), "", clearTextB$Year))
  
  # Generate full DOB for salting
  clearTextB$DOB <-  paste(clearTextB$Day,clearTextB$Month,clearTextB$Year, sep="")
  #clearTextB$DOB <- substr(clearTextB$DOB, 3,6)
  
  ## Deduplicate by ID
  clearTextA <- clearTextA[!duplicated(clearTextA$ID),]
  clearTextB <- clearTextB[!duplicated(clearTextB$ID),]
  
  
  ## File sizes, overlap
  filesizeA <- (nrow(clearTextA))
  filesizeB <- (nrow(clearTextB))
  
  pairs <- as.numeric(filesizeA) * as.numeric(filesizeB)
  truepositives <- sum(clearTextA$ID %in% clearTextB$ID)
  positives <- truepositives
  
  clearTextA$LinkKey <- paste(clearTextA$Vorname, clearTextA$Nachname, clearTextA$DOB, sep="")
  clearTextB$LinkKey <- paste(clearTextB$Vorname, clearTextB$Nachname, clearTextB$DOB, sep="")
  
  cat("\n\n\nClear-Text files ready, analyzing files", format(Sys.time()," (%a, %d-%m-%Y, %H:%M:%S)"),"\n")
  
  
  
  # Loop over every combination of parameters
  for (v in variant){
    
    ### Encryption phase ###
    cat("\n\n\nEncryption start!", format(Sys.time()," (%a, %d-%m-%Y, %H:%M:%S)"),"\n")
    
    if (v == "Standard CLK k = 20"){
      
      # Parameters
      k <- 20
      q <- 2
      RMSE <- NA
      
    } else   if (v == "Standard CLK k = 10"){
      
      # Parameters
      k <- 10
      q <- 2
      RMSE <- NA
      
    } else   if (v == "Standard CLK k = 5"){
      
      # Parameters
      k <- 5
      q <- 2
      RMSE <- NA 
      
    } else   if (v == "Standard CLK k = 15"){
      
      # Parameters
      k <- 15
      q <- 2
      
      RMSE <- NA
      
    } else   if (v == "Optimal Parameter choice (LM)"){
      
      # Calculate Root Mean Square Error (RMSE)
      RMSE = (sum((newfitted$Fmeasure - newfitted$predicted)^2) / nrow(newfitted)) ^ 0.5
      
      # Parameters
      k <- read.csv("Optimal_Parameters.csv", sep="\t")$k
      q <- read.csv("Optimal_Parameters.csv", sep="\t")$q
      
    } else   if (v == "Optimal Parameter choice (RF)"){
      
      # Calculate Root Mean Square Error (RMSE)
      RMSE = (sum((newfittedRF$Fmeasure - newfittedRF$predicted)^2) / nrow(newfittedRF)) ^ 0.5
      
      # Parameters
      k <- read.csv("Optimal_Parameters_RF.csv", sep="\t")$k
      q <- read.csv("Optimal_Parameters_RF.csv", sep="\t")$q    
      
      
    } else   if (v == "Optimal Parameter choice (BIC)"){
      
      # Parameters
      k <- read.csv("Optimal_Parameters_BIC.csv", sep="\t")$k
      q <- read.csv("Optimal_Parameters_BIC.csv", sep="\t")$q
      
      # Calculate Root Mean Square Error (RMSE)
      RMSE = (sum((newfittedBIC$Fmeasure - newfittedBIC$predicted)^2) / nrow(newfittedBIC)) ^ 0.5
      
    } else   if (v == "Optimal Parameter choice (RSM)"){
      
      # Parameters
      k <-  kRSM
      q <-  qRSM
      
      # Calculate Root Mean Square Error (RMSE)
      RMSE = (sum((RSM$Fmeasure - RSM$predicted)^2) / nrow(RSM)) ^ 0.5
      
      
    } else   if (v == "Optimal Parameter choice (NN)"){
      
      # Parameters
      k <-  kNN
      q <-  qNN
      
      # Calculate Root Mean Square Error (RMSE)
      RMSE = (sum((newfitNN$Fmeasure -  newfitNN$predict_testNN)^2) / nrow(newfitNN)) ^ 0.5
      
    } else   if (v == "50%-Rule Parameter choice"){
      
      # Parameters
      q <- 2
      kA <- round((l / round(mean(nchar(clearTextA$LinkKey)-(q-1))))  * log(2))
      kB <- round((l / round(mean(nchar(clearTextB$LinkKey)-(q-1))))  * log(2))
      k <- max(c(kA,kB))
      
    } else {
      stop("Illegal variant: ", v)
    }
    
    # Encrypt to CLKs
    encryptedA <- CreateCLK(ID = clearTextA$ID, clearTextA[, c(2,3,7)], k = k, padding = rep(0,3), q = c(q,q,q), l = l, password = c("12341231", "41412412", "415212512"))
    encryptedB <- CreateCLK(ID = clearTextB$ID, clearTextB[, c(2,3,7)], k = k, padding = rep(0,3), q = c(q,q,q), l = l, password = c("12341231", "41412412", "415212512"))
    

    
    cat("\n\n\nEncryption finished!", format(Sys.time()," (%a, %d-%m-%Y, %H:%M:%S)"),"\n")
    # Encrypt
    
    # Write Tree data
    write.table(encryptedA, "TreeA.csv", sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
    write.table(encryptedB, "TreeB.csv", sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
    
    
    ## Loop over thresholds
    for (Tani in tanimotoThresholds){
      
      # Run MBT
      multibitTree.load("TreeB.csv",threads=cores, leafLimit = leaflimit)
      result <- multibitTree.searchFile("TreeA.csv", Tani)
      
      
      # Classify candidate pairs
      tp <- sum((result$fingerprint) == (result$query))
      fp <- sum((result$fingerprint) != (result$query))
      fn <- positives - tp
      tn <- pairs - (fn + tp + fp)
      
      # Calculate central measures
      ReductionRate <-  1 - ((tp+fp)/(tn+fn))
      precision <- tp / (tp+fp)
      recall <- tp / (tp+fn)
      Fmeasure =   (recall + precision) / 2
      
      rm(result)
      #gc()
      
      # Set counter plus one
      counter <- counter + 1
      
      cat("\nClassification", counter, "of", numCombinations,"Finished (", format((counter / numCombinations)*100, digits=3),"%)!", format(Sys.time(), " (%a, %d-%m-%Y, %H:%M:%S)"), ": \nFilesizes:\t", format(filesizeA,scientific=FALSE), "\nData source:\t", d, "\nMethod:\t\t", v, "\nk/q:\t\t", k, "/", q, "\nTanimoto:\t", Tani, "\n\nPositives:\t", truepositives, "\nTrue Positives:\t", tp, "\nRecall:\t\t", recall, "\nFalse Positives:", fp, "\nPrecision:\t", precision, "\nFalse Negatives:", fn, "\nF-Score:\t", Fmeasure, "\n\n")
      
      # Append to result
      meta  <- rbind(meta, cbind(d, v, q, k, l, Tani, filesizeA,filesizeB,positives, tp, fp, fn, tn, recall, precision,Fmeasure, RMSE))
      write.table(meta, "tempResult.csv", sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)
      
    }
  }
}


# Format Output column names and write final result file
colnames(meta) <- c("data", "variant", "q", "k", "l", "Tani", "filesizea", "filesizeb", "realpositives", "tp", "fp", "tn", "fn", "recall", "precision", "Fmeasure", "RMSE")
write.table(meta,"./results/R_Metainfo_Testbett_CompareOPC_vs_Train.csv", sep="\t", row.names=FALSE, col.names=TRUE)

cat("\n\n\nAll Classifications Finished!", format(Sys.time()," (%a, %d-%m-%Y, %H:%M:%S)"),"\n")

