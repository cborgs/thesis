

# Clear workspace, set working dir
#rm(list=ls())
#setwd("E:/Priv/Dropbox/Dissertation/Programme")


# Libs used
library(PPRL)
library(fastdigest)
library(stringdist)
library(entropy)
library(data.table)
library(fastLink)
library(ggplot2)



#### File names ####
ewodataFN <- "Daten/Abgleich_UKE_SCHUFA_EMA.csv"
ncvoterAFN <- "Daten/ncvoter-20140619-temporal-balanced-ratio-1to1-a.csv"
ncvoterBFN <- "Daten/ncvoter-20140619-temporal-balanced-ratio-1to1-b.csv"
telCDAFN <- "Daten/A_10000_python_R0_C0_O100_percent.csv"
telCDBFN <- "Daten/B_10000_python_R20_C50_O100_percent.csv"
FebrlAFN <- "Daten/FEBRL_10000_20_A.csv"
FebrlBFN <- "Daten/FEBRL_10000_20_B.csv"

#### Loop vars ####
# Set data sets for loop
datasets <- c("German Mortality", "German Telephone CD 20% errors", "German Telephone CD 0% errors", "FEBRL")  
# Set to standard identifier set for CLK
v <- "Standard CLK"
# Amount of hash functions to loop over
K <- seq(from = 1, to = 40, by = 1) # Amount of Hash functions
# q-grams to loop over
Q <- c(1,2,3,4)
# Length of the BF
l <- 500 
# Size of NC Voter data subsample
nNCV <- 50000
# No loopable errors in data
fehler <-  c(0)
# Set padding to false
padding <- FALSE
# Thresholds to loop over
tanimotoThresholds <- seq(from= 1.0, to = 0.8, by = -0.01)
# Set cores and leaf limit
cores <- 7
leaflimit <- 3


# Init Counter to count progress
counter <- 0

# Calc number of combinations in total
numCombinations <- length(tanimotoThresholds) * length(fehler) * 
  length(datasets) * length(K) * length(Q)


# Init result matrix
meta <- NULL

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
    if (nchar(string) < 2){
      string = paste0("0",string)
    }
    return(string)
  }
}




#### Loop over datasets ####

for (d in datasets){
  
 
  #### Preprocess depending on data ####  
  if (d =="NC Voter"){
    
    # Read
    clearTextB <- fread(ncvoterBFN, sep=",")
    clearTextA <- fread(ncvoterAFN, sep=",")
    
    # Subset identifiers
    clearTextA <- clearTextA[,c(1,2,4,5)] 
    clearTextB <- clearTextB[,c(1,2,4,5)] 
    
    # Draw sample
    set.seed(69)
    clearTextA <- clearTextA[sample(1:nrow(clearTextA), nNCV, replace=FALSE),]
    set.seed(69)
    clearTextB <- clearTextB[sample(1:nrow(clearTextB), nNCV, replace=FALSE),]
    
    clearTextA$Day <- ""
    clearTextA$Month <- ""
    clearTextA$Year <- as.character(2014-clearTextA$age)
    
    clearTextB$Day <- ""
    clearTextB$Month <- ""
    clearTextB$Year <- as.character(2014-clearTextB$age)
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
    
    clearTextB <- subset(clearText, quelle=="schufa")
    clearTextA <- subset(clearText, quelle=="einwo")
    clearTextA <- clearTextA[,c(14,2,3,5,6,7)] 
    clearTextB <- clearTextB[,c(14,2,3,5,6,7)] 
    
  } else {
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

  
  # Preprocess B
  
  ## Standardize data
  clearTextB$Vorname <- normalize_string(clearTextB$Vorname)
  clearTextB$Nachname <- normalize_string(clearTextB$Nachname)
  
  clearTextB$Day <- sapply(clearTextB$Day, normalize_dobs)
  clearTextB$Month <- sapply(clearTextB$Month, normalize_dobs)
  clearTextB$Year <- as.character(ifelse(is.na(clearTextB$Year), "", clearTextB$Year))
  
  # Generate full DOB for salting
  clearTextB$DOB <-  paste(clearTextB$Day,clearTextB$Month,clearTextB$Year, sep="")
  
  
  ## Deduplicate by ID
  clearTextA <- clearTextA[!duplicated(clearTextA$ID),]
  clearTextB <- clearTextB[!duplicated(clearTextB$ID),]
  

  ## File sizes, overlap
  filesizeA <- (nrow(clearTextA))
  filesizeB <- (nrow(clearTextB))
  
  pairs <- as.numeric(filesizeA) * as.numeric(filesizeB)
  truepositives <- sum(clearTextA$ID %in% clearTextB$ID)
  positives <- truepositives
  
  
  cat("\n\n\nClear-Text files ready, analyzing files", format(Sys.time()," (%a, %d-%m-%Y, %H:%M:%S)"),"\n")
  
  ## Loop over combos ##
  for (q in Q){
    for (k in K){
    
      # Placeholder measures
      ngrams <- q * k
      meanentropy <- 1
      meanengrams <- 1
      entropyIdentifier <- 1
      skewness <- 1
      ginicoefficient <- 1
      errorestimation  <- 1 
      missingamount <- 1
      m <- 1
      u <- 1
      
      cat("\n\n\nEncrypting...", format(Sys.time()," (%a, %d-%m-%Y, %H:%M:%S)"),"\n")
      
      # Encrypt
      encryptedA <- CreateCLK(ID = clearTextA$ID, clearTextA[, c(2,3,7)], k = k, padding = rep(0,3), q = c(q,q,q), l = l, password = c("(H]$6Uh*-Z204q", "asd", "afsd"))
      encryptedB <- CreateCLK(ID = clearTextB$ID, clearTextB[, c(2,3,7)], k = k, padding = rep(0,3), q = c(q,q,q), l = l, password = c("(H]$6Uh*-Z204q", "asd", "afsd"))
      
      # Calculate hamming weight and no of unique patterns
      hammingweight <- mean(c(nchar(gsub("0", "", encryptedA$CLKs)),nchar(gsub("0", "", encryptedB$CLKs))))
      uniquepatternsA <- length(unique(encryptedA$CLKs))
      uniquepatternsB <- length(unique(encryptedB$CLKs))
      
      # Tree data
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
  
        # Set counter plus one
        counter <- counter + 1
        
        cat("\nClassification", counter, "of", numCombinations,"Finished (", format((counter / numCombinations)*100, digits=3),"%)!", format(Sys.time()," (%a, %d-%m-%Y, %H:%M:%S)"),": \nFilesizes:\t", format(filesizeA,scientific=FALSE), "\nData source:\t", d, "\nMethod:\t\t", v, "\nk/q:\t\t", k, "/", q, "\nTanimoto:\t", Tani, "\n\nPositives:\t", truepositives, "\nTrue Positives:\t", tp,"\nRecall:\t\t", recall, "\nFalse Positives:", fp, "\nPrecision:\t", precision, "\nFalse Negatives:", fn, "\nF-Score:\t", Fmeasure, "\n\n")
        
        # Append to result
        meta  <- rbind(meta, cbind(d, v, q, k, l, Tani, filesizeA,filesizeB,positives, tp, fp, fn, tn, recall, precision,Fmeasure, ngrams, meanentropy, meanengrams, entropyIdentifier, skewness, ginicoefficient, errorestimation , missingamount,hammingweight,uniquepatternsA,uniquepatternsB,m,u))
        write.table(meta, "tempResult.csv", sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)
        
      }
    }
  }
}


# format result
result <- as.data.frame(meta)
colnames(result) <- c("data", "variant", "q", "k", "l", "Tani", "filesizeA", "filesizeB", "positives", "tp", "fp", "fn", "tn", "recall", "precision", "Fmeasure", "ngrams", "meanentropy", "meanengrams", "entropyIdentifier", "skewness", "ginicoefficient", "errorestimation", "missingamount", "hammingweight", "uniquepatternsA", "uniquepatternsB", "m", "u")

# Write final result file
write.table(result, "fullResult.csv", sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)

rm(list=ls())
gc()


#### Merging metadata to result data ####

# Read result file
result <- read.csv("fullResult.csv", sep="\t")


result2 <- read.csv(paste0("fullResult_Metadata.csv"), sep="\t")
names(result2)[1] <- "data"

result <- result[,!names(result) %in% names(result2)[7:22]]
result$entropyIdentifier <- NULL
result$skewness <- NULL
result$ngrams <- NULL
result2$l <- NULL
result2$filesizeA <- NULL
result2$filesizeB <- NULL
result2$positives <- NULL
result <- subset(result, Tani >= 0.7)


# Merge metadata with simulation results
plotdata <- merge(result, result2, by=c("data", "q"), all.x=TRUE)
rm(result, result2)

write.table(plotdata, "Complete_Results.csv", sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)

