
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
library(moments)
library(ineq)

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
datasets <- c("German Mortality","German Telephone CD 20% errors","German Telephone CD 0% errors", "FEBRL")

# Set to standard identifier set for CLK
v <- "Standard CLK"
# Amount of hash functions to loop over
K <- seq(from = 1, to = 40, by = 1) # Amount of Hash functions
# q-grams to loop over
Q <- c(1,2,3,4)
# Length of the BF
l <- 500 
# No loopable errors in data
fehler <-  c(0)
# Set padding to false
padding <- FALSE
# Size of NC Voter data subsample
nNCV <- 50000
# Thresholds to loop over
tanimotoThresholds <- seq(from= 1.0, to = 0.8, by = -0.01)
# Set cores and leaf limit
cores <- 7
leaflimit <- 3



# identifiers <- c("FN","LN", "SEX", "YEAR", "DAY", "MONTH","PLACE")
# 
# # Build Identifier Sets from all combinations (except for nonsense two-ID cases, only FNLN FNYEAR LNYEAR count)
# identifierSets <- c(paste(identifiers[c(1,2)],  collapse=""), paste(identifiers[c(1,4)],  collapse=""), 
#                     paste(identifiers[c(2,4)],  collapse=""))
# identifierSets <- c(identifierSets, apply(t(combn(identifiers,3)), 1, paste, collapse=""))
# identifierSets <- c(identifierSets, apply(t(combn(identifiers,4)), 1, paste, collapse=""))
# identifierSets <- c(identifierSets, apply(t(combn(identifiers,5)), 1, paste, collapse=""))
# identifierSets <- c(identifierSets, apply(t(combn(identifiers,6)), 1, paste, collapse=""))
# identifierSets <- c(identifierSets, paste(identifiers,  collapse=""))

# Init Counter to count progress
counter <- 0

# Calc number of combinations in total
numCombinations <-  length(datasets) * length(Q) 


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
  string=gsub("[]$*+.?[^{|(\\#%&~_/<=>'!,:;`\")}@]","",string)
  string=gsub("-","",string)
  string =chartr(paste(names(unwanted_array), collapse=''),
                 paste(unwanted_array, collapse=''),
                 string)
  string=gsub("[\t\n\r\f\v]","",string)
  string=gsub("[0-9]","",string)
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
    
    clearTextB <- fread(ncvoterBFN, sep=",")
    clearTextA <- fread(ncvoterAFN, sep=",")
    
    clearTextA <- clearTextA[,c(1,2,4,5)] 
    clearTextB <- clearTextB[,c(1,2,4,5)] 
    
    # Draw sample
    set.seed(69)
    clearTextA <- clearTextA[sample(1:nrow(clearTextA), nNCV, replace=FALSE),]
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
  names(clearTextA) <- c("ID","Vorname","Nachname","Day","Month","Year")
  names(clearTextB) <- c("ID","Vorname","Nachname","Day","Month","Year")
  
  
  
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
  
  
  
  ###############
  
  # Get EM weights
  ests <- fastLink(clearTextA, clearTextB, varnames =  c("Vorname","Nachname","DOB"), 
                   stringdist.match=c("Vorname", "Nachname"),  stringdist.method="jw",  n.cores = 7, estimate.only = TRUE)
  
  m <- ests$p.m
  u <- ests$p.u
  
  
  for (q in Q){
    
    cat("\n\n\nLine-by-line checking...", format(Sys.time()," (%a, %d-%m-%Y, %H:%M:%S)"),"\n")
    
    
    # Generate full linkage key
    keys <- apply(clearTextA, 1, paste0, collapse ="")
    keys <- gsub(" ","", keys)
    
    # Split into bigrams
    qgrams <- (sapply(keys,function(key)substring(key,first=seq(1,nchar(key)),last=seq(q, nchar(key)+q))))
    # Count uniques and save measures
    uniquengrams <- as.numeric(unlist(lapply(qgrams, function(x)length(unique(x)))))
    # ngram Count
    ngramsC <- as.numeric(unlist(lapply(qgrams, function(x)length(x))))
    qgrams <- unlist(qgrams)
    qgramlist <- as.character(na.omit(ifelse(qgrams == "" | nchar(qgrams)!= q,NA,qgrams)))
    
    cat("\n\n\nLine-by-line checking complete, calculating measures...", format(Sys.time()," (%a, %d-%m-%Y, %H:%M:%S)"),"\n")
    
    # Sum missings
    missings <- sum(apply(clearTextA, 2, function(x)sum(is.na(x)))) + sum(apply(clearTextA, 2, function(x)sum(x=="")))
    
    # Mean entropy calculation
    meanentropy <- mean(apply(clearTextA, 2, function(x) entropy(table(x))))
    
    # q90/q10-ratio
    q90 <- round(quantile(as.numeric(table(qgramlist)), probs=seq(0,1,0.1))["90%"])
    q10 <- round(quantile(as.numeric(table(qgramlist)), probs=seq(0,1,0.1))["10%"])
    
    # Measures
    meanengrams <- mean(ngramsC)
    sdgrams <- sd(ngramsC)
    
    skew <-   skewness(table(qgramlist))
    ginicoefficientqgrams <- Gini(table(qgramlist))
    
    ginicoefficient <- mean(apply(clearTextA, 2, function(x) Gini(table(x))))
    
    errorestimation  <- q90/q10 # amount of rare combinations in data
    missingamount <- missings
    meanmissing <- missings/nrow(clearTextA)
    uniqueness <- mean(uniquengrams/ngramsC)
    sduniqueness <- sd(uniquengrams/ngramsC) 
  
    
    
    # Set counter plus one
    counter <- counter + 1
    
    #gc()
    cat("\nClassification", counter, "of", numCombinations,"Finished (",format((counter / numCombinations)*100, digits=3),"%)!",format(Sys.time()," (%a, %d-%m-%Y, %H:%M:%S)"),": \nFilesizes:\t",format(filesizeA,scientific=FALSE),"\nData source:\t",d,"\n q:\t\t",q,"\n\n")
    
    # Append to result
    meta  <- rbind(meta, cbind(d, q, l, filesizeA, filesizeB, positives, meanentropy, meanengrams, sdgrams, skew, ginicoefficient,ginicoefficientqgrams, errorestimation , uniqueness,sduniqueness, meanmissing, missingamount,m,u))
    write.table(meta, "tempResult2.csv", sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)
 }
}


# format result
result <- as.data.frame(meta)

# Write final result file
write.table(result, "fullResult_Metadata.csv", sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)

