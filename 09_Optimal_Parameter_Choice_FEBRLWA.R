#########################################################################
#
# Estimating optimal parameter choices for the FEBRL WA data
# 
#########################################################################


#### User settings ####

# Set working dir for script
#setwd("E:/Priv/Dropbox/Dissertation/Programme")
# File to run the optimization over
inputFile <- "Daten/Adrian_Perfect_100k.csv"
# Identifiers
idents <- c("Given Name","Family Name","Date Of Birth")
# ID-column. If none, set to FALSE (WARNING: This will be detrimental to the model quality!)
IDcol <- "Group Id"
# String-Dinstance columns. For EM weighting
Stringidents <- c("Given Name","Family Name")

# Desired BF length
l <- 500
# Type: BF or CLK
type <- "CLK"
# Bootstrap repetitions
bootstraps <- 1 
# Desired Tanimoto threshold
tanimotoGoal <- 0.9
# Cores of the machine for multithreading
cores <- 7
# Mode: "Train" for using a subset of the original data set as training data, "Link" if both files are available
mode <- "Link"
# If mode "Link": relative size of the linking subsample (for large files)
linkSample <- 1
# If mode "Train": relative size of the training subsample
trainSample <- 0.2
# If mode "Link": file name of file to link to
inputFileB <- "Daten/Adrian_10percent_100k.csv"

#######################

# Libs used
library(PPRL)
library(fastdigest)
library(stringdist)
library(entropy)
library(data.table)
library(fastLink)
library(moments)
library(ineq)

#### Static vars, change only if you know what you are doing ####

# Training data BF length
lTraining <- 500
# Training data file
trainFN <- "Complete_Results.csv"
# q-grams to test
Q <- c(1,2,3,4)
# Hash function test space
K <- seq(from = 1, to = 40, by = 1)
# MBT leaf limit
leaflimit <- 3
# Calc number of combinations in total
numCombinations <- length(K) * length(Q) * bootstraps
# Set counter for progress report
counter <- 0
# Result vector initialization
meta <- NULL

#####################


# Read training file for model building
traindata <- read.csv(trainFN, sep="\t")

# Subset to fixed tanimoto threshold
fixed <- subset(traindata, Tani == tanimotoGoal)

# Fixed model calls

BIC <- glm(Fmeasure ~ k + q + errorestimation + skew + meanentropy + 
  log(u) + ginicoefficient + poly(k * q * meanengrams * hammingweight, 
                                  2) + meanmissing + uniquepatternsA,
  data =fixed)
LM <- lm(
  Fmeasure ~  k + q + 
       errorestimation*skew + uniqueness + ginicoefficient + 
    poly(k*q*meanengrams*hammingweight, 2) + 
    uniquepatternsA + m, data=fixed)

#### Work with data

# Read input
clearTextA <- fread(inputFile, colClasses = "character", stringsAsFactors = FALSE)

# If second file, read it. Otherwise, sample
if (mode == "Link"){
  clearTextB <- fread(inputFileB, colClasses = "character", stringsAsFactors = FALSE)
  
  clearTextA <- clearTextA[sample(nrow(clearTextA),round(nrow(clearTextA)*linkSample)),]
  
  clearTextB <- clearTextB[sample(nrow(clearTextB),round(nrow(clearTextB)*linkSample)),]
  
} else {
  clearTextB <- clearTextA[sample(nrow(clearTextA),round(nrow(clearTextA)*trainSample)),]
}

## Deduplicate by ID
clearTextA <- clearTextA[!duplicated(clearTextA$`Group Id`),]
clearTextB <- clearTextB[!duplicated(clearTextB$`Group Id`),]



# Check file sizes
filesizeA <- nrow(clearTextA)
filesizeB <- nrow(clearTextB)
# Number of pairs
pairs <- as.numeric(filesizeA) * as.numeric(filesizeB)

# Pairs and positives
if (IDcol != FALSE){
  # Calc true matches
  truepositives <- sum(unlist(clearTextA[,..IDcol]) %in% unlist(clearTextB[,..IDcol]))

  # Generate ID Col
  IDA <- as.character(unlist(clearTextA[,..IDcol]))
  IDB <- as.character(unlist(clearTextB[,..IDcol]))
  
  # subset the data
  clearTextA <- clearTextA[, ..idents]
  clearTextB <- clearTextB[, ..idents]
} else {
  
  # Generate crude TP estimate by using keys
  clearTextA$key <- apply(clearTextA[,..idents],1,paste0, collapse="")
  clearTextB$key = apply(clearTextB[,..idents],1,paste0, collapse="")

  #Estimate TP
  truepositives <- sum(clearTextA$key %in% clearTextB$key)
  
  # Generate ID Col by generating same IDs for exact matches, and SHA1 random IDs for all other pairs
  clearTextA$IDA <- ""
  set.seed(42)
  clearTextA[which(clearTextA$key %in% clearTextB$key), "IDA"] <- replicate(truepositives,fastdigest(sample(LETTERS, 1000, replace = TRUE)))
  set.seed(42)
  clearTextB$IDB <-  ""
  clearTextB[which(clearTextB$key %in% clearTextA$key), "IDB"] <- replicate(length(clearTextB$key %in% clearTextA$key),fastdigest(sample(LETTERS, 1000, replace = TRUE)))
  clearTextA$IDA[clearTextA$IDA == ""] <-   sha1(as.character(1:length(clearTextA$IDA[clearTextA$IDA == ""] )), key = "A")
  clearTextB$IDB[clearTextB$IDB == ""] <-   sha1(as.character(1:length(clearTextA$IDA[clearTextA$IDA == ""] )), key = "B")

   # Generate ID Col
  IDA <- as.character(clearTextA$IDA)
  IDB <- as.character(clearTextB$IDB)
  
  # subset the data
  clearTextA <- clearTextA[, ..idents]
  clearTextB <- clearTextB[, ..idents]
}


####### Start working ########

# Get EM weights
ests <- fastLink(clearTextA, clearTextB, varnames =  idents, 
                 stringdist.match=Stringidents,  stringdist.method="jw",  n.cores = cores, estimate.only = TRUE)


save(ests, file="EM_ests.Rdata") # In case of crash
load("EM_ests.Rdata")

# Weights
m <- ests$p.m
u <- ests$p.u

# Loop over parameter space
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
    missings <- sum(apply(clearTextA, 2, function(x)sum(is.na(x)))) + sum(apply(clearTextA, 2, function(x)sum(x[!is.na(x)]=="")))
    
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
    
  
    
    for (k in K){
      
      for (i in seq(1:bootstraps)){
      
        cat("\n\n\nEncrypting...", format(Sys.time()," (%a, %d-%m-%Y, %H:%M:%S)"),"\n")
        
        ## Generate password vector for BF encoding, resample for boot-strapping
        PWs <- replicate(length(idents), fastdigest(sample(LETTERS, 10000, replace = TRUE)))
        
  
        encryptedA <- CreateCLK(ID = IDA, clearTextA[, ..idents], k = k, padding = rep(0,length(idents)), q = rep(q,length(idents)), l = l, password = PWs)#,"13124","124124"))
        encryptedB <- CreateCLK(ID = IDB, clearTextB[, ..idents], k = k, padding = rep(0,length(idents)), q = rep(q,length(idents)), l = l, password = PWs)#,"13124","124124"))
  
        
        # Calculate hamming weight and NO of unique patterns
        hammingweight <- mean(c(nchar(gsub("0","", encryptedA$CLKs)),nchar(gsub("0","", encryptedB$CLKs))))
        uniquepatternsA <- length(unique(encryptedA$CLKs))
        uniquepatternsB <- length(unique(encryptedB$CLKs))
        # 
        # # Tree data
        write.table(encryptedA, "TreeA_OPC.csv", sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
        write.table(encryptedB, "TreeB_OPC.csv", sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
        # 
        # Run MBT
        multibitTree.load("TreeB_OPC.csv",threads=cores, leafLimit = leaflimit)
        result <- multibitTree.searchFile("TreeA_OPC.csv", Tani)
        # 
        # 
        # # Classify candidate pairs
        tp <- sum((result$fingerprint) == (result$query))
        fp <- sum((result$fingerprint) != (result$query))
        fn <- truepositives - tp
        tn <- pairs - (fn + tp + fp)
        # 
        # # Calculate central measures
         ReductionRate <-  1 - ((tp+fp)/(tn+fn))
         precision <- tp / (tp+fp)
         recall <- tp / (tp+fn)
         Fmeasure =   (recall + precision) / 2
        
        #tp <- fp <- fn <- tn <- precision <- recall <- Fmeasure <- 0
        
        # Save RAM
        gc()
        
        # Set counter plus one
        counter <- counter + 1
        
        
        cat("\nClassification", counter, "of", numCombinations,"Finished (",format((counter / numCombinations)*100, digits=3),"%)!", format(Sys.time()," (%a, %d-%m-%Y, %H:%M:%S)"), ": \nFilesizes:\t", format(filesizeA,scientific=FALSE),"\nData source and mode:\t" , inputFile, " ", mode,"\nBootstrap:\t\t", i, "/", bootstraps, "\nk/q:\t\t", k, "/",q ,"\nTanimoto:\t", tanimotoGoal, "\n\nPositives:\t", truepositives, "\nTrue Positives:\t", tp, "\nRecall:\t\t", recall, "\nFalse Positives:", fp, "\nPrecision:\t", precision, "\nFalse Negatives:", fn, "\nF-Score:\t", Fmeasure, "\n\n")
        
        # Append to result
        meta  <- rbind(meta, data.frame(data="FEBRL WA", type, q, k, i, l, Tani=tanimotoGoal,Fmeasure, filesizeA, filesizeB, meanentropy, meanengrams, sdgrams, skew, ginicoefficient, ginicoefficientqgrams, errorestimation , hammingweight, uniqueness, sduniqueness, meanmissing, missingamount, uniquepatternsA,uniquepatternsB,m,u))
        write.table(meta, "tempResult.csv", sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)
    }
  }
}


# Prepare resulting data to predict from it
newdata <- data.frame(meta, stringsAsFactors = FALSE, row.names = NULL)
names(newdata)[names(newdata) == "tanimotoGoal"] <- "Tani"
write.table(newdata, "Trainresults.csv", sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)
newdata <- read.csv("Trainresults.csv", sep = "\t", stringsAsFactors = FALSE)

# Predict
newfitted <- data.frame(newdata, predicted = predict(LM, newdata, type = "response"))
newfittedBIC <- data.frame(newdata, predicted = predict(BIC, newdata, type = "response"))

# Get optimal prameter estimate
parameters <- newfitted[which.max((newfitted$predicted)),c("Tani","q","k","l","predicted")]
parameters2 <- newfittedBIC[which.max((newfittedBIC$predicted)),c("Tani","q","k","l","predicted")]

# Write out parameters for use in loop
write.table(parameters, "Optimal_Parameters.csv",sep="\t", row.names = FALSE)
write.table(parameters2, "Optimal_Parameters_BIC.csv",sep="\t", row.names = FALSE)

cat("\n All done!")


