#########################################################################
#
# Estimating optimal identifier choices for input files specified
# 
# Input: At least one readily preprocessed input file. If two files, both must have the same column names. 
# Recommended: One or two files with IDs.
# Packages required: PPRL, fastdigest, stringdist, entropy, data.table, fastLink, moments, ineq
# User-defined settings in lines 21 to 45.
#
# This script reads the file(s), generates all measures required for the estimation, 
# and outputs the optimal identifier choice.
#
# Output: 
#    IdentifierChoice.csv 
#
# Christian Borgs, 2019
#
#########################################################################


#### User settings ####

# Set working dir for script
setwd("E:/Priv/Dropbox/Dissertation/Programme")
#setwd("~/Dropbox/Dissertation/Programme")

# File to run the optimization over
inputFile <- "Daten/A_10000_python_R0_C0_O100_percent.csv"
# file name of file to link to
inputFileB <- "Daten/B_10000_python_R20_C50_O100_percent.csv"

# Identifiers: Smalles set to test
identsStandard <- c("V7","V8")
# Identifiers: All possible
identsMax <- c("V7","V8","V3","V4","V5","V6","V2")
#identsMax <- c("V7","V8","V5","V4","V3","V2","V6","V9")
# ID-column. If none, set to FALSE (WARNING: This will be detrimental to the model quality!)
IDcol <- "V1"
# String-Dinstance columns. For EM weighting
Stringidents <- c("V7","V8")

# Percent of file to sample fromm if file large (1 = no sampling)
trainfraction <- 1
# Cores of the machine for multithreading
cores <- 6



#######################

# Libs used
library(PPRL)
library(data.table)
library(fastLink)
library(rlist)

#####################

## Automate list of identifiers
xyc <- list(identsStandard)
lenMin <- length(identsStandard)
for (j in lenMin:length(identsMax)) {
  id <- c(identsStandard, identsMax[(lenMin + 1): (j + 1)])
  id <- c(na.omit(id))
  xyc <- list.append(xyc, id)
}
identsets <- unique(xyc)
# 
# # Temporary
# identsets <- list(
#   c("Vorname","Nachname"),
#   c("Vorname","Nachname", "Day"),
#   c("Vorname","Nachname", "Day", "Month"),
#   c("Vorname","Nachname", "Day", "Month", "Year"),
#   c("Vorname","Nachname", "Day", "Month", "Year", "DOB"),
#   c("Vorname","Nachname", "Day", "Month", "Year", "RAND"),
#   c("Vorname","Nachname", "Day", "Month", "Year", "RAND2")
# )

#### Work with data

# Read input
clearTextA <- fread(inputFile, colClasses = "character", stringsAsFactors = FALSE)
clearTextB <- fread(inputFileB, colClasses = "character", stringsAsFactors = FALSE)


clearTextA <- clearTextA[sample(nrow(clearTextA),round(trainfraction*nrow(clearTextA)), replace = FALSE),]
clearTextB <- clearTextB[sample(nrow(clearTextB),round(trainfraction*nrow(clearTextB)), replace = FALSE),]


# Check file sizes
filesizeA <- as.numeric(nrow(clearTextA))
filesizeB <- as.numeric(nrow(clearTextB))


  # Generate ID Col
  IDA <- as.character(unlist(clearTextA[,..IDcol]))
  IDB <- as.character(unlist(clearTextB[,..IDcol]))

  
  clearTextAbak <- clearTextA
  clearTextBbak <- clearTextB
  
  bestIdentifier <- 0
  #bestIdentifierF <- 0
  #bestIdentifierFName <- ""
  bestIdentifierName <- ""
  
  tests <- NULL
  
for (idents in identsets){
  
  cat("\n\n")
  print(idents)

  # subset the data
  clearTextA <- clearTextAbak[, ..idents]
  clearTextB <- clearTextBbak[, ..idents]

  cat("\nActual overlap:", sum(IDA %in% IDB),"\n")
  
  ####### Start working ########

  # Get EM weights (silently)
  capture.output(ests <- fastLink(clearTextA, clearTextB, varnames =  idents, 
                 stringdist.match=Stringidents[Stringidents %in% idents],  stringdist.method="jw",  n.cores = cores, estimate.only = FALSE), file='NULL')
  # Weights
  #m <- ests$p.m
  #u <- ests$p.u
  
  # Weights
  m <- ests$EM$p.m
  u <- ests$EM$p.u
  
  meanAW <- log2(m/u)
  meanDAW <- log2((1-m)/(1-u))
  
  print("Agree/Disagree-Weights")
  print(cbind(meanAW,meanDAW))

  
  
  print("Global m and u")
  print(cbind(m,u))
  
  print("Matrix for m")
  #print(matrix(unlist(ests$p.gamma.k.m), ncol=2, nrow=length(idents), byrow=TRUE, dimnames=list(idents)))
  print(matrix(unlist(ests$EM$p.gamma.k.m), ncol=2, nrow=length(idents), byrow=TRUE, dimnames=list(idents)))
  
  print("Matrix for u")
  #print(matrix(unlist(ests$p.gamma.k.u), ncol=2, nrow=length(idents), byrow=TRUE, dimnames=list(idents)))
  print(matrix(unlist(ests$EM$p.gamma.k.u), ncol=2, nrow=length(idents), byrow=TRUE, dimnames=list(idents)))
  
  #
   positives <- sum(IDA %in% IDB)
   pairs <- filesizeA * filesizeB
  #
  # # Calc measures
   tp <- sum(IDA[ests$matches[,1]] == IDB[ests$matches[,2]])
   fp <- sum(IDA[ests$matches[,1]] != IDB[ests$matches[,2]])
   fn <- positives - tp
   tn <- pairs - positives + fp
   precision <- tp / (tp + fp)
   sens <- recall <- tp /(tp + fn)
   spec <- tn / (tn + fp)
   Fmeasure <- mean(precision,recall)
  #
   print("F-Measure, TP, FP")
   print(cbind(Fmeasure, tp, fp))

  #if(Fmeasure > bestIdentifierF){
  #  bestIdentifierF <- Fmeasure
  #  bestIdentifierNameF <- idents
  #}
  if(abs(meanAW) > bestIdentifier){
    bestIdentifier <- abs(meanAW)
    bestIdentifierName <- idents
  }
  tests <- rbind(tests, data.frame(idents = paste0(idents, collapse = ","),Fmeasure, recall, precision, m,u, meanAW, meanDAW))
  #tests <- rbind(tests, data.frame(idents = paste0(idents, collapse = ","), m,u, meanAW, meanDAW))
  write.table(tests,"IdentifierChoice.csv",sep="\t",row.names=FALSE)
}


#print("Final Matrix for m")
#print(matrix(unlist(ests$p.gamma.k.m), ncol=2, nrow=length(idents), byrow=TRUE, dimnames=list(idents)))

#print("Final Matrix for u")
#print(matrix(unlist(ests$p.gamma.k.u), ncol=2, nrow=length(idents), byrow=TRUE, dimnames=list(idents)))

cbind(aggregate(abs(tests$meanAW), by=list(tests$idents), mean),aggregate(abs(tests$Fmeasure), by=list(tests$idents), mean)$x)



cat("\n Best identifier combination:", bestIdentifierName, "\tAgreement weight: ",bestIdentifier)  

