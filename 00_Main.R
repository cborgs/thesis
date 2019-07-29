####################
#
# Master file to reproduce results from the thesis
# All required packages were stored by packrat
#
####################

# Set working dir ----------------
#setwd("E:/Priv/Dropbox/Dissertation/Abgabe")


# Load libs ---------------------------

require(packrat)
# Only needed once
#packrat::init()
#devtools::install_github("germanrecordlinkage/multibitTree")
#packrat::snapshot()
#packrat::clean() 
#packrat::bundle(file = "Diss_CB_Source.tar.gz",include.lib = TRUE) 

# Load packrat
#packrat::unbundle(file = "Diss_CB_Source.tar.gz") 
packrat::packify()

# Prepare data ---------------------------

print("Make 100k sample of FEBRL WA for evaluation")
source("06_Make_FEBRLWA_100kFile.R")


# Preliminaries and misc. plots -------

print("Plotting Big data and Record linkage plots (Fig. 1.1 and 1.2.)")
source("01_BigData_Scholar_RecordLinkage_medline_soc.R")

print("Plotting all preliminary plots: Variance of repeated BF generation, k/fill/f Pretests with sim and real data and random forest vs. linear regression plot")
source("02_Preliminaries.R")

print("Plot Bigram frequency plot (Fig 4.3)")
source("03_Plot_BigramFreqs.R")

print("Plot time by tanimoto threshold, l and tree type; Plot 1 mio vs 1 mio linkage time by threshold")
print("Warning: Requires compilation of the mbtSearch_2.0.1 folder: Navigate there and execute 'make'. Make sure the resulting program is executable (chmod -X).")
source("04_Plot_Times.R")


# Run main analysis files ---------------------------

# Meta data evaluation of training data
source("08_Datasets_Metadata.R")

# Training data full evaluation, including joining meta data
source("07_TestSetup.R")

# OPC of the three test data sets
source("09_Optimal_Parameter_Choice_FEBRLWA.R")
source("09b_Optimal_Parameter_Choice_NCVoter.R")
source("09c_Optimal_Parameter_Choice_Mortality.R")

# Evaluation of the parameter choices
source("11_OPTC_vs_TrainingResults.R")


# Plot results ---------------------------

# Plot all result plots
source("12_Plots_Chap5.R")


