########################################################
## Script Figure1_FPCA_Metastases_3FPC.R  ###
########################################################
# Updated 07/28/2024 

rm(list = ls())

# require(devtools)  # to install PartiallyFD 
# require(mclust)    # for PartiallyFD
# require("fda")     # for PartiallyFD
# library("devtools") # for PartiallyFD
# library("mclust") # for PartiallyFD
# library("fda")    # for PartiallyFD
# install_github("stefanrameseder/PartiallyFD") # just for a function (PartiallyFD:::addAlpha) 
# for nice semi-transparent curves in the row curve plot of part 3 (optional if you accept non-semi transp. curves look)
# library("PartiallyFD") 
library(fdadensity) # main package for FPCA compatible with CT density distributions

########################
##############
# read the files with all smooth CT density curves extracted from the segmentation of the metastases regions.
metastases_io_files <- "C:\\Users\\psalazar\\Desktop\\metastases_io_files"
originfile <- paste(metastases_io_files,"\\All_metas_fda_out_v3.csv",sep="")

metascurves0 <- read.csv(file=originfile,sep=",",header = TRUE,fileEncoding ="UTF-8-BOM")
head(metascurves0)

metascurves <- (t(metascurves0[,-1]))## transpose matrix to get matrix like Bacteria example. We remove the ctHU values ranging from -1000HU to 499HU to have a matrix with curves only.
## head(metascurves)
## tail(metascurves)
dim(metascurves)  ## rows columns

metasrow <- nrow(metascurves)## 111 cases
metascol <- ncol(metascurves) ## 1500 since we have [-1000HU to 499HU] range

## rownames(metascurves)## case names(subject number)- diagnostic
## colnames(metascurves) ## diagnostic - NULL is ok

################################################################################
## Part 1 - Normalize/Check the constraint of integral to 1 for each row (case).
## This is very important that the density curves after smoothing have an area strictly = 0.
for (i in 1:metasrow) { 
  metascurves[i,] = metascurves[i,]/fdadensity:::trapzRcpp(1:ncol(metascurves), metascurves[i,])
}
###############################################################################################
# Part 2: Perform density FPCA for metas CT densities
alpha = 0.00003  ##
pointsNumber <- c(1:metascol)
pointsNumber <- c(-1000:499)

X = FPCAdens(dmatrix = metascurves, dSup = pointsNumber, useAlpha = TRUE,alpha=alpha,
             optns = list(dataType = 'Dense', error = FALSE, methodSelectK =4,verbose=TRUE))

## Comments options FPCAdens function ##
# added a very small alpha (0.00003) to replace zero density values with very small ones for regularization.
# dsup: support grid for Density domain. Here pointsNumber: point for measures (768) 
# useAlpha: regularization should be performed.
# optns: list of options for FPCA. See FPCA of the fdapace package
# Option 1: Dense data (not sparse). 
# Option 2: We assume no errors. 
# Option 3: methodSelectK=4 (we need 4 functional principal components=4 (we don't even report the 4th component in the metastase study) based on pilot study and previous papers))
# FVE: Fraction-of-Variance-Explained - below we compute the cumulative FVE and diaply it for each FPC plot
##################################################
## Tittle 4 Plots FPCA ##############################################################
textplot1 <- paste("cumul.FVE:",round(100*X$cumFVE[1],1),"%"," Metastases -1st Mode")
textplot2 <- paste("cumul.FVE:",round(100*X$cumFVE[2],1),"%"," Metastases -2nd Mode")
textplot3 <- paste("cumul.FVE:",round(100*X$cumFVE[3],1),"%"," Metastases -3rd Mode")
## End tittle 4 Plots FPCA ###########################################################

# Plot the  4 first functional principal components (FPCs) curves: for 10 percentile, 25th percentile, mean curve, 75th percentile and 90th percentile.
#####################
# Call the pdf command to start the plot
Fig1filename <- paste(metastases_io_files,"\\Figure1_MetastasesFPCA.pdf",sep="")

pdf(file = Fig1filename,   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 14) # The height of the plot in inches

# plot the 3 sets of curves 

par(mfrow=c(3,1))
# Plot Modes
Qvec1 = quantile(X$xiEst[,1], probs = c(0.1, 0.25, 0.75, 0.9))/sqrt(X$lambda[1])
CreateModeOfVarPlotLQ2D(X, k = 1, dSup = pointsNumber, Qvec = Qvec1, main = textplot1,xlab="CT density (HU)",ylab="Frequency")
## take X returned by FPCA on the log quantile function and plot in the density space D (default)
## Nu is the mean of LQD space and Psi is the LQD transform. Psi-1 is the back transform from log quantile back to densities
## Same plot before backtransform - in the LQD space

## Second principal component
Qvec2 = quantile(X$xiEst[,2], probs = c(0.1, 0.25, 0.75, 0.9))/sqrt(X$lambda[2])
CreateModeOfVarPlotLQ2D(X, k = 2, dSup = pointsNumber, Qvec = Qvec2, main = textplot2,xlab="CT density (HU)",ylab="Frequency")

## Third principal component
Qvec3 = quantile(X$xiEst[,3], probs = c(0.1, 0.25, 0.75, 0.9))/sqrt(X$lambda[3])
CreateModeOfVarPlotLQ2D(X, k = 3, dSup = pointsNumber, Qvec = Qvec3, main = textplot3,xlab="CT density (HU)",ylab="Frequency")


# Run dev.off() to create the file!
dev.off()
##################################################################################
