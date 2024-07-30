########################################################
## Script Figure_Suppl_1_FPCA_metasPeritumor_3FPC.R  ###
########################################################
# Updated 07/28/2024 

rm(list = ls())

require(devtools)  # to install PartiallyFD 
require(mclust)    # for PartiallyFD
require("fda")     # for PartiallyFD
library("devtools") # for PartiallyFD
library("mclust") # for PartiallyFD
library("fda")    # for PartiallyFD
install_github("stefanrameseder/PartiallyFD") # just for a function (PartiallyFD:::addAlpha) 
# for nice semi-transparent curves in the row curve plot of part 3 (optional if you accept non-semi transp. curves look)
library("PartiallyFD") 
library(fdadensity) # main package for FPCA compatible with CT density distributions

##############
# read the files with all smooth CT density curves extracted from the segmentation of the peri-metastases regions.
metastases_io_files <- "C:\\Users\\psalazar\\Desktop\\metastases_io_files"
originfile <- paste(metastases_io_files,"\\All_metasPeritumor_fda_out_v3.csv",sep="")

metascurves0 <- read.csv(file=originfile,sep=",",header = TRUE,fileEncoding ="UTF-8-BOM")
head(metascurves0)

metascurves <- (t(metascurves0[,-1]))## transpose matrix to get matrix like example in fdadensity package example. We remove the ctHU values ranging from -1000HU to 499HU to have a matrix with curves only.
## head(metascurves)
## tail(metascurves)
dim(metascurves)  ## rows columns

metasrow <- nrow(metascurves)## OK  currently 111 (TGH metastase cases)
metascol <- ncol(metascurves) ## OK 1500 since we have [-1000HU to 499HU] range

## rownames(metascurves)## case names(subject number)- diagnostic
## colnames(metascurves) ## diagnostic - NULL is ok

################################################################################
## Part 1 - Normalize/ check the constraint of integral to 1 for each row - case
# Each CT density curve area should be strictly equal to 1 for fdadensity FPCA algorithms.
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
# Option 1: Dense data (not sparse). option 2: We assume no errors. Option 3: methodSelectK=4 (we need 4 functional principal components=4 (we don't even report the 4th component in the metastase study) based on pilot study and previous papers))
# FVE: Fraction-of-Variance-Explained
# str(X) ## Diagnostic result matrix

## Tittle 4 Plots FPCA #########################################################################
textplot1 <- paste("cumul.FVE:",round(100*X$cumFVE[1],1),"%"," Peritumor -1st Mode")
textplot2 <- paste("cumul.FVE:",round(100*X$cumFVE[2],1),"%"," Peritumor -2nd Mode")
textplot3 <- paste("cumul.FVE:",round(100*X$cumFVE[3],1),"%"," Peritumor -3rd Mode")
textplot4 <- paste("cumul.FVE:",round(100*X$cumFVE[4],1),"%"," Peritumor -4th Mode")
## End tittle 4 Plots FPCA #####################################################################

# Plot the  4 first functional principal components (FPCs) curves: for 10 percentile, 25th percentile, mean curve, 75th percentile and 90th percentile.
#####################
# Step 1: Call the pdf command to start the plot


# metastases_io_files 
sup_Fig1filename <- paste(metastases_io_files,"\\Suppl_Figure1_PerimetastasesFPCA.pdf",sep="")

pdf(file = sup_Fig1filename,   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

# plot the 4 curves

par(mfrow=c(2,2))
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

## Fourth principal component (not used in the publication - low specific FEV)
Qvec4 = quantile(X$xiEst[,4], probs = c(0.1, 0.25, 0.75, 0.9))/sqrt(X$lambda[4])
CreateModeOfVarPlotLQ2D(X, k = 4, dSup = pointsNumber, Qvec = Qvec4, main = textplot4,xlab="CT density (HU)",ylab="Frequency")

# Run dev.off() to create the file!
dev.off()


##################################################################################
### Part 3: Visualize all smooth row CT density curves for peri-metastase regions 
# together with the classic Euclidean L2 cross sectional mean curve
# and the Wasserstein-Frechet mean. See fdadensity reference for more details.

## No use of alpha for L2 based FEV - suspicious result but this is not the object of our study - only unpublished sup data here
fve.L2 <- GetFVE(fpcaObj=X,dmatrix = metascurves, dSup = pointsNumber, useAlpha = FALSE, metric="L2")
fve.L2 <- round(100*fve.L2,2)
cat(paste("Fract. Variance Explained (L2):",fve.L2,"\n"))

## wasserstein-Frechet based FEV -  reported results
fve.Wass <- GetFVE(fpcaObj=X,dmatrix = metascurves, dSup = pointsNumber, useAlpha = FALSE, metric="W")## we Cannot apply regularization ("does not result in valid density")
fve.Wass <- round(100*fve.Wass,2)
cat(paste("Fract. Variance Explained (Wasserstein):",fve.Wass,"\n"))

## Plot W Frechet mean - Good row plot for curves - alpha = 0.00003 
# Regularize CT densities before computing the Wasserstein mean CT density curve.
metascurvesReg <- t(apply(metascurves, 1, function(u) RegulariseByAlpha(u, x = pointsNumber, alpha = alpha))) 

# alpha = 0.00003  ##
pointsNumber <- c(1:metascol)
pointsNumber <- c(-1000:499)

min(pointsNumber)
wfmean <- getWFmean(dmatrix = metascurves, dSup = pointsNumber,useAlpha = TRUE,alpha = 0.0001)## get Wasserstein Frechet mean through quantile density averaging
# GetWFmean() computes the mean in the quantile space (using dens2qd() and colMeans()) then backtransform to density space using qd2dens() function


# peri-metastases CT density smooth row curves
sup_Fig1B_filename <- paste(metastases_io_files,"\\Suppl_Figure1B_Perimetastases_RowCurvesandMean.pdf",sep="")

pdf(file = sup_Fig1B_filename,   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches


# Display plot with all the peri-metastase CT density curves (smooth) and the L2 mean curve and Wasserstein-Frechet mean curve.
par(mfrow=c(1,1)) ## Single plot
matplot(pointsNumber,t(metascurvesReg),xlab="Peri-Metastases CT Density (HU)", ylab="Voxel Frequency",col = PartiallyFD:::addAlpha("black", 0.25), type = "l", lwd = 1, xaxt = "n", yaxt = "n", lty = "solid")
axis(1, xaxp=c(min(pointsNumber), max(pointsNumber)+1, 10), las=2)
lines(pointsNumber,wfmean,lwd=3,col="red")
lines(pointsNumber, colMeans(metascurves), lwd=3,col="blue")
# Add a legend
legend('topright',col=c('black','red','blue'),lty=1,lwd=c(1,2,2),legend=c('Densities (prob.)','Wasserstein-Frechet Mean','Euclidean L2 mean'))
## str(metascurves) ## Diagnostic

# Run dev.off() to create the file!
dev.off()

#############################################################################################

