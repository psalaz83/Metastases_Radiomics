#################################################
# Script - SmoothDensityCurve_Metas_fda.R
# This R-script provides a curve smoothing adapted to density distribution - forced to 1 for sigma
# TGH metastases dataset with lambda=2000 and 101 spline bases.

# Generate new smooth histogram using as input the CT histograms called xxx_output2  (single row with CT HU observations: example: 1024,1024,...,45,45,45,46,48,48,50,..499) 
# Output file: _output_Smooth.csv to be used in FPCA or other stat. analysis
# Accept single number as filename key (_output2 files). Generate a density smooth curve file _output_smooth.
# smoothing: define curve from -1000HU to 499HU
##################


rm(list = ls())
require(fda) # fda R package for the smoothing curve methods with constraint (CT density curves like any frequency distributions sum up to 1 (area under the curve)
# this requires a special method to smooth the curves - see reference fda in metastase paper and supplementary materials.

library(fda)

# read exported histogram from Vitrea - Small region (example: metas segmented in lung nodule analysis).

# A dictionary matches the cases number withthe case name from the oncologists
path_diction <- "C:\\Users\\psalazar\\Desktop\\metastases_io_files\\ThreeSmoothCurves\\Alldata_111tumor_filenames_v3.csv"
#######################
# Enter here the vector of all the case numbers. 
dict_casenumber_name <- read.csv(file= path_diction,sep=",",header=TRUE,fileEncoding ="UTF-8-BOM" )
dict_casenumber_name$case_num # list of cases number

 metas.cases <- dict_casenumber_name$case_num # we use a small list of three cases
 metas.cases <- c(109:111) # To demonstrate the algorithm, we use a small list of only three cases: #109, #110 and #111.

val <- 0
for (val in metas.cases){ 
  cat("Case:",val,"\n")
  inputFilenumber <- as.character(val)
  inputFileName <- dict_casenumber_name$casename_tumor[val]
  
  #######################
  filename1 <- paste("C:\\Users\\psalazar\\Desktop\\metastases_io_files\\ThreeSmoothCurves\\",inputFileName, sep="")
  filename.source <- paste(filename1,"_output2.csv",sep="")

##############################################################
proc_values <- read.csv(file=filename.source,sep=",",header=TRUE,fileEncoding ="UTF-8-BOM")
##############################################################
# the algorithm is presented in fda R package reference and supplementary material of the metastase paper.
  
x <- proc_values$x
minx <- min(x)
maxx <- max(x)

par(mfrow = c(3,1))  # plots with 3 rows 2 columns

#  set up range for density
rangeval <- c(-1000,499) # Typical useful range for the metastases CT HU densities

# set up basis for W(x)
Wbasis <- create.bspline.basis(rangeval, 155) # 305 (1000+500+5)/10 is the number of bases. 4 is implicit order of the splines
plot(Wbasis, main=paste("Case:",inputFilenumber," B-Spline bases: 155"))

#  set up initial value for Wfdobj
Wfd0 <- fd(matrix(0,155,1), Wbasis)
WfdParobj <- fdPar(Wfd0,lambda= 5*10^4)

#  estimate density - key function
denslist <- density.fd(x, WfdParobj,dbglev=2,iterlim=9)
# Gives the optimal function to be minimized ($f) and the norm of the gradient vector at the optimal solution ($norm)

###### Plot the  convergence criterium - diagnostic plot ##
# Check values 
denslist$Flist$f
str(denslist)

iternum <- denslist$iterhist[,1]
maxx <- max(iternum)

itercriter <- denslist$iterhist[,2]
miny <- floor(min(itercriter))
maxy <- ceiling(max(itercriter))+ 100

plot(iternum,itercriter,xlim=c(0,maxx),ylim =c(miny,maxy), col="red",main=paste("Criterium:",miny) )
lines(iternum,itercriter,col="blue")


### plot ######
# We backtransform the density function to the usual space with exp()
# plot density
xval <- seq(-1000,499,1)#from -1000HU to 499HU
wval <- eval.fd(xval, denslist$Wfdobj)
pval <- exp(wval)/denslist$C

# length(x)
# length(pval)
# length(wval)

plot(xval, pval, type="l", ylim=c(0,0.005),main=paste("Case:",inputFilenumber)) # we will choose the optimum density max later
points(x,rep(0,length(x)))

# write a new file _output_Smooth.csv with the smoothed values 2 columns and length(xval) rows with headers X Y#

mat <- matrix(c(xval,pval),nrow=length(xval)) 

filename3 <- paste(filename1,"_output_smoothC.csv",sep="")
case.numb <- paste("case_",inputFilenumber,sep="")


#######################################################################
write.table(mat, file=filename3,row.names=FALSE, na="",col.names=c("ctHU",case.numb), sep = ",",append = FALSE)
#######################################################################

}


