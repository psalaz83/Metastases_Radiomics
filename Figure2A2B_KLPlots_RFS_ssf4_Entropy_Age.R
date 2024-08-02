
## R-Script for computing Figure2_KLPlots_RFS_ssf4_Entropy_Age.R

## works with file: Metas111_v3_all_variables.csv
## The analysis performs Kaplan-Meier plots, computes median survival time with CI95%, log rank tests for difference in survival curves.
# ======================================
rm(list = ls())

## Load required packages
library(survival)
library(ggplot2)
library(survminer) ## to plot the KM plots

#######################################################################################################################
## Function for the p-values of the log-rank test
# Log-rank test:to better isolate the p-value
pval_km_test  <- function(fit_variable_km,iter,datasurv=surv1){ 
pval_km_var<- survminer::surv_pvalue(fit_variable_km, datasurv)
vec_res <-  c(iter,pval_km_var$pval)## isolate the p-value
return(vec_res)   
}
####end function ############################################################################################################
## Load metastases data
surv0 <- read.csv("C:\\Users\\psalazar\\Desktop\\metastases_io_files\\Metas111_v3_all_variables.csv",sep=",",header=TRUE)
head(surv0)
surv1 <- surv0[,-c(24,31,56:60)] ## we exclude the 3 useless variables "_total"
head(surv1)

listpredictorskm <- colnames(surv1[,3:54]) ## All predictors testing for logrank test and Kaplan-Meier plots
length(listpredictorskm)
# Use months as time scale
surv1$time=surv1$rfs_days/30.417
surv1$status=surv1$rfs_event
head(surv1)
## Add survival object. status == 1 is not rfs
surv1$SurvObj <- with(surv1, Surv(time, status == 1))

############################
## Collect table with results KM analysis - P-values

maxiter=length(listpredictorskm)  ## currently 52 - the number of variables tested for KM analysis
# Since we are showing the KM plots and corresponding result variables only for Figures 2A 2B, that is Age and ssf_4 - RFS, we set maxiter to 2.
maxiter <- 2
res.names <- c("Val","P-value")
mat.results <- matrix(data=0.0, nrow = maxiter, ncol=2, dimnames= list(c(1:maxiter),res.names))

str(mat.results)
mat.results[,2]
# 
#######################################################################################################
# age - get quantiles
##########
iter=1

quant.age <- quantile(surv1$age)
quant0 <- round(quant.age[[1]],2)
quant50 <- round(quant.age[[3]],2)
quant100 <- round(quant.age[[5]]+0.01,2)

# First we create a categorical variable for age (approximately according to median)
surv1$ageGroup=cut(surv1$age,breaks=c(quant0,quant50,quant100), labels=1:2)

############
## Km plot
fit.age.km=survfit(Surv(time,status)~ageGroup,data=surv1, conf.type="plain",conf.int=T)

ggsurvplot(fit.age.km,conf.int=T,pval=T,pval.size=4,xlab="Months",ylab="RFS",legend.title="age",legend.labs=c("Below median","median and above"),main="age")
ggsave("KMPlot_RFS_age.pdf",width = 8, height = 6, path="C:\\Users\\psalazar\\Desktop\\metastases_io_files") # save KM plot

# Estimates of two years survival probabilities with "plain" confidence intervals 
summary(fit.age.km,time=24) 
# Estimates of median survival time 
print(fit.age.km)

######################################################
## call function for p-value of fit km for each variable
respval <- pval_km_test(fit_variable_km=fit.age.km ,iter=iter,datasurv=surv1)
# # Log-rank test:to better isolate the p-value
mat.results[iter,1]<- respval[1] ## we fill the result matrix with iter and p-value (log-rank test)
mat.results[iter,2]<- respval[2]
#######################################################################################################
# ssf4_entropy
##########
iter=iter+1

quant.ssf4_entropy <- quantile(surv1$ssf4_entropy)
quant0 <- round(quant.ssf4_entropy[[1]],2)
quant50 <- round(quant.ssf4_entropy[[3]],2)
quant100 <- round(quant.ssf4_entropy[[5]]+0.01,2)

# First we create a categorical variable for ssf4_entropy (approximately according to median)
surv1$ssf4_entropyGroup=cut(surv1$ssf4_entropy,breaks=c(quant0,quant50,quant100), labels=1:2)

############
## Km plot
fit.ssf4_entropy.km=survfit(Surv(time,status)~ssf4_entropyGroup,data=surv1, conf.type="plain",conf.int=T)

ggsurvplot(fit.ssf4_entropy.km,conf.int=T,pval=T,pval.size=4,xlab="Months",ylab="RFS",legend.title="ssf4_entropy",legend.labs=c("Below median","median and above"),main="ssf4_entropy")

ggsave("KMPlot_RFS_ssf4_entropy.pdf",width = 8, height = 6, path="C:\\Users\\psalazar\\Desktop\\metastases_io_files")

# Estimates of two years survival probabilities with "plain" confidence intervals 
summary(fit.ssf4_entropy.km,time=24) 

# Estimates of median survival time 
print(fit.ssf4_entropy.km)

######################################################
## call function for p-value of fit km for each variable
respval <- pval_km_test(fit_variable_km=fit.ssf4_entropy.km ,iter=iter,datasurv=surv1)
# # Log-rank test:to better isolate the p-value
mat.results[iter,1]<- respval[1] ## we fill the result matrix with iter and p-value (log-rank test)
mat.results[iter,2]<- respval[2]
#######################################################################################################


