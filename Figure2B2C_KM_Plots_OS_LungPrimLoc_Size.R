
# R-Script: Figure2B-C_KM_Plots_OS_LungPrimLoc_Size.R
# Km plots for 'lung - primary lesion origin' and 'Size'. 

rm(list = ls())

## Load required packages
library(survival)
library(ggplot2)
library(survminer) ## to plot the KM plots

## Load metastases data - version with LungPrimLoc binary predictor (..v4..)
surv1 <- read.csv("C:\\Users\\psalazar\\Desktop\\metastases_io_files\\Metas_v4_all_variables.csv",sep=",",header=TRUE)
head(surv1)

##################################################################################
set.seed(825)
size <- surv1$size
clinic_indic_bin <- as.factor(surv1$clinic_indic_bin)
lung_primorigin <- as.factor(surv1$lung_primorigin)
class(lung_primorigin)

os_days<- surv1$os_days
os_months <- surv1$os_days/30.417 # use a more precise os_month than original
os_event <- surv1$os_event
table(os_event) # check the correctness of the survival event values
S <- Surv(os_months,os_event)

#############################################################################################################
## Function for the p-values of the log-rank test
# Log-rank test:to better isolate the p-value
pval_km_test  <- function(fit_variable_km,iter,datasurv=surv1){   ##pval_km_ssf6_entropy is the test variable 
  pval_km_var<- survminer::surv_pvalue(fit_variable_km, datasurv)
  vec_res <-  c(iter,pval_km_var$pval)## isolate the p-value
  return(vec_res)   
}
####end function ############################################################################################################
surv1$time=surv1$os_days/30.417
surv1$status=surv1$os_event

table(surv1$os_event)
summary(surv1$os_event)
## Add survival object. status == 1 is not os
surv1$SurvObj <- with(surv1, Surv(time, status == 1))

#KM survival Plot
fit.lungprim.km= survival::survfit(Surv(time,status)~lung_primorigin,data=surv1, conf.type = "log-log")
ggsurvplot(fit.lungprim.km,conf.int=T,pval=T,pval.size=4,xlab="Months", ylab="OS",legend.title="Primary Lesion Origin.",legend.labs=c("Other","Lung"),main="Prim. lesion site")
ggsave("KM_Plot_OS_lungOrigin.pdf",width = 8, height = 6, path="C:\\Users\\psalazar\\Desktop\\metastases_io_files")

# Estimates of one year survival probabilities with "plain" confidence intervals 
sumtest <- summary(fit.lungprim.km,time=1) 
str(sumtest)
# Estimates of median survival time - If one class does not reach 50% OS, median = NA.
print(fit.lungprim.km)


#######################################################################################################
# KM analysis for predictor: size 
##########
quant.size <- quantile(surv1$size)
quant0 <- round(quant.size[[1]],2)
quant50 <- round(quant.size[[3]],2)
quant100 <- round(quant.size[[5]]+0.01,2)

# First we create a categorical variable for size (approximately according to median)
surv1$sizeGroup=cut(surv1$size,breaks=c(quant0,quant50,quant100), labels=1:2)

############
## Km plot
fit.size.km=survival::survfit(Surv(time,status)~sizeGroup,data=surv1, conf.type="plain",conf.int=T)


ggsurvplot(fit.size.km,conf.int=T,pval=T,pval.size=4,xlab="Follow-up time (Months)", ylab="Survival - OS",legend.title="size",legend.labs=c("Below median","median and above"),main="size")
ggsave("KM_Plot_OS_size.pdf", width = 8, height = 6, path="C:\\Users\\psalazar\\Desktop\\metastases_io_files")


# Estimates of two years survival probabilities with "plain" confidence intervals 
summary(fit.size.km,time=2) 

# Estimates of median survival time 
print(fit.size.km)

# Log-rank test:
survdiff(Surv(time,status)~surv1$sizeGroup,data=surv1)

#######################################################################################################


