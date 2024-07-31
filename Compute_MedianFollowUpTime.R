# Script: Compute_MedianFollowUpTime.R
# Compute median follow-up time in months and days taking into account the survival effect
# we use the 'reverse'Kaplan-Meier method: from D.Moore's book: 'Applied survival analysis using R'. Chap. 3.3.
# DOI:10.1007/978-3-319-31245-3

rm(list = ls())

# Load packages
library(survival) # for survfit()
library(dplyr) # for glimpse()

# Load metastases data
surv_data <- read.csv("C:\\Users\\psalazar\\Desktop\\metastases_io_files\\Metas111_v3_all_variables.csv",sep=",",header=TRUE)

glimpse(surv_data)

#########################################
set.seed(1962)

FPC1 <- surv_data$FPC1
FPC2 <- surv_data$FPC2
peri_FPC1 <- surv_data$peri_FPC1
peri_FPC3 <- surv_data$peri_FPC3
size <- surv_data$size
age <- surv_data$age
ssf0_entropy <- surv_data$ssf0_entropy
ssf4_entropy <- surv_data$ssf4_entropy
ssf0_skewness <- surv_data$ssf0_skewness
ssf4_skewness <- surv_data$ssf4_skewness


rfs_days<- surv_data$rfs_days
rfs_months <- rfs_days/30.417 # use a more precise rfs_month than original
rfs_event<- surv_data$rfs_event
table(rfs_event) # to check that 1=recurrence or death

# we use the 'reverse'Kaplan-meier method: from D.Moore's book: 'Applied survival analysis'. Chap. 3.3. 

event_followup <- 1-rfs_event
result_surv <- survfit(Surv(rfs_months,event_followup)~1)
resFU <- quantile(result_surv)
# Results in months
median_FU_Months <- round(resFU$quantile[2],1)
median_low_FU_Months <- round(resFU$lower[2],1)
median_upp_FU_Months <- round(resFU$upper[2],1)
# Results in months
median_FU_TimeDays <- resFU$quantile[2]*30.417 # To report in Metastases_Radiomics results - Descriptive.
median_low_FU_Days <- resFU$lower[2]*30.417
median_upp_FU_Days <- resFU$upper[2]*30.417

# Print result median follow-up time - days and months: 
cat("Median follow-up time",median_FU_TimeDays,"days [",median_low_FU_Days,"to",median_upp_FU_Days,"]","\n")
cat("Median follow-up Months",median_FU_Months,"days [",median_low_FU_Months,"to",median_upp_FU_Months,"]")

# Compare with plain median for rfs_months - our median follow-up time is much longer than the median RFS
median(rfs_months)
median(rfs_days)



