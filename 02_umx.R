setwd("~/GitHub/twin-data-models")
rm(list = ls())

# Load libraries
library(umx)
library(tidyverse)

# Set up data
famData <- readRDS("QTAB_familywise.RDS")
mzData <- famData %>% filter(zyg == 1 | zyg == 2)
dzData <- famData %>% filter(zyg == 3 | zyg == 4 | zyg == 5)

# Session 1, NIH Toolbox Cognition Battery, Processing Speed (raw)
## No covariates
m1 <- umxACE(selDVs = "ProcSpeed_raw", sep = "_0", dzData = dzData, mzData = mzData)
m1$top$a_std$result^2 # A (square to get % variance explained)
m1$top$c_std$result^2 # C (square to get % variance explained)
m1$top$e_std$result^2 # E (square to get % variance explained)
# Estimates also available from umxSummaryACE
umxSummaryACE(m1, std = T)

## Covariates (age, sex) regressed out prior to twin modelling
famData.resid = umx_residualize(c("ProcSpeed_raw"), cov = c("ses01_age_months", "sex_female"), data = famData, suffixes = c("_01", "_02")) # regress out age & sex
mzData.resid <- famData.resid %>% filter(zyg == 1 | zyg == 2)
dzData.resid <- famData.resid %>% filter(zyg == 3 | zyg == 4 | zyg == 5)
m2 <- umxACE(selDVs = "ProcSpeed_raw", sep = "_0", dzData = dzData.resid, mzData = mzData.resid)
m2$top$a_std$result^2 # A (square to get % variance explained)
m2$top$c_std$result^2 # C (square to get % variance explained)
m2$top$e_std$result^2 # E (square to get % variance explained)
# Estimates also available from umxSummaryACE
umxSummaryACE(m2, std = T)
