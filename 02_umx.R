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
umxSummaryACE(m1, std = T)
0.44^2 # A (square to get % variance explained)
0.59^2 # C (square to get % variance explained)
0.68^2 # E (square to get % variance explained)

## Covariates (age, sex) included in model
m2 <- umxACE(selDVs = "ProcSpeed_raw", selCovs = c("ses01_age_months", "sex_female"), sep = "_0", dzData = dzData, mzData = mzData)
umxSummaryACE(m2, std = T)
0.48^2 # A (square to get % variance explained)
0.48^2 # C (square to get % variance explained)
0.73^2 # E (square to get % variance explained)

## Covariates (age, sex) regressed out prior to twin modelling
famData.resid = umx_residualize(c("ProcSpeed_raw"), cov = c("ses01_age_months", "sex_female"), data = famData, suffixes = c("_01", "_02")) # regress out age & sex
mzData.resid <- famData.resid %>% filter(zyg == 1 | zyg == 2)
dzData.resid <- famData.resid %>% filter(zyg == 3 | zyg == 4 | zyg == 5)
m3 <- umxACE(selDVs = "ProcSpeed_raw", sep = "_0", dzData = dzData.resid, mzData = mzData.resid)
umxSummaryACE(m3, std = T)
0.49^2 # A (square to get % variance explained)
0.48^2 # C (square to get % variance explained)
0.73^2 # E (square to get % variance explained)
