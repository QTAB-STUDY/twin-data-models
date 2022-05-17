setwd("~/GitHub/twin-data-models")
rm(list = ls())

# Load libraries
library(umx)
library(tidyverse)

#### Session 1 ####
# Set up data
famData.ses01 <- readRDS("QTAB_familywise_ses01.RDS")
mzData.ses01 <- famData.ses01 %>% filter(zyg == 1 | zyg == 2)
dzData.ses01 <- famData.ses01 %>% filter(zyg == 3 | zyg == 4 | zyg == 5)

# NIH Toolbox Cognition Battery, Processing Speed (raw)
## No covariates
m1 <- umxACE(selDVs = "ProcSpeed_rawZ", sep = "_0", mzData = mzData.ses01, dzData = dzData.ses01)
m1$top$a_std$result^2 # A (square to get % variance explained)
m1$top$c_std$result^2 # C (square to get % variance explained)
m1$top$e_std$result^2 # E (square to get % variance explained)
# Estimates also available from umxSummaryACE
umxSummaryACE(m1, std = T)

## Covariates (age, sex) regressed out prior to twin modelling
famData.ses01.resid = umx_residualize(c("ProcSpeed_rawZ"), cov = c("ses01_age_months", "sex_female"), data = famData.ses01, suffixes = c("_01", "_02")) # regress out age & sex
mzData.ses01.resid <- famData.ses01.resid %>% filter(zyg == 1 | zyg == 2)
dzData.ses01.resid <- famData.ses01.resid %>% filter(zyg == 3 | zyg == 4 | zyg == 5)
m2 <- umxACE(selDVs = "ProcSpeed_rawZ", sep = "_0", mzData = mzData.ses01.resid, dzData = dzData.ses01.resid)
m2$top$a_std$result^2 # A (square to get % variance explained)
m2$top$c_std$result^2 # C (square to get % variance explained)
m2$top$e_std$result^2 # E (square to get % variance explained)
# Estimates also available from umxSummaryACE
umxSummaryACE(m2, std = T)

#### Session 2 ####
# Set up data
famData.ses02 <- readRDS("QTAB_familywise_ses02.RDS")
mzData.ses02 <- famData.ses02 %>% filter(zyg == 1 | zyg == 2)
dzData.ses02 <- famData.ses02 %>% filter(zyg == 3 | zyg == 4 | zyg == 5)

# NIH Toolbox Cognition Battery, Processing Speed (raw)
## No covariates
m3 <- umxACE(selDVs = "ProcSpeed_rawZ", sep = "_0", mzData = mzData.ses02, dzData = dzData.ses02)
m3$top$a_std$result^2 # A (square to get % variance explained)
m3$top$c_std$result^2 # C (square to get % variance explained)
m3$top$e_std$result^2 # E (square to get % variance explained)
# Estimates also available from umxSummaryACE
umxSummaryACE(m3, std = T)

## Covariates (age, sex) regressed out prior to twin modelling
famData.ses02.resid = umx_residualize(c("ProcSpeed_rawZ"), cov = c("ses02_age_months", "sex_female"), data = famData.ses02, suffixes = c("_01", "_02")) # regress out age & sex
mzData.ses02.resid <- famData.ses02.resid %>% filter(zyg == 1 | zyg == 2)
dzData.ses02.resid <- famData.ses02.resid %>% filter(zyg == 3 | zyg == 4 | zyg == 5)
m4 <- umxACE(selDVs = "ProcSpeed_rawZ", sep = "_0", mzData = mzData.ses02.resid, dzData = dzData.ses02.resid)
m4$top$a_std$result^2 # A (square to get % variance explained)
m4$top$c_std$result^2 # C (square to get % variance explained)
m4$top$e_std$result^2 # E (square to get % variance explained)
# Estimates also available from umxSummaryACE
umxSummaryACE(m4, std = T)
