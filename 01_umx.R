setwd("/qtab/twin-data-models")
rm(list = ls())

# Load libraries
library(umx)
library(tidyverse)

# Setup datafile
qtab.participants <- read.table("/qtab/participants.tsv", sep = "\t", header = T, na.strings = "n/a")
ses01.cog <- read.table("/qtab/phenotype/2_cognition_ses-01.tsv", sep = "\t", header = T, na.strings = "n/a")
qtab.data <- left_join(qtab.participants, ses01.cog, "participant_id")
# Recode sex from M/F to 0/1
qtab.data <- qtab.data %>% mutate(sex_num = if_else(sex=="F", 1, 0))
head(qtab.data %>% select(sex, sex_num))

##### Convert to familywise dataset (i.e. one family per row) ####
# Create an individual ID based on birth order (i.e. Twin 1, Twin 2)
qtab.data$indID <- qtab.data$birth_order

# Recode pairs from triplets so that indID is 1 & 2 (where twin 1 is the first born, and twin 2 is the second born)
qtab.data <- qtab.data %>%                               
  mutate(indID = replace(indID, participant_id == "sub-0081", 2)) # fam-0085
qtab.data <- qtab.data %>%                               
  mutate(indID = replace(indID, participant_id == "sub-0386", 1)) # fam-0210
qtab.data <- qtab.data %>%                               
  mutate(indID = replace(indID, participant_id == "sub-0138", 2)) # fam-0210
qtab.data <- qtab.data %>%                               
  mutate(indID = replace(indID, participant_id == "sub-0355", 2)) #fam-0217
# Check that pairs now have indID 1 & 2
qtab.data[which(qtab.data$family_id=="fam-0085"), c("participant_id", "birth_order", "indID")]
qtab.data[which(qtab.data$family_id=="fam-0210"), c("participant_id", "birth_order", "indID")]
qtab.data[which(qtab.data$family_id=="fam-0217"), c("participant_id", "birth_order", "indID")]

# Reshape from long to wide
famData <- reshape(qtab.data, timevar = "indID", idvar = "family_id", direction = "wide", sep = "_0") 
# Create single zygosity variable
famData$zyg <- famData$zyg_01
# Check to make sure there are no families without zygosity
famData %>% filter(is.na(zyg))

#### Twin Modelling - umx ####
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
m2 <- umxACE(selDVs = "ProcSpeed_raw", selCovs = c("ses01_age_months", "sex_num"), sep = "_0", dzData = dzData, mzData = mzData)
umxSummaryACE(m2, std = T)
0.48^2 # A (square to get % variance explained)
0.48^2 # C (square to get % variance explained)
0.73^2 # E (square to get % variance explained)

## Covariates (age, sex) regressed out prior to twin modelling
famData.resid = umx_residualize(c("ProcSpeed_raw"), cov = c("ses01_age_months", "sex_num"), data = famData, suffixes = c("_01", "_02")) # regress out age & sex
mzData.resid <- famData.resid %>% filter(zyg == 1 | zyg == 2)
dzData.resid <- famData.resid %>% filter(zyg == 3 | zyg == 4 | zyg == 5)
m3 <- umxACE(selDVs = "ProcSpeed_raw", sep = "_0", dzData = dzData.resid, mzData = mzData.resid)
umxSummaryACE(m3, std = T)
0.49^2 # A (square to get % variance explained)
0.48^2 # C (square to get % variance explained)
0.73^2 # E (square to get % variance explained)
