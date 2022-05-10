#### Notes ####
# Author Lachlan T Strike
# Example code for running linear mixed effects models on QTAB data (OpenNeuro DOI)
# Code is provided as a example/starting point. It should not be considered the best practice/approach

#### Setup ####
library(lmerTest)
library(tidyverse)
library(performance)
setwd("~/GitHub/twin-data-models")
rm(list = ls())

#### Create long datafile ####
qtab.data <- read.table("C:/Documents and Settings/uqlstri1/OneDrive - The University of Queensland/3_Papers/03_scientific_data/02_data/01_demographics/participants.tsv", sep = "\t", header = T, na.strings = "n/a")
# Recode sex from M/F to 0/1
qtab.data <- qtab.data %>% mutate(sex_female = if_else(sex=="F", 1, 0))
head(qtab.data %>% select(sex, sex_female))

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

# Create Pair & M variables for linear-mixed models
# See Visscher et al. (2004). The Use of Linear Mixed Models to Estimate Variance Components from Data on Twin Pairs by Maximum Likelihood
qtab.data <- qtab.data %>% arrange(family_id)
qtab.fams <- unique(qtab.data$family_id)
qtab.fams.MZ <- unique(qtab.data[which(qtab.data$zyg<=2), "family_id"])
qtab.fams.DZ <- unique(qtab.data[which(qtab.data$zyg>=3), "family_id"])
qtab.data$Pair <- NA
qtab.data$M <- NA

length(qtab.fams) # 211 families
for(i in 1:212){
  famID <- qtab.fams[i]
  qtab.data[which(qtab.data$family_id==famID & !is.na(qtab.data$zyg)), "Pair"] <- paste0("Pair", i)
}

length(qtab.fams.MZ) #104 MZ families
for(i in 1:104){
  famID <- qtab.fams.MZ[i]
  qtab.data[which(qtab.data$family_id==famID & !is.na(qtab.data$zyg)), "M"] <- paste0("M", i)
}

length(qtab.fams.DZ) # 107 DZ families
M.DZ.start <- 105 # 1 + number of MZ families
M.DZ.end <- 318 # number of MZ families (104) + number of DZ twins (107*2). 104 + 214 = 318
qtab.data[which(qtab.data$family_id %in% qtab.fams.DZ  & !is.na(qtab.data$zyg)), "M"] <- paste0("M", M.DZ.start:M.DZ.end)

head(qtab.data[, c("participant_id", "family_id", "zyg", "indID", "M", "Pair")], n = 30)

ses01.cog <- read.table("C:/Documents and Settings/uqlstri1/OneDrive - The University of Queensland/3_Papers/03_scientific_data/02_data/08_phenotypes/phenotype/02_cognition_ses-01.tsv", sep = "\t", header = T, na.strings = "n/a")
ses01.cog <- left_join(qtab.data, ses01.cog, "participant_id")
ses01.cog <- ses01.cog %>% rename(age_months = ses01_age_months)

ses02.cog <- read.table("C:/Documents and Settings/uqlstri1/OneDrive - The University of Queensland/3_Papers/03_scientific_data/02_data/08_phenotypes/phenotype/02_cognition_ses-02.tsv", sep = "\t", header = T, na.strings = "n/a")
ses02.cog <- left_join(qtab.data, ses02.cog, "participant_id")
ses02.cog <- ses02.cog %>% rename(age_months = ses02_age_months)

ses01.cog <- ses01.cog %>% select(participant_id, ProcSpeed_raw, CrystallizedComposite_uncorr_stan, sex_female, age_months, indID, Pair, M)
ses02.cog <- ses02.cog %>% select(participant_id, ProcSpeed_raw, CrystallizedComposite_uncorr_stan, sex_female, age_months, indID, Pair, M)

ses01.cog$session <- 1
ses02.cog$session <- 2

ses01.02 <- rbind(ses01.cog, ses02.cog)

saveRDS(ses01.02, "QTAB_mixed_effects_datafile.RDS")
