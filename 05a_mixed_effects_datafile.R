#### Notes ####
# Author: Lachlan T Strike
# Example code for running linear mixed effects models on QTAB data:
# Strike, L. T. et al. QTAB Non-Imaging Phenotypes. UQ eSpace https://doi.org/10.48610/e891597 (2022).
# Strike, L. T. et al. Queensland Twin Adolescent Brain (QTAB). OpenNeuro https://doi.org/10.18112/openneuro.ds004146.v1.0.0 (2022).
# Code is provided as a example/starting point (may contain errors). It should not be considered the best practice/approach!!!

#### Setup ####
library(tidyverse)
setwd("C:/GitHub/twin-data-models")
rm(list = ls())

# Create Pair & M variables for linear-mixed models
# See Visscher et al. (2004). The Use of Linear Mixed Models to Estimate Variance Components from Data on Twin Pairs by Maximum Likelihood
qtab.data <- read.table("non-imaging_phenotypes/participants_restricted.tsv", sep = "\t", header = T, na.strings = "n/a")
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

# Create Pair & M
qtab.data <- qtab.data %>% arrange(family_id)
qtab.fams <- unique(qtab.data$family_id)
qtab.fams.MZ <- unique(qtab.data[which(qtab.data$zyg<=2), "family_id"])
qtab.fams.DZ <- unique(qtab.data[which(qtab.data$zyg>=3), "family_id"])
qtab.data$Pair <- NA
qtab.data$M <- NA

length(qtab.fams) # 211 families
for(i in 1:211){
  famID <- qtab.fams[i]
  qtab.data[which(qtab.data$family_id==famID & !is.na(qtab.data$zyg)), "Pair"] <- paste0("Pair", i)
}

length(qtab.fams.MZ) #111 MZ families
for(i in 1:111){
  famID <- qtab.fams.MZ[i]
  qtab.data[which(qtab.data$family_id==famID & !is.na(qtab.data$zyg)), "M"] <- paste0("M", i)
}

length(qtab.fams.DZ) # 100 DZ families
M.DZ.start <- 112 # 1 + number of MZ families
M.DZ.end <- 311 # number of MZ families (111) + number of DZ twins (100*2). 111 + 200 = 311
qtab.data[which(qtab.data$family_id %in% qtab.fams.DZ  & !is.na(qtab.data$zyg)), "M"] <- paste0("M", M.DZ.start:M.DZ.end)

#### Save the Pair and M variables
qtab.data.lmm <- 
  qtab.data %>% select(participant_id, M, Pair) 
write.table(qtab.data.lmm, "qtab_Pair_M.txt", sep = '\t', quote = F, row.names = F)
