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
qtab.data <- qtab.data %>% mutate(sex_female = if_else(sex=="F", 1, 0))
head(qtab.data %>% select(sex, sex_female))

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

# Z-score (Useful for OpenMx)
qtab.data$ProcSpeed_rawZ <- scale(qtab.data$ProcSpeed_raw)
qtab.data$CJOLOZ <- scale(qtab.data$CJOLO)

# Reshape from long to wide
famData <- reshape(qtab.data, timevar = "indID", idvar = "family_id", direction = "wide", sep = "_0") 
# Create single zygosity variable
famData$zyg <- famData$zyg_01
# Check to make sure there are no families without zygosity
famData %>% filter(is.na(zyg))

# Save familywise datafile
saveRDS(famData, "QTAB_familywise.RDS")
