setwd("C:/GitHub/twin-data-models")
rm(list = ls())

# Load libraries
library(tidyverse)

#### Setup participants datafile ####
qtab.participants <- read_delim("non_imaging_phenotypes/participants_restricted.tsv", delim = "\t",  na = "n/a")
# Recode sex from M/F to sex_female 0/1
qtab.participants <- qtab.participants %>% mutate(sex_female = if_else(sex=="F", 1, 0))
head(qtab.participants %>% select(sex, sex_female))
# Create an individual ID based on birth order (i.e. Twin 1, Twin 2)
qtab.participants$indID <- qtab.participants$birth_order
# Maximise the twin dataset - recode pairs from triplets so that indID is 1 & 2 
# (with 1 & 2 reflecting the relative birth order)
qtab.participants <- qtab.participants %>%                               
  mutate(indID = replace(indID, participant_id == "8IxEjY", 2))
qtab.participants <- qtab.participants %>%                               
  mutate(indID = replace(indID, participant_id == "LIPs0p", 1))
qtab.participants <- qtab.participants %>%                               
  mutate(indID = replace(indID, participant_id == "DPfWmQ", 2))
qtab.participants <- qtab.participants %>%                               
  mutate(indID = replace(indID, participant_id == "KQb9f9", 2))
# Check that pairs now have indID 1 & 2
select(qtab.participants %>% filter(family_id=="9vHiUW"), c("participant_id", "birth_order", "indID"))
select(qtab.participants %>% filter(family_id=="CuoAOr"), c("participant_id", "birth_order", "indID"))
select(qtab.participants %>% filter(family_id=="siiJJD"), c("participant_id", "birth_order", "indID"))
# DZOS recode so females are 01, males are 02
qtab.participants <- qtab.participants %>%                               
  mutate(indID = replace(indID, zyg==5 & sex_female==1, 1))
qtab.participants <- qtab.participants %>%                               
  mutate(indID = replace(indID, zyg==5 & sex_female==0, 2))

#### Session 01 data ####
ses01.pub <- read_delim("non_imaging_phenotypes/01_puberty_ses-01.tsv", delim = "\t", na = "n/a")
ses01.cog <- read_delim("non_imaging_phenotypes/02_cognition_ses-01.tsv", delim = "\t", na = "n/a")
ses01.anxdep <- read_delim("non_imaging_phenotypes/03_anxiety_depression_ses-01.tsv", delim = "\t", na = "n/a")
ses01.emotsoc <- read_delim("non_imaging_phenotypes/04_emot_soc_behav_ses-01.tsv", delim = "\t", na = "n/a")
ses01.socsupp <- read_delim("non_imaging_phenotypes/05_social_support_family_functioning_ses-01.tsv", delim = "\t", na = "n/a")
ses01.stress <- read_delim("non_imaging_phenotypes/06_stress_ses-01.tsv", delim = "\t", na = "n/a")
ses01.sleep <- read_delim("non_imaging_phenotypes/07_sleep_physical_health_ses-01.tsv", delim = "\t", na = "n/a")
ses01.diet <- read_delim("non_imaging_phenotypes/09_dietary_behaviour_ses-01.tsv", delim = "\t", na = "n/a")
earlylifedemogs <- read_delim("non_imaging_phenotypes/08_early_life_family_demographics.tsv", delim = "\t", na = "n/a")

qtab.data.ses01 <- left_join(qtab.participants, ses01.pub, "participant_id")
qtab.data.ses01 <- left_join(qtab.data.ses01, ses01.cog, "participant_id")
qtab.data.ses01 <- left_join(qtab.data.ses01, ses01.anxdep, "participant_id")
qtab.data.ses01 <- left_join(qtab.data.ses01, ses01.emotsoc, "participant_id")
qtab.data.ses01 <- left_join(qtab.data.ses01, ses01.socsupp, "participant_id")
qtab.data.ses01 <- left_join(qtab.data.ses01, ses01.stress, "participant_id")
qtab.data.ses01 <- left_join(qtab.data.ses01, ses01.sleep, "participant_id")
qtab.data.ses01 <- left_join(qtab.data.ses01, ses01.diet, "participant_id")
qtab.data.ses01 <- left_join(qtab.data.ses01, earlylifedemogs, "participant_id")

# Scale (z-score) phenotypes of interest
# Useful for OpenMx model convergence & identifying outliers
ses.01.variable_list <- c("ProcSpeed_raw", "CJOLO", "FAS_60", "DSf", "SCAS_score", "SMFQ_score", "SPHERE_anxdep_score", "SPHERE_fat_score", "CRSQ_rum_score", "CASQ_comp_score", "MSPSS_score", "DLSS_score", "GBSvictim_score", "PDSS_score", 
                          "pPRMQ_score", "pSCAS_score", "pSMFQ_score", "pSDQ_total_score", "pEATQ_act_score", "kJwithDF")
for (i in ses.01.variable_list){
  qtab.data.ses01[, paste0(i, "Z")] <- as.numeric(scale(qtab.data.ses01[, i]))
}

##### Convert to familywise dataset (i.e. one family per row) ####
# _01 appended to twin 1, _02 appended to twin 2
famData.ses01 <- reshape(as.data.frame(qtab.data.ses01), timevar = "indID", 
                   idvar = "family_id", direction = "wide", sep = "_0")
# Create single zygosity variable based on twin 1
famData.ses01$zyg <- famData.ses01$zyg_01
# Check to make sure there are no families without zygosity
famData.ses01 %>% filter(is.na(zyg))

# Save familywise datafile
saveRDS(famData.ses01, "QTAB_familywise_ses01.RDS")

#### Session 02 data ####
ses02.pub <- read_delim("phenotype/01_puberty_ses-02.tsv", delim = "\t", na = "n/a")
ses02.cog <- read_delim("phenotype/02_cognition_ses-02.tsv", delim = "\t", na = "n/a")
ses02.anxdep <- read_delim("phenotype/03_anxiety_depression_ses-02.tsv", delim = "\t", na = "n/a")
ses02.emotsoc <- read_delim("phenotype/04_emot_soc_behav_ses-02.tsv", delim = "\t", na = "n/a")
ses02.socsupp <- read_delim("phenotype/05_social_support_family_functioning_ses-02.tsv", delim = "\t", na = "n/a")
ses02.stress <- read_delim("phenotype/06_stress_ses-02.tsv", delim = "\t", na = "n/a")
ses02.sleep <- read_delim("phenotype/07_sleep_physical_health_ses-02.tsv", delim = "\t", na = "n/a")
earlylifedemogs <- read_delim("phenotype/08_early_life_family_demographics.tsv", delim = "\t", na = "n/a")
covid19 <- read_delim("phenotype/10_covid-19.tsv", delim = "\t", na = "n/a")

qtab.data.ses02 <- left_join(qtab.participants, ses02.pub, "participant_id")
qtab.data.ses02 <- left_join(qtab.data.ses02, ses02.cog, "participant_id")
qtab.data.ses02 <- left_join(qtab.data.ses02, ses02.anxdep, "participant_id")
qtab.data.ses02 <- left_join(qtab.data.ses02, ses02.emotsoc, "participant_id")
qtab.data.ses02 <- left_join(qtab.data.ses02, ses02.socsupp, "participant_id")
qtab.data.ses02 <- left_join(qtab.data.ses02, ses02.stress, "participant_id")
qtab.data.ses02 <- left_join(qtab.data.ses02, ses02.sleep, "participant_id")
qtab.data.ses02 <- left_join(qtab.data.ses02, earlylifedemogs, "participant_id")
qtab.data.ses02 <- left_join(qtab.data.ses02, covid19, "participant_id")

# Scale (z-score) phenotypes of interest
# Useful for OpenMx model convergence & identifying outliers
ses.02.variable_list <- c("ProcSpeed_raw", "SCAS_score", "SMFQ_score", "SPHERE_anxdep_score", "SPHERE_fat_score", "CRSQ_rum_score", "CASQ_comp_score", "MSPSS_score", "DLSS_score", "GBSvictim_score", "PDSS_score", 
                          "pPRMQ_score", "pSCAS_score", "pSMFQ_score", "pSDQ_total_score", "pEATQ_act_score", "BRS_score")
for (i in ses.02.variable_list){
  qtab.data.ses02[, paste0(i, "Z")] <- as.numeric(scale(qtab.data.ses02[, i]))
}

##### Convert to familywise dataset (i.e. one family per row) ####
# _01 appended to twin 1, _02 appended to twin 2
famData.ses02 <- reshape(as.data.frame(qtab.data.ses02), timevar = "indID", 
                   idvar = "family_id", direction = "wide", sep = "_0")
# Create single zygosity variable based on twin 1
famData.ses02$zyg <- famData.ses02$zyg_01
# Check to make sure there are no families without zygosity
famData.ses02 %>% filter(is.na(zyg))

# Save familywise datafile
saveRDS(famData.ses02, "QTAB_familywise_ses02.RDS")
