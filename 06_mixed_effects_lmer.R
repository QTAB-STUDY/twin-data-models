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

ses01.02 <- readRDS("QTAB_mixed_effects_datafile.RDS")
ses01.02.unrelated <- ses01.02 %>% filter(indID==1)
ses01 <- ses01.02 %>% filter(session==1)

#### Model controlling for repeat observations & relatedness ####
mod1 <- lmer(ProcSpeed_raw ~ age_months + sex_female + (1 | participant_id) + (1 | Pair) + (1 | M), REML = FALSE, data = ses01.02)
summary(mod1)

#### Model controlling for repeat observations only (i.e. if using unrelated dataset) ####
mod2 <- lmer(ProcSpeed_raw ~ age_months + sex_female + (1 | participant_id), REML = FALSE, data = ses01.02.unrelated)
summary(mod2)

#### Model controlling for relatedness only (i.e. not using multi-session data ####
mod3 <- lmer(ProcSpeed_raw ~ age_months + sex_female + (1 | Pair) + (1 | M), REML = FALSE, data = ses01)
summary(mod3)

#### Compare models with different covariates ####
mod4 <- lmer(ProcSpeed_raw ~ age_months + (1 | participant_id) + (1 | Pair) + (1 | M), REML = FALSE, data = ses01.02)
mod5 <- lmer(ProcSpeed_raw ~ age_months + sex_female + (1 | participant_id) + (1 | Pair) + (1 | M), REML = FALSE, data = ses01.02)
compare_performance(mod4, mod5)
