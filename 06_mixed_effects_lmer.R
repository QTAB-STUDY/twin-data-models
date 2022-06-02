#### Notes ####
# Author: Lachlan T Strike
# Example code for running linear mixed effects models on QTAB data:
# Strike, L. T. et al. QTAB Non-Imaging Phenotypes. UQ eSpace https://doi.org/10.48610/ec1585d (2022).
# Strike, L. T. et al. Queensland Twin Adolescent Brain (QTAB). OpenNeuro (2022).
# Code is provided as a example/starting point (may contain errors). It should not be considered the best practice/approach!!!

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

#### Model controlling for repeat observations only (i.e., if using unrelated dataset) ####
mod2 <- lmer(ProcSpeed_raw ~ age_months + sex_female + (1 | participant_id), REML = FALSE, data = ses01.02.unrelated)
summary(mod2)

#### Model controlling for relatedness only (i.e., not using multi-session data ####
mod3 <- lmer(ProcSpeed_raw ~ age_months + sex_female + (1 | Pair) + (1 | M), REML = FALSE, data = ses01)
summary(mod3)

#### Compare models with different covariates ####
mod4 <- lmer(ProcSpeed_raw ~ age_months + (1 | participant_id) + (1 | Pair) + (1 | M), REML = FALSE, data = ses01.02)
mod5 <- lmer(ProcSpeed_raw ~ age_months + sex_female + (1 | participant_id) + (1 | Pair) + (1 | M), REML = FALSE, data = ses01.02)
mod6 <- lmer(ProcSpeed_raw ~ age_months * sex_female + (1 | participant_id) + (1 | Pair) + (1 | M), REML = FALSE, data = ses01.02)
compare_performance(mod4, mod5, mod6)

#### Plot ####
library(viridis)
library(scales)

theme_lachlan <- function() {
  theme_minimal(base_size = 12) +
    theme(
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      legend.position="bottom",
      legend.title = element_blank(),
      text=element_text(family="sans")
    )
}

show_col(viridis(n = 4, option = "A"))
myColors <- c("#721F81FF", "#F1605DFF")
names(myColors) <- levels(as.factor(ses01.02$sex))
colScale <- scale_colour_manual(name = "sex", values = myColors)

ses01.02 <- ses01.02 %>% mutate(sex = recode(sex_female, `0` = "Boys", `1` = "Girls"))

# Age on x-axis
g1 <- ggplot(data = ses01.02, aes(x = age_months, y = ProcSpeed_raw, color = sex))+
  geom_point(size = 1) + 
  xlab("Age (months)") +
  ylab("Processing speed (raw)") +
  ggtitle("QTAB Processing Speed") + 
  geom_line(size = 0.5, aes(group = participant_id)) +
  theme_lachlan() +
  colScale
  
g1 
